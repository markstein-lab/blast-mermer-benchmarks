use rand::{thread_rng, Rng};

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

type Chromosome = (String, usize);
type Genome = (Vec<Chromosome>, Vec<u8>);
type IterResult = (Vec<String>, Duration, Duration);

/// Read the reference sequence at `f`.
///
/// Returns a Vec of chromosome markers, and a flat Vec of nucleotides.
pub fn read_fasta(f: &File) -> io::Result<Genome> {
    let mut genome: Vec<u8> = Vec::new();
    let mut chromosomes: Vec<(String, usize)> = Vec::new();

    let reader = BufReader::new(f);

    for line in reader.lines() {
        let line = line?;
        if line.contains('>') {
            chromosomes.push((
                String::from(if let Some(n) = line.find(' ') {
                    &line[1..n]
                } else {
                    &line[1..]
                }),
                genome.len(),
            ));
        } else {
            for nucleotide in line.chars() {
                genome.push(nucleotide as u8);
            }
        }
    }

    Ok((chromosomes, genome))
}

/// Return a random `length`-mer in `genome`.
fn random_subsequence(genome: &[u8], length: usize) -> &[u8] {
    let low = 0;
    let high = genome.len() - length;
    let start = thread_rng().gen_range(low, high);
    &genome[start..start + length]
}

fn measure_iter<F>(
    genome: &Genome,
    iter: F,
    seq_len: usize,
    seq_count: usize,
) -> io::Result<IterResult>
where
    F: Fn(&[&[u8]]) -> IterResult,
{
    let (chromosomes, genome) = genome;

    let sequences = &(0..seq_count)
        .map(|_| random_subsequence(&genome, seq_len))
        .collect::<Vec<&[u8]>>();
    let results = iter(sequences);

    // Validate results against known values.

    let actual_hits = sequences
        .iter()
        .flat_map(|seq| find_hits(seq, &genome))
        .collect::<Vec<usize>>();

    let mut reported_hits = results
        .0
        .iter()
        .filter(|line| line.len() > 0)
        .map(|line| {
            let components: Vec<&str> = line.split(",").collect();
            let offset = chromosome_offset(components[0], &chromosomes);
            offset.unwrap() + components[1].parse::<usize>().unwrap() - 1
        })
        .collect::<Vec<usize>>();
    reported_hits.sort();

    assert!(actual_hits.len() == reported_hits.len());
    reported_hits
        .iter()
        .zip(actual_hits)
        .for_each(|(a, b)| assert!(*a == b));

    Ok(results)
}

/// Naïve motif search algorithm.
///
/// Return a Vec of offsets at which `motif` can be found in `genome`.
fn find_hits(motif: &[u8], genome: &[u8]) -> Vec<usize> {
    let mut result = Vec::new();

    for i in 0..genome.len() - motif.len() {
        if *motif == genome[i..i + motif.len()] {
            result.push(i);
        }
    }

    result
}

/// Return the offset at which `chromosome` begins.
///
/// Offset in a 0-indexed system, where 0 is the first nucleotide in the
/// reference sequence.
fn chromosome_offset(chromosome: &str, chromosomes: &[Chromosome]) -> Option<usize> {
    for (name, position) in chromosomes {
        if name == chromosome {
            return Some(*position);
        }
    }
    return None;
}

fn mermer_iter(reads: &[&[u8]]) -> Duration {
    let query = format!(
        "\n{}\n\n1\n\n1000\n",
        reads.iter().fold("".to_string(), |a, b| a
            + &std::str::from_utf8(b).unwrap()
            + "\n")
    );
    let start = Instant::now();
    let mut cmd = Command::new("scanner")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .stdin(Stdio::piped())
        .spawn()
        .expect("mermer iteration failed");
    if let Some(ref mut fd) = cmd.stdin {
        writeln!(fd, "{}", query).unwrap();
    } else {
        panic!("No stdin");
    }
    cmd.wait().unwrap();
    Instant::now().duration_since(start)
}

fn blast_iter(reads: &[&[u8]]) -> IterResult {
    let query = reads
        .iter()
        .enumerate()
        .map(|pair| {
            let (i, read) = pair;
            format!(">seq_{}\n{}", i, std::str::from_utf8(read).unwrap())
        })
        .collect::<Vec<String>>()
        .join("\n");
    let word_size = reads[0].len() - 1;
    let start = Instant::now();
    let mut cmd = Command::new("blastn")
        .stderr(Stdio::null())
        .stdout(Stdio::piped())
        .stdin(Stdio::piped())
        .arg("-db")
        .arg("dm6.db")
        .arg("-task")
        .arg("blastn-short")
        .arg("-qcov_hsp_perc")
        .arg("100")
        .arg("-outfmt")
        .arg("10 sseqid sstart")
        .arg("-evalue")
        .arg("0.001")
        .arg("-word_size")
        .arg(word_size.to_string())
        .arg("-dust")
        .arg("no")
        .arg("-strand")
        .arg("plus")
        .spawn()
        .expect("blast iteration failed");

    if let Some(ref mut fd) = cmd.stdin {
        writeln!(fd, "{}", query).unwrap();
    } else {
        panic!("No stdin");
    }

    cmd.wait().unwrap();
    let stop = Instant::now();

    let mut output = Vec::new();
    cmd.stdout.take().unwrap().read_to_end(&mut output).unwrap();

    let output = String::from_utf8_lossy(&output[..]);
    let lines: Vec<&str> = output.split('\n').collect();

    let time_to_search = lines[0][20..].parse::<u64>().unwrap();

    (
        lines[1..].iter().map(|s| s.trim().to_string()).collect(),
        stop.duration_since(start),
        Duration::from_millis(time_to_search),
    )
}

fn main() -> Result<(), std::io::Error> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() != 4 {
        eprintln!("usage: {} [program] [query_length] [query_count]", args[0]);
        std::process::exit(1);
    }

    let f = File::open("dm6.fa").unwrap();
    let genome = read_fasta(&f).unwrap();

    let query_length = args[2].parse::<usize>().unwrap();
    let query_count = args[3].parse::<usize>().unwrap();
    let program = match args[1].as_ref() {
        "blast" => blast_iter,
        // "mermer" => mermer_iter,
        _ => panic!("Invalid program name"),
    };

    for _ in 0..query_count {
        let (_, canonical, reported) = measure_iter(&genome, program, query_length, 1)?;
        println!("{},{}", canonical.as_millis(), reported.as_millis());
    }

    Ok(())
}
