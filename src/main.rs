use rand::{thread_rng, Rng};

use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::path::Path;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

pub fn read_fasta(f: &File) -> io::Result<Vec<u8>> {
    let reader = BufReader::new(f);
    let mut genome: Vec<u8> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if !line.contains('>') {
            for nucleotide in line.chars() {
                genome.push(nucleotide as u8);
            }
        }
    }

    Ok(genome)
}

fn random_subsequence(genome: &[u8], length: usize) -> &[u8] {
    let low = 0;
    let high = genome.len() - length;
    let start = thread_rng().gen_range(low, high);
    &genome[start..start + length]
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

fn blast_iter(reads: &[&[u8]]) -> Duration {
    let query = reads
        .iter()
        .enumerate()
        .map(|pair| {
            let (i, read) = pair;
            format!(">seq_{}\n{}", i, std::str::from_utf8(read).unwrap())
        })
        .collect::<Vec<String>>()
        .join("\n");
    let start = Instant::now();
    let mut cmd = Command::new("blastn")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .stdin(Stdio::piped())
        .arg("-db")
        .arg("/tmp/dm6.db")
        .arg("-task")
        .arg("blastn-short")
        .arg("-qcov_hsp_perc")
        .arg("100")
        .arg("-outfmt")
        .arg("10 sstart")
        .spawn()
        .expect("blast iteration failed");
    if let Some(ref mut fd) = cmd.stdin {
        writeln!(fd, "{}", query).unwrap();
    } else {
        panic!("No stdin");
    }
    cmd.wait().unwrap();
    Instant::now().duration_since(start)
}

fn measure_iter<F>(
    genome: &[u8],
    iter: F,
    seq_len: usize,
    seq_count: usize,
) -> Result<Duration, std::io::Error>
where
    F: Fn(&[&[u8]]) -> Duration,
{
    Ok(iter(
        &(0..seq_count)
            .map(|_| random_subsequence(&genome, seq_len))
            .collect::<Vec<&[u8]>>(),
    ))
}

fn main() -> Result<(), std::io::Error> {
    let f = File::open("/home/jakob/University/BIOL 396/dm6/dm6.fa").unwrap();
    let genome = read_fasta(&f).unwrap();

    for _ in 0..100 {
        println!("{:?}", measure_iter(&genome, blast_iter, 10, 1)?);
    }

    Ok(())
}
