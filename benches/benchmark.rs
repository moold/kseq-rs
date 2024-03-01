use criterion::{criterion_group, criterion_main, Criterion};
use kseq;
use needletail::parser::{FastaReader, FastqReader, FastxReader};
use std::{io::Cursor, iter};

fn simulate_fastq(total: usize) -> Vec<u8> {
    let mut data: Vec<u8> = vec![];
    let mut n = 0;
    let mut sum = 0;
    let mut seq_len = 100;
    loop {
        n += 1;
        data.push(b'@');
        data.extend(n.to_string().as_bytes());
        data.push(b'\n');
        if sum + seq_len > total {
            seq_len = total - sum;
        }
        data.extend(iter::repeat(b'A').take(seq_len));
        data.extend([b'\n', b'+', b'\n']);
        data.extend(iter::repeat(b'!').take(seq_len));
        data.push(b'\n');
        sum += seq_len;
        seq_len += 2;
        if sum >= total {
            break;
        }
    }
    // println!("{}", str::from_utf8(&data).unwrap());
    data
}

fn simulate_fasta(total: usize) -> Vec<u8> {
    let mut data: Vec<u8> = vec![];
    let mut n = 0;
    let mut sum = 0;
    let mut seq_len = 100;
    loop {
        n += 1;
        data.push(b'>');
        data.extend(n.to_string().as_bytes());
        data.push(b'\n');
        if sum + seq_len > total {
            seq_len = total - sum;
        }
        for _ in 0..seq_len / 100 {
            data.extend(iter::repeat(b'A').take(100));
            data.push(b'\n');
        }
        data.extend(iter::repeat(b'A').take(seq_len % 100));
        data.push(b'\n');
        sum += seq_len;
        seq_len += 2;
        if sum >= total {
            break;
        }
    }
    // println!("{}", str::from_utf8(&data).unwrap());
    data
}

fn bench_fasta_file(c: &mut Criterion) {
    let n_total = 1_000_000_000;
    let data = simulate_fasta(n_total);

    let mut group = c.benchmark_group("FASTA parsing(1GB)");
    group.sample_size(30);

    group.bench_function("kseq", |bench| {
        bench.iter(|| {
            let mut n_bases = 0;
            let mut records = kseq::parse_reader(Cursor::new(&data)).unwrap();
            while let Ok(Some(record)) = records.iter_record() {
                n_bases += record.seq().len() as u64;
            }
            assert_eq!(n_bases, n_total as u64);
        });
    });

    group.bench_function("needletail", |bench| {
        bench.iter(|| {
            let mut n_bases = 0;
            let mut records = FastaReader::new(Cursor::new(&data));
            while let Some(Ok(record)) = records.next() {
                n_bases += record.seq().len() as u64;
            }
            assert_eq!(n_bases, n_total as u64);
        });
    });

    group.finish();
}

fn bench_fastq_file(c: &mut Criterion) {
    let n_total = 1_000_000_000;
    let data = simulate_fastq(n_total);

    let mut group = c.benchmark_group("FASTQ parsing(1GB)");
    group.sample_size(30);

    group.bench_function("kseq", |bench| {
        bench.iter(|| {
            let mut n_bases = 0;
            let mut records = kseq::parse_reader(Cursor::new(&data)).unwrap();
            while let Ok(Some(record)) = records.iter_record() {
                n_bases += record.seq().len() as u64;
            }
            assert_eq!(n_bases, n_total as u64);
        });
    });

    group.bench_function("needletail", |bench| {
        bench.iter(|| {
            let mut n_bases = 0;
            let mut records = FastqReader::new(Cursor::new(&data));
            while let Some(Ok(record)) = records.next() {
                n_bases += record.seq().len() as u64;
            }
            assert_eq!(n_bases, n_total as u64);
        });
    });

    group.finish();
}

criterion_group!(io, bench_fastq_file, bench_fasta_file);
criterion_main!(io);
