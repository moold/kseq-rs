[![Crates.io](https://img.shields.io/crates/d/kseq?logo=rust)](https://github.com/moold/kseq/archive/refs/heads/main.zip)
[![Crates.io](https://img.shields.io/crates/v/kseq)](https://crates.io/crates/kseq)
[![docs.rs](https://img.shields.io/docsrs/kseq)](https://docs.rs/kseq/)
# kseq
`kseq` is a simple fasta/fastq (**fastx**) format parser library for [Rust](https://www.rust-lang.org/), its main function is to iterate over the records from fastx files (similar to [kseq](https://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser) in `C`). It uses shared buffer to read and store records, so the speed is very fast. It supports a plain or gz fastx file or [`io::stdin`](https://doc.rust-lang.org/std/io/fn.stdin.html), as well as a fofn (file-of-file-names) file, which contains multiple plain or gz fastx files (one per line).

Using `kseq` is very simple. Users only need to call `parse_path` to parse the path, and then use `iter_record` method to get each record.

- `parse_path` This function takes a path (`Option<String>` or `String` or `&str`) as input, a path can be a `fastx` file, `None` or `-` for [`io::stdin`](https://doc.rust-lang.org/std/io/fn.stdin.html), or a `fofn` file. It returns a `Result` type:
	- `Ok(T)`: A struct `T` with the `iter_record` method.
	- `Err(E)`: An error `E` including missing input, can't open or read, wrong fastx format or invalid path or file errors.

- `iter_record` This function can be called in a loop, it returns a `Result<Option<Record>>` type:
	- `Ok(Some(Record))`: A struct `Record` with methods:
		- `head -> &str`: get sequence id/identifier
		- `seq -> &str`:  get sequence
		- `des -> &str`:  get sequence description/comment
		- `sep -> &str`:  get separator
		- `qual -> &str`: get quality scores
		- `len -> usize`: get sequence length

		***Note:*** call `des`, `sep` and `qual` will return `""` if `Record` doesn't have these attributes.
	- `Ok(None)`: Stream has reached `EOF`.
	- `Err(ParseError)`: An error [`ParseError`](https://docs.rs/kseq/0.3.0/kseq/record/enum.ParseError.html) including `IO`, `TruncateFile`, `InvalidFasta` or `InvalidFastq` errors.

## Example
```no_run 
use std::env::args;
use kseq::parse_path;

fn main(){
	let path: Option<String> = args().nth(1);
	let mut records = parse_path(path).unwrap();
	while let Some(record) = records.iter_record().unwrap() {
		println!("head:{} des:{} seq:{} qual:{} len:{}", 
			record.head(), record.des(), record.seq(), 
			record.qual(), record.len());
	}
}
```

## Installation
```text 
[dependencies]
kseq = "0.3"
```
