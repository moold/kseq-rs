// Note: kseq is inspired by fastq-rs and kseq in C
//! # kseq
//! `kseq` is a simple fasta/fastq (**fastx**) format parser library for [Rust](https://www.rust-lang.org/), its main function is to iterate over the records from fastx files (similar to [kseq](https://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser) in `C`). It uses shared buffer to read and store records, so the speed is very fast. It supports a plain or gz fastx file or [`io::stdin`](https://doc.rust-lang.org/std/io/fn.stdin.html), as well as a fofn (file-of-file-names) file, which contains multiple plain or gz fastx files (one per line).
//! 
//! Using `kseq` is very simple. Users only need to call `parse_path` to parse the path, and then use `iter_record` method to get each record.
//! 
//! - `parse_path` This function takes a path (`Option<String>`) as input, a path can be a fastx file, `None` or `-` for [`io::stdin`](https://doc.rust-lang.org/std/io/fn.stdin.html), or a fofn file. It returns a [`Result`](https://doc.rust-lang.org/std/result/) type:
//! 	- `Ok(T)`: An struct `T` with the `iter_record` method.
//! 	- `Err(E)`: An error `E` including can't open or read, wrong fastx format or invalid path or file errors.
//! 
//! - `iter_record` This function can be called in a loop, it return an [`Option`](https://doc.rust-lang.org/std/option/index.html) type:
//! 	- `Some(Record)`: An struct `Record` with methods:
//! 		- `head -> &str`: get sequence id/identifier
//! 		- `seq -> &str`:  get sequence
//! 		- `des -> &str`:  get sequence description/comment
//! 		- `sep -> &str`:  get separator
//! 		- `qual -> &str`: get quality scores
//! 		- `len -> usize`:  get sequence length
//! 
//! 		***Note:*** call `des`, `sep` and `qual` will return `""` if `Record` doesn't have these attributes. `head`, `seq`, `des`, `sep` and `qual` use unsafe code [`str::from_utf8_unchecked`](https://doc.rust-lang.org/std/str/fn.from_utf8_unchecked.html).
//! 	- `None`: Stream has reached `EOF`.
//! 
//! # Examples
//! ```
//! use std::env::args;
//! use kseq::path::parse_path;
//! 
//! fn main(){
//! 	let path: Option<String> = args().nth(1);
//! 	let mut records = parse_path(path).unwrap();
//! 	while let Some(record) = records.iter_record() {
//! 		println!("head:{} des:{} seq:{} qual:{} len:{}", 
//! 			record.head(), record.des(), record.seq(), record.qual(), record.len());
//! 	}
//! }
//! ```

pub mod path;
pub mod record;
