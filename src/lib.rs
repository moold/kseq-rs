#![doc = include_str!("../README.md")]
// Note: kseq is inspired by fastq-rs and kseq in C

use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader, Cursor, Error, ErrorKind, Read, Result};
use std::path::Path;

pub mod record;
use record::{Fastx, Reader, Readers, Result as ParseResult};

// A wrapper that allows parse_path to accept Option<String>, String and &str as parameter
pub struct Wrapper(String);

impl From<Option<String>> for Wrapper {
    fn from(path: Option<String>) -> Self {
        Wrapper(path.unwrap_or_else(|| "-".into()))
    }
}

impl From<String> for Wrapper {
    fn from(path: String) -> Self {
        Wrapper(path)
    }
}

impl From<&str> for Wrapper {
    fn from(path: &str) -> Self {
        Wrapper(path.into())
    }
}

impl From<&&str> for Wrapper {
    fn from(path: &&str) -> Self {
        Wrapper((*path).into())
    }
}

/// Reader for a single path or Readers for multiple paths
pub enum Paths {
    Reader(Reader),
    Readers(Readers),
}

impl Paths {
    /// iterate a fatsx record for a Reader or Readers
    pub fn iter_record(&mut self) -> ParseResult<Option<Fastx>> {
        match self {
            Paths::Reader(t) => t.iter_record_check(),
            Paths::Readers(t) => t.iter_record(),
        }
    }
}

/// parse path to a Reader or Readers
pub fn parse_path<T>(path: T) -> Result<Paths>
where
    T: Into<Wrapper>,
{
    let path = path.into().0;
    let mut reader: Box<dyn BufRead> = if path == "-" {
        if atty::is(atty::Stream::Stdin) {
            return Err(Error::new(ErrorKind::InvalidInput, "Missing input"));
        }
        Box::new(BufReader::with_capacity(65536, std::io::stdin()))
    } else {
        Box::new(BufReader::with_capacity(65536, std::fs::File::open(&path)?))
    };

    let mut format_bytes = [0u8; 4];
    reader.read_exact(&mut format_bytes)?;
    reader = Box::new(Cursor::new(format_bytes.to_vec()).chain(reader));
    if &format_bytes[..2] == b"\x1f\x8b" {
        // for gz foramt
        reader = Box::new(BufReader::with_capacity(65536, MultiGzDecoder::new(reader)));
        format_bytes.iter_mut().for_each(|m| *m = 0);
        reader.read_exact(&mut format_bytes)?;
        reader = Box::new(Cursor::new(format_bytes.to_vec()).chain(reader));
    }

    match format_bytes[0] {
        b'@' | b'>' => Ok(Paths::Reader(Reader::new(reader))),
        _ => {
            // for a fofn file
            let mut paths = Readers::new();
            let _path = if path == "-" { "" } else { &path };
            let parent = Path::new(&_path).parent().unwrap_or_else(|| Path::new(""));

            for _line in reader.lines().map(|l| l.unwrap()) {
                let line = _line.trim();
                if line.starts_with('#') || line.is_empty() {
                    continue;
                }
                let _path = parent.join(line); // convert to an absolute path
                if _path.exists() {
                    match parse_path(Some(_path.to_str().unwrap().to_string()))? {
                        Paths::Reader(reader) => paths.readers.push(reader),
                        _ => unreachable!(),
                    }
                } else {
                    return Err(Error::new(
                        ErrorKind::InvalidData,
                        format!("{:?} is not a valid fastq/fasta/fofn file", _path),
                    ));
                }
            }
            Ok(Paths::Readers(paths))
        }
    }
}
