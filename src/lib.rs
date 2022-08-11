#![doc = include_str!("../README.md")]
// Note: kseq is inspired by fastq-rs and kseq in C

use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{stdin, BufRead, BufReader, Cursor, Error, ErrorKind, Read, Result},
    path::Path,
};

pub mod record;
use record::{Fastx, Reader, Readers, Result as ParseResult};

/// Reader for a single path or Readers for multiple paths
pub enum Paths {
    Reader(Reader),
    Readers(Readers),
}

impl Paths {
    // parse a reader to a Reader or Readers
    fn new(mut reader: Box<dyn BufRead>, path: &Path) -> Result<Self> {
        let mut format_bytes = [0u8; 4];
        reader.read_exact(&mut format_bytes)?;
        reader = Box::new(Cursor::new(format_bytes.to_vec()).chain(reader));
        if &format_bytes[..2] == b"\x1f\x8b" {
            // for gz format
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
                let parent = path.parent().unwrap_or_else(|| Path::new(""));

                for _line in reader.lines() {
                    let _line = _line?;
                    let line = _line.trim();
                    if line.starts_with('#') || line.is_empty() {
                        continue;
                    }
                    let path = parent.join(line); // convert to an absolute path
                    if path.exists() {
                        match parse_path(path)? {
                            Paths::Reader(reader) => paths.readers.push(reader),
                            Paths::Readers(readers) => paths.readers.extend(readers.readers),
                        }
                    } else {
                        return Err(Error::new(
                            ErrorKind::InvalidData,
                            format!("{:?} is not a valid fastq/fasta/fofn file", path),
                        ));
                    }
                }
                Ok(Paths::Readers(paths))
            }
        }
    }

    /// iterate a fatsx record for a Reader or Readers
    pub fn iter_record(&mut self) -> ParseResult<Option<Fastx>> {
        match self {
            Paths::Reader(t) => t.iter_record_check(),
            Paths::Readers(t) => t.iter_record(),
        }
    }
}

/// parse path to a Reader or Readers
pub fn parse_path<P: AsRef<Path>>(path: P) -> Result<Paths> {
    let path = path.as_ref();
    let reader: Box<dyn BufRead> = if path == Path::new("-") {
        if atty::is(atty::Stream::Stdin) {
            return Err(Error::new(ErrorKind::InvalidInput, "Missing input"));
        }
        Box::new(BufReader::with_capacity(65536, stdin()))
    } else {
        Box::new(BufReader::with_capacity(65536, File::open(path)?))
    };
    Paths::new(reader, path)
}

/// parse reader to a Reader or Readers
pub fn parse_reader<R: Read + 'static>(reader: R) -> Result<Paths> {
    Paths::new(
        Box::new(BufReader::with_capacity(65536, reader)),
        Path::new(""),
    )
}
