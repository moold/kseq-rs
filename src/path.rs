use std::io::{Result, Read, BufRead, BufReader, Error, ErrorKind, Cursor};
use flate2::read::MultiGzDecoder;

use crate::record::{Reader, Paths};

pub fn parse_path(path: Option<String>) -> Result<Reader>{
    let mut reader: Box<dyn BufRead> = match path.as_ref().map(String::as_str) {
        None | Some("-") => {
            Box::new(BufReader::with_capacity(65536, std::io::stdin()))
        },
        Some(path) => {
            Box::new(BufReader::with_capacity(65536, std::fs::File::open(path)?))
        }
    };
    let mut magic_bytes = [0u8; 4];
    reader.read_exact(&mut magic_bytes)?;
    reader = Box::new(Cursor::new(magic_bytes.to_vec()).chain(reader));
   if &magic_bytes[..2] == b"\x1f\x8b" {
        reader = Box::new(BufReader::with_capacity(65536, MultiGzDecoder::new(reader)));
        // magic_bytes.iter_mut().for_each(|m| *m = 0);
        reader.read_exact(&mut magic_bytes)?;
        reader = Box::new(Cursor::new(magic_bytes.to_vec()).chain(reader));
    }

    match magic_bytes[0] {
        b'@' => Ok(Reader::new(reader, true)),
        b'>' => Ok(Reader::new(reader, false)),
        _ => Err(Error::new(ErrorKind::InvalidData, "Not a gzip or plain fastq/fasta file"))
    }
}

pub fn parse_paths(path: Option<String>) -> Result<Paths>{
    let reader: Box<dyn BufRead> = match path.as_ref().map(String::as_str) {
        None | Some("-") => {
            Box::new(BufReader::with_capacity(65536, std::io::stdin()))
        },
        Some(path) => {
            Box::new(BufReader::with_capacity(65536, std::fs::File::open(path)?))
        }
    };
    
    let mut result = Paths::new();
    for line in reader.lines() {
        let path = line.unwrap();
        if path.trim().starts_with('#') || path.trim().is_empty(){
            continue;
        }
        result.paths.push(parse_path(Some(path)));
    }
    Ok(result)
}