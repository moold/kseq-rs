use std::{
    error, fmt,
    io::{self, ErrorKind},
    str,
};

pub type Result<T> = std::result::Result<T, ParseError>;

/// The type of error that returned during parsing fastx files
#[derive(Debug)]
pub enum ParseError {
    /// IO error
    Io(io::Error),
    /// A truncated record was found
    TruncateFile(String),
    /// Not a valid fastx record, the record doesn't start with `>` and '@'
    InvalidFastx(String),
    /// Not a valid fasta record, the record starts with `>` but the sequence length is 0
    InvalidFasta(String),
    /// Not a valid fastq record, the record start with `@` but the sequence and quality lengths are not equal or 0
    InvalidFastq(String),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ParseError::Io(err) => write!(f, "IO error: {}", err),
            ParseError::TruncateFile(record) => {
                write!(f, "Truncate file, problematic record: {}", record)
            }
            ParseError::InvalidFastx(record) => {
                write!(f, "Not a valid fastx record: {}", record)
            }
            ParseError::InvalidFasta(record) => {
                write!(f, "Not a valid fasta record: {}", record)
            }
            ParseError::InvalidFastq(record) => {
                write!(f, "Not a valid fastq record: {}", record)
            }
        }
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError::Io(err)
    }
}

impl error::Error for ParseError {}

/// a structure representing the sequence in a fastx file
pub struct Fastx<'a> {
    _head: usize,
    _des: usize,
    _seq: usize,
    _sep: usize,
    _qual: usize,
    _data: &'a Vec<u8>,
}

impl Fastx<'_> {
    /// get sequence id/identifier
    #[inline]
    pub fn head(&self) -> &str {
        unsafe { str::from_utf8_unchecked(&self._data[1..self._head]) }
    }

    /// get sequence
    #[inline]
    pub fn seq(&self) -> &str {
        unsafe { str::from_utf8_unchecked(&self._data[self._des..self._seq]) }
    }

    /// get sequence description/comment
    #[inline]
    pub fn des(&self) -> &str {
        if self._head < self._des {
            unsafe { str::from_utf8_unchecked(&self._data[self._head..self._des]) }
        } else {
            ""
        }
    }

    /// get separator
    #[inline]
    pub fn sep(&self) -> &str {
        if self._seq < self._sep {
            unsafe { str::from_utf8_unchecked(&self._data[self._seq..self._sep]) }
        } else {
            ""
        }
    }

    /// get quality scores
    #[inline]
    pub fn qual(&self) -> &str {
        if self._sep < self._qual {
            unsafe { str::from_utf8_unchecked(&self._data[self._sep..self._qual]) }
        } else {
            ""
        }
    }

    /// get sequence length
    #[inline]
    pub fn len(&self) -> usize {
        self._seq - self._des
    }

    /// check whether a fastx record is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// check whether a fastx record is a fasta record
    pub fn is_fasta(&self) -> bool {
        (!self._data.is_empty()) && self._data[0] == b'>'
    }

    /// check whether a fastx record is a fastq record
    pub fn is_fastq(&self) -> bool {
        (!self._data.is_empty()) && self._data[0] == b'@'
    }

    /// check a fastq record is valid
    fn validate_fastq(&self) -> bool {
        self.is_fastq() && !self.is_empty() && self._seq - self._des == self._qual - self._sep && self._head > 1
    }

    /// check a fasta record is valid
    fn validate_fasta(&self) -> bool {
        self.is_fasta() && !self.is_empty() && self._head > 1
    }
}

/// a reader with shared buffer
pub struct Reader<'a> {
    reader: Box<dyn io::BufRead + 'a>,
    data: Vec<u8>,
}

impl<'a> Reader<'a> {
    // Create a new Reader
    pub(crate) fn new(r: Box<dyn io::BufRead + 'a>) -> Self {
        Reader {
            reader: r,
            data: Vec::with_capacity(1024),
        }
    }

    // Check if this reader has any data left to be read.
    fn has_data_left(&mut self) -> Result<bool> {
        loop{
            let available = self.reader.fill_buf().map_err(ParseError::Io)?;
            if available.iter().any(|&x| !char::is_whitespace(x as char)){
                return Ok(true);
            }else if available.is_empty() {
                return Ok(false);
            }
            let len = available.len();
            self.reader.consume(len);
        }
    }

    // Return the next byte of the internal buffer
    #[allow(dead_code)]
    fn next_byte(&mut self) -> Result<Option<u8>> {
        loop {
            match self.reader.fill_buf() {
                Ok(n) => return Ok(n.first().copied()),
                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                Err(e) => return Err(ParseError::Io(e)),
            };
        }
    }

    // Read all non-newline bytes into data until the newline byte or EOF is reached,
    // the newline byte (if found) will not be appended to data.
    fn read_line(&mut self, skip_blank_line: bool) -> Result<usize> {
        let delim = b'\n';
        loop {
            let mut n = self.reader.read_until(delim, &mut self.data)?;
            // reached EOF
            if n == 0 {
                return Ok(n);
            }

            if self.data.last() == Some(&delim) {
                self.data.pop();
                n -= 1;
            }
            if n != 0 || !skip_blank_line {
                return Ok(n);
            }
        }
    }

    // Read all non-newline bytes into data until the delimiter byte or EOF is reached,
    // the delimiter (if found) will not be appended to data.
    fn read_until(&mut self, delim: u8) -> Result<usize> {
        let mut read = 0;
        loop {
            let (done, used) = {
                let available = match self.reader.fill_buf() {
                    Ok(n) => n,
                    Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                    Err(e) => return Err(ParseError::Io(e)),
                };
                let mut s = 0;
                let mut mch = memchr::memchr2_iter(delim, b'\n', available);
                loop {
                    match mch.next() {
                        Some(i) => {
                            self.data.extend_from_slice(&available[s..i]);
                            read += i - s;
                            s = i + 1;
                            if available[i] == delim {
                                break (true, i);
                            }
                        }
                        None => {
                            self.data.extend_from_slice(&available[s..]);
                            read += available.len() - s;
                            break (false, available.len());
                        }
                    }
                }
            };
            self.reader.consume(used);
            if done || used == 0 {
                return Ok(read);
            }
        }
    }

    // Read the exact number of non-newline bytes into data.
    fn read_exact(&mut self, len: usize) -> Result<usize> {
        let mut read = 0;
        loop {
            let (done, used) = {
                let available = match self.reader.fill_buf() {
                    Ok(n) => n,
                    Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                    Err(e) => return Err(ParseError::Io(e)),
                };
                let mut s = 0;
                let mut mch = memchr::memchr_iter(b'\n', available);
                loop {
                    match mch.next() {
                        Some(i) => {
                            if read + i - s >= len {
                                let e = len - read + s;
                                self.data.extend_from_slice(&available[s..e]);
                                read = len;
                                break (true, e);
                            } else {
                                self.data.extend_from_slice(&available[s..i]);
                                read += i - s;
                            }
                            s = i + 1;
                        }
                        None => {
                            if available.len() - s + read >= len {
                                let e = len - read + s;
                                self.data.extend_from_slice(&available[s..e]);
                                read = len;
                                break (true, e);
                            } else {
                                self.data.extend_from_slice(&available[s..]);
                                read += available.len() - s;
                                break (false, available.len());
                            }
                        }
                    }
                }
            };
            self.reader.consume(used);
            if done || used == 0 {
                return Ok(read);
            }
        }
    }

    /// iterate over a record from this Reader
    pub fn iter_record(&mut self) -> Result<Option<Fastx>> {
        // clean the last record
        self.data.clear();
        // read sequence head
        let des = self.read_line(true)?;
        if des == 0 {
            // reach the EOF
            return Ok(None);
        } else if self.data[0] != b'>' && self.data[0] != b'@' {
            // safely unwrap
            return Err(ParseError::InvalidFastx(
                String::from_utf8(self.data.to_owned()).unwrap(),
            ));
        }

        let head = self
            .data
            .iter()
            .position(|&x| char::is_whitespace(x as char))
            .map_or_else(|| des, |x| x);
        let mut seq = des;
        let mut sep = seq;
        let mut qual = sep;

        let is_fasta = self.data[0] == b'>';
        if is_fasta {
            seq += self.read_until(b'>')?;
        } else {
            seq += self.read_until(b'+')?;
            sep = seq + self.read_line(true)?;
            qual = sep + self.read_exact(seq - des)?;
        }

        if !self.has_data_left()? && (head == 1 || seq == des || (!is_fasta && (sep == seq || qual == sep))){
            // safely unwrap
            return Err(ParseError::TruncateFile(
                String::from_utf8(self.data.to_owned()).unwrap(),
            ));
        }
        // println!("head:{head} des {des} seq {seq} sep {sep} qual {qual}");
        let fastx = Fastx {
            _head: head,
            _des: des,
            _seq: seq,
            _sep: sep,
            _qual: qual,
            _data: &self.data,
        };

        if is_fasta && !fastx.validate_fasta() {
            return Err(ParseError::InvalidFasta(fastx.head().to_string()));
        } else if !(is_fasta || fastx.validate_fastq()) {
            return Err(ParseError::InvalidFastq(fastx.head().to_string()));
        }
        Ok(Some(fastx))
    }
}

/// multiple readers for a fofn file
pub struct Readers<'a> {
    index: usize,
    pub(crate) readers: Vec<Reader<'a>>,
}

impl<'a> Default for Readers<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> Readers<'a> {
    /// create a new Readers
    pub(crate) fn new() -> Self {
        Readers {
            index: 0,
            readers: Vec::new(),
        }
    }

    /// iterate over a record from this Readers
    pub(crate) fn iter_record(&mut self) -> Result<Option<Fastx>> {
        for idx in self.index..self.readers.len() {
            if self.readers[idx].has_data_left()? {
                return self.readers[idx].iter_record();
            }
            self.index += 1;
        }
        Ok(None)
    }
}
