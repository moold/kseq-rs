use std::{error, fmt, io, str};

pub type Result<T> = std::result::Result<T, ParseError>;

/// The type of error that returned during parsing fastx files
#[derive(Debug)]
pub enum ParseError {
	/// IO error
	Io(io::Error),
	/// A truncated record was found
	TruncateFile(String),
	/// Not a validate fasta record, a record didn't start with `>` or sequence length is 0
	InvalidFasta(String),
	/// Not a validate fastq record, a record didn't start with `@` or sequence and quality lengths are not equal or 0 
	InvalidFastq(String)
}

impl fmt::Display for ParseError {
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		match self {
			ParseError::Io(err) => write!(f, "IO error: {}", err),
			ParseError::TruncateFile(record) => write!(f, "Truncate file, problematic record: {}", record),
			ParseError::InvalidFasta(record) => write!(f, "Not a validate fasta record: {}", record),
			ParseError::InvalidFastq(record) => write!(f, "Not a validate fastq record: {}", record),
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
	_data: &'a Vec<u8>
}

impl Fastx<'_> {
	/// get sequence id/identifier 
	#[inline]
	pub fn head(&self) -> &str {
		unsafe {
			str::from_utf8_unchecked(&self._data[1 .. self._head - 1])
		}
	}

	/// get sequence
	#[inline]
	pub fn seq(&self) -> &str {
		unsafe {
			str::from_utf8_unchecked(&self._data[self._des .. self._seq - 1])
		}
	}

	/// get sequence description/comment
	#[inline]
	pub fn des(&self) -> &str {
		if self._des != self._head {
			unsafe {
				str::from_utf8_unchecked(&self._data[self._head .. self._des - 1])
			}
		}else {
			""
		}
	}

	/// get separator
	#[inline]
	pub fn sep(&self) -> &str {
		if self._seq != self._sep{
			unsafe {
				str::from_utf8_unchecked(&self._data[self._seq .. self._sep - 1])
			}
		}else{
			""
		}
	}

	/// get quality scores
	#[inline]
	pub fn qual(&self) -> &str {
		if self._sep != self._qual {
			unsafe {
				str::from_utf8_unchecked(&self._data[self._sep .. self._qual - 1])
			}
		}else{
			""
		}
	}
	
	/// get sequence length
	#[inline]
	pub fn len(&self) -> usize {
		self._seq - self._des - 1
	}
	
	/// check a fastq record is valid
	fn validate_fastq(&self)-> bool {
		self._seq - self._des == self._qual - self._sep && self._data[0] == b'@' && self.len() != 0
	}

	/// check a fasta record is valid
	fn validate_fasta(&self)-> bool {
		self._data[0] == b'>' && self.len() != 0
	}
}

/// a Reader with shared buffer
pub struct Reader {
	reader: Box<dyn io::BufRead>,
	data: Vec<u8>,
}

impl Reader {
	/// create a new Reader
	pub fn new(r: Box<dyn io::BufRead>) -> Self {
		Reader {
			reader: r,
			data: Vec::with_capacity(1024),
		}
	}

	/// check whether this Reader has reached EOF
	fn reach_eof(&mut self) -> Result<bool> {
		loop {
			self.data.clear();
			let head = self.reader.read_until(b'\n', &mut self.data)?;
			if head == 0 {return Ok(true);}
			// skip blank line, assume a valid line doesn't start with a space
			if !char::is_whitespace(self.data[0] as char) {break;}
		}
		Ok(false)
	}

	/// iterate over a record from this Reader
	pub fn iter_record(&mut self) -> Result<Option<Fastx>> {

		let head = self.data.iter().position(|&x| char::is_whitespace(x as char)).unwrap() + 1;
		let des = self.data.len();
		let mut seq = des + self.reader.read_until(b'\n', &mut self.data)?;
		let mut sep = seq;
		let mut qual = seq;
		let mut is_fasta = true;
		if seq != des {
			if self.data[0] == b'@' {
				is_fasta = false;
				let mut seq_line = 1;
				loop {
					let _data = self.reader.fill_buf()?; // for multiline fastq
					if !_data.is_empty() && _data[0] != b'+' {
						self.data.pop();
						seq += self.reader.read_until(b'\n', &mut self.data)? - 1;
						seq_line += 1;
					}else {
						break;
					}
				}
				sep = seq + self.reader.read_until(b'\n', &mut self.data)?;
				qual = sep + self.reader.read_until(b'\n', &mut self.data)?;
				if qual != sep {
					while seq_line > 1 {
						self.data.pop();
						let r = self.reader.read_until(b'\n', &mut self.data)?;
						if r == 0 {
							qual -= 1;
							break;
						}
						qual += r - 1;
						seq_line -= 1;
					}
					if self.data[qual - 1] != b'\n' { // for files not ending with '\n'
						qual += 1;
					}
				}
			}else{
				loop {
					let _data = self.reader.fill_buf()?; // for multiline fasta
					if !_data.is_empty() && _data[0] != b'>' {
						self.data.pop();
						seq += self.reader.read_until(b'\n', &mut self.data)? - 1;
					}else {
						break;
					}
				}
				if self.data[seq - 1] != b'\n' { // for files not ending with '\n'
					seq += 1;
				}
				sep = seq; //reset sep & qual
				qual = sep;
			}
		}

		if seq == des || !is_fasta && (sep == seq || qual == sep) {
			return Err(ParseError::TruncateFile(String::from_utf8(self.data.to_owned()).unwrap()));
		}
		//println!("head:{:?} des {:?} seq {:?} sep {:?} qual {:?}", head,des,seq,sep,qual);
		let fastx = Fastx {
			_head: head,
			_des: des,
			_seq: seq,
			_sep: sep,
			_qual: qual,
			_data: &self.data
			};

		if is_fasta {
			if !fastx.validate_fasta(){
				return Err(ParseError::InvalidFasta(fastx.head().to_string()));
			}
		}else if !fastx.validate_fastq() {
			return Err(ParseError::InvalidFastq(fastx.head().to_string()));
		}
		Ok(Some(fastx))
	}

	pub fn iter_record_check(&mut self) -> Result<Option<Fastx>> {
		if !self.reach_eof()? {
			return self.iter_record();
		}
		Ok(None)
	}
}

/// multiple readers for a fofn file
pub struct Readers {
	index: usize,
	pub readers: Vec<Reader>
}

impl Readers {

	/// create a new Readers
	pub fn new() -> Self {
		Readers {
			index: 0,
			readers: Vec::new()
		}
	}

	/// iterate over a record from this Readers
	pub fn iter_record(&mut self) -> Result<Option<Fastx>> {
		for idx in self.index..self.readers.len() {
			if !self.readers[idx].reach_eof()? {
				return self.readers[idx].iter_record();
			}
			self.index += 1;
		}
		Ok(None)
	}
}
