use std::{io::{BufRead, Result}, str::from_utf8_unchecked};

pub trait Seq {

    fn head(&self) -> &str;

    fn seq(&self) -> &str;

    fn sep(&self) -> &str;

    fn qual(&self) -> &str;

    fn len(&self) ->usize;

    fn validate_format(&self)-> bool;

    fn validate_dna(&self) -> bool {
        self.seq().bytes().all(|x| x == b'A' || x == b'C' || x == b'T' || x == b'G')
    }

    fn validate_dnan(&self) -> bool {
        self.seq().bytes().all(|x| x == b'A' || x == b'C' || x == b'T' || x == b'G' || x == b'N')
    }
}

pub struct Fastq<'a> {
    _head: usize,
    _seq: usize,
    _sep: usize,
    _qual: usize,
    _data: &'a [u8]
}

impl Seq for Fastq<'_> {
    #[inline]
    fn head(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[1 .. self._head - 1]) }
    }

    #[inline]
    fn seq(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[self._head .. self._seq - 1]) }
    }

    #[inline]
    fn sep(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[self._seq .. self._sep - 1]) }
    }

    #[inline]
    fn qual(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[self._sep .. self._qual - 1]) }
    }

    fn len(&self) -> usize{
        self._seq - self._head - 1
    }
    
    #[inline]
    fn validate_format(&self)-> bool {
        if self._seq - self._head != self._qual - self._sep || self._data[0] != b'@' {
            false
        }else {true}
    }
}

pub struct Fasta<'a> {
    _head: usize,
    _seq: usize,
    _data: &'a [u8]
}

impl Seq for Fasta<'_> {
    #[inline]
    fn head(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[1 .. self._head - 1]) }
    }

    #[inline]
    fn seq(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[self._head .. self._seq - 1]) }
    }


    #[inline]
    fn sep(&self) -> &str {
        "\0"
    }

    #[inline]
    fn qual(&self) -> &str {
        "\0"
    }

    fn len(&self) -> usize{
        self._seq - self._head - 1
    }
    
    #[inline]
    fn validate_format(&self)-> bool {
        if self._data[0] != b'>' {
            false
        }else {true}
    }
}

pub enum Record<'a> {
    Fastq(Fastq<'a>),
    Fasta(Fasta<'a>)
}

impl Seq for Record<'_> {
    #[inline]
    fn head(&self) -> &str {
        match self {
            Record::Fasta(t) => t.head(),
            Record::Fastq(t) => t.head(),
        }
    }

    #[inline]
    fn seq(&self) -> &str {
        match self {
            Record::Fasta(t) => t.seq(),
            Record::Fastq(t) => t.seq(),
        }
    }


    #[inline]
    fn sep(&self) -> &str {
        match self {
            Record::Fasta(t) => t.sep(),
            Record::Fastq(t) => t.sep(),
        }
    }

    #[inline]
    fn qual(&self) -> &str {
        match self {
            Record::Fasta(t) => t.qual(),
            Record::Fastq(t) => t.qual(),
        }
    }

    fn len(&self) -> usize{
        match self {
            Record::Fasta(t) => t.len(),
            Record::Fastq(t) => t.len(),
        }
    }
    
    #[inline]
    fn validate_format(&self)-> bool {
        match self {
            Record::Fasta(t) => t.validate_format(),
            Record::Fastq(t) => t.validate_format(),
        }
    }
}

pub struct Reader {
    reader: Box<dyn BufRead>,
    data: String,
    is_fastq: bool
}

impl Reader {
    pub fn new(r: Box<dyn BufRead>, is_fastq: bool) -> Self {
        Reader {
            reader: r,
            data: String::with_capacity(1024),
            is_fastq: is_fastq
        }
    }

    #[inline]
    pub fn iter_record(&mut self) -> Option<Record> {
        if self.is_fastq {
            Some(Record::Fastq(self.iter_fastq()))
        }else {
            Some(Record::Fasta(self.iter_fasta()))
        }
    }

    fn is_end(&mut self) -> bool {
        self.data.clear();
        loop {
            let head = self.reader.read_line(&mut self.data).expect("fail read record!");
            if head == 0 {return true;}
            if self.data.trim().is_empty() {continue;}
            break;
        }
        return false;
    }

    #[inline]
    fn iter_fastq(&mut self) -> Fastq {

        let head = self.data.len();
        let seq = head + self.reader.read_line(&mut self.data).expect("fail read record!");
        let sep = seq + self.reader.read_line(&mut self.data).expect("fail read record!");
        let qual = sep + self.reader.read_line(&mut self.data).expect("fail read record!");
        if seq == head || sep == seq || qual == sep {
            panic!("truncated fastq file");
        }
        let fastq = Fastq {
            _head: head,
            _seq: seq,
            _sep: sep,
            _qual: qual,
            _data: self.data.as_bytes()
            };

        println!("{:?} {:?}", head, seq);
        
        if !fastq.validate_format() {
            panic!("Not a validate fastq record {}.", fastq.head());
        }
        fastq
    }

    #[inline]
    fn iter_fasta(&mut self) -> Fasta {
        let head = self.data.len();
        let seq = head + self.reader.read_line(&mut self.data).expect("fail read record!");
        if seq == head {
            panic!("truncated fasta file");
        }
        let fasta = Fasta {
            _head: head,
            _seq: seq,
            _data: self.data.as_bytes()
            };
        if !fasta.validate_format() {
            panic!("Not a validate fastq record {}.", fasta.head());
        }
        fasta
    }
}

// impl<'a> Iterator for Record {
//     type Item = Fastq<'a>;

//     #[inline]
//     fn next(&mut self) -> Option<Self::Item> {
//         self.data.clear();
//         let head = self.reader.read_line(&mut self.data).expect("fail read record!");
//         let seq = head + self.reader.read_line(&mut self.data).expect("fail read record!");
//         let sep = seq + self.reader.read_line(&mut self.data).expect("fail read record!");
//         let qual = sep + self.reader.read_line(&mut self.data).expect("fail read record!");
//         Some(Fastq {
//             head: head,
//             seq: seq,
//             sep: sep,
//             qual: qual,
//             data: self.data.as_bytes()
//         })
//     }
// }

pub struct Paths {
    index: usize,
    pub paths: Vec<Result<Reader>>
}

impl Paths {
    pub fn new() -> Self {
        Paths {
            index: 0,
            paths: Vec::new()
        }
    }

    pub fn iter_record(&mut self) -> Option<Record> {
        for idx in self.index..self.paths.len() {
            if !self.paths[idx].as_mut().unwrap().is_end() {
                match self.paths[idx].as_mut().unwrap().iter_record() {
                    Some(r) => return Some(r),
                    None => unreachable!(),
                }
            }else {
                self.index += 1;
            }
        }
        None
    }
}