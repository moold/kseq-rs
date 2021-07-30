use std::{io::BufRead, str::from_utf8_unchecked};

pub struct Fastx<'a> {
    _head: usize,
    _des: usize,
    _seq: usize,
    _sep: usize,
    _qual: usize,
    _data: &'a [u8]
}

impl Fastx<'_> {
    #[inline]
    pub fn head(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[1 .. self._head - 1]) }
    }

    #[inline]
    pub fn seq(&self) -> &str {
        unsafe{ from_utf8_unchecked(&self._data[self._des .. self._seq - 1]) }
    }

    #[inline]
    pub fn des(&self) -> &str {
        if self._des != self._head {
            unsafe{ from_utf8_unchecked(&self._data[self._head .. self._des - 1]) }
        }else {
            ""
        }
    }

    #[inline]
    pub fn sep(&self) -> &str {
        if self._seq != self._sep{
            unsafe{ from_utf8_unchecked(&self._data[self._seq .. self._sep - 1]) }
        }else{
            ""
        }
    }

    #[inline]
    pub fn qual(&self) -> &str {
        if self._sep != self._qual {
            unsafe{ from_utf8_unchecked(&self._data[self._sep .. self._qual - 1]) }
        }else{
            ""
        }
    }
    
    #[inline]
    pub fn len(&self) -> usize {
        self._seq - self._des - 1
    }
    
    fn validate_fastq(&self)-> bool {
        if self._seq - self._des != self._qual - self._sep || self._data[0] != b'@' {
            false
        }else {true}
    }

    fn validate_fasta(&self)-> bool {
        if self._data[0] != b'>' {
            false
        }else {true}
    }

    pub fn validate_dna(&self) -> bool {
        self.seq().bytes().all(|x| x == b'A' || x == b'C' || x == b'T' || x == b'G')
    }

    pub fn validate_dnan(&self) -> bool {
        self.seq().bytes().all(|x| x == b'A' || x == b'C' || x == b'T' || x == b'G' || x == b'N')
    }
}

pub struct Reader {
    reader: Box<dyn BufRead>,
    data: String,
}

impl Reader {
    pub fn new(r: Box<dyn BufRead>) -> Self {
        Reader {
            reader: r,
            data: String::with_capacity(1024),
        }
    }

    fn reach_eof(&mut self) -> bool {
        loop {
            self.data.clear();
            let head = self.reader.read_line(&mut self.data).expect("fail read record!");
            if head == 0 {return true;}
            if self.data.trim().is_empty() {continue;}
            break;
        }
        return false;
    }

    pub fn iter_record(&mut self) -> Option<Fastx> {

        let head = self.data.find(char::is_whitespace).unwrap() + 1;
        let des = self.data.len();
        let seq = des + self.reader.read_line(&mut self.data).expect("fail read record!");
        let mut sep = seq;
        let mut qual = seq;
        let mut is_fasta = true;
        if self.data.starts_with('@'){
            is_fasta = false;
            sep = seq + self.reader.read_line(&mut self.data).expect("fail read record!");
            qual = sep + self.reader.read_line(&mut self.data).expect("fail read record!");
            if seq == des || sep == seq || qual == sep {
                panic!("truncated fastq file, problematic record: {:?}", self.data);
            }
        }else if seq == des {
            panic!("truncated fasta file, problematic record: {:?}", self.data);
        }
        // println!("head:{:?} des {:?} seq {:?} qual {:?}", head,des,seq,qual);
        let fastx = Fastx {
            _head: head,
            _des: des,
            _seq: seq,
            _sep: sep,
            _qual: qual,
            _data: self.data.as_bytes()
            };

        if is_fasta {
            if !fastx.validate_fasta(){
                panic!("Not a validate fasta record {}.", fastx.head());
            }
        }else if !fastx.validate_fastq() {
            panic!("Not a validate fastq record {}.", fastx.head());
        }
        Some(fastx)
    }

    pub fn iter_record_check(&mut self) -> Option<Fastx> {
        if !self.reach_eof() {
            return self.iter_record();
        }
        None
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

pub struct Readers {
    index: usize,
    pub readers: Vec<Reader>
}

impl Readers {
    pub fn new() -> Self {
        Readers {
            index: 0,
            readers: Vec::new()
        }
    }

    pub fn iter_record(&mut self) -> Option<Fastx> {
        for idx in self.index..self.readers.len() {
            if !self.readers[idx].reach_eof() {
                return self.readers[idx].iter_record();
            }
            self.index += 1;
        }
        None
    }
}