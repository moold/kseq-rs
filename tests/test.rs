use kseq;
use std::io::Cursor;
use std::result::Result;

// https://stackoverflow.com/questions/53124930/how-do-you-test-for-a-specific-rust-error
macro_rules! assert_err {
    ($expression:expr, $($pattern:tt)+) => {
        match $expression {
            $($pattern)+ => (),
            ref e => panic!("expected `{}` but got `{:?}`", stringify!($($pattern)+), e),
        }
    }
}

fn count_base(input: Vec<u8>) -> Result<usize, kseq::record::ParseError> {
    let mut len_bases = 0;
    let mut seq_bases = 0;
    let mut records = kseq::Paths::Reader(kseq::record::Reader::new(Box::new(Cursor::new(input))));
    while let Some(record) = records.iter_record()? {
        len_bases += record.len();
        seq_bases += record.seq().len();
    }
    assert_eq!(len_bases, seq_bases);
    Ok(len_bases)
}

static BASE_SEQ: &str = "ATGCATGCATGC";
static BASE_QUAL: &str = "@@@@@@@@@@@@";

#[test]
fn test_normal_one_line_fasta() {
    let data: Vec<u8> =
        format!(">1 record1\n{seq}\n>2 record2\n{seq}", seq = BASE_SEQ).into_bytes();
    assert_eq!(BASE_SEQ.len() * 2, count_base(data).unwrap());
}

#[test]
fn test_normal_one_line_fastq() {
    let data: Vec<u8> = format!(
        "@1 record1\n{seq}\n+\n{qual}\n@2 record2\n{seq}\n+\n{qual}",
        seq = BASE_SEQ,
        qual = BASE_QUAL
    )
    .into_bytes();
    assert_eq!(BASE_SEQ.len() * 2, count_base(data).unwrap());
}

#[test]
fn test_normal_multi_line_fasta() {
    let data: Vec<u8> = format!(
        ">1 record1\n{seq}{seq}\n>2 record2\n{seq}{seq}",
        seq = BASE_SEQ
    )
    .into_bytes();
    assert_eq!(BASE_SEQ.len() * 4, count_base(data).unwrap());
}

#[test]
fn test_normal_multi_line_fastq() {
    let data: Vec<u8> = format!(
        "@1 record1\n{seq}{seq}\n+\n{qual}{qual}\n@2 record2\n{seq}{seq}\n+\n{qual}{qual}\n",
        seq = BASE_SEQ,
        qual = BASE_QUAL
    )
    .into_bytes();
    assert_eq!(BASE_SEQ.len() * 4, count_base(data).unwrap());
}

#[test]
fn test_truncate_fasta() {
    let data: Vec<u8> = format!(">1 record1\n{seq}{seq}\n>2 record2", seq = BASE_SEQ).into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::TruncateFile(_))
    );
}

#[test]
fn test_truncate_fastq_miss_seq() {
    let data: Vec<u8> = format!("@1 record1").into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::TruncateFile(_))
    );
}

#[test]
fn test_truncate_fastq_miss_sep() {
    let data: Vec<u8> = format!("@1 record1\n{seq}{seq}", seq = BASE_SEQ).into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::TruncateFile(_))
    );
}

#[test]
fn test_truncate_fastq_miss_qual() {
    let data: Vec<u8> = format!("@1 record1\n{seq}{seq}\n+", seq = BASE_SEQ).into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::TruncateFile(_))
    );
}

#[test]
fn test_invalid_fasta() {
    let data: Vec<u8> =
        format!("+1 record1\n{seq}\n>2 record2\n{seq}", seq = BASE_SEQ).into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::InvalidFasta(_))
    );
}

#[test]
fn test_invalid_fasta_with_seq_len_is_0() {
    let data: Vec<u8> = format!(">1 record1\n{seq}\n>2 record2\n{seq}", seq = "").into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::InvalidFasta(_))
    );
}

#[test]
fn test_invalid_fastq_seq_has_diff_len_with_qual() {
    let data: Vec<u8> = format!(
        "@1 record1\n{seq}{seq}\n+\n{qual}{qual}\n@2 record2\n{seq}{seq}\n+\n{qual}{qual}\n",
        seq = BASE_SEQ,
        qual = &BASE_QUAL[0..BASE_SEQ.len() - 1]
    )
    .into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::InvalidFastq(_))
    );
}

#[test]
fn test_invalid_fastq_with_seq_len_is_0() {
    let data: Vec<u8> = format!(
        "@1 record1\n{seq}\n+\n{qual}\n@2 record2\n{seq}\n+\n{qual}\n",
        seq = "",
        qual = ""
    )
    .into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::InvalidFastq(_))
    );
}

#[test]
fn test_invalid_fastq_miss_sep() {
    let data: Vec<u8> = format!(
        "@1 record1\n{seq}{seq}\n{qual}{qual}\n@2 record2\n{seq}{seq}\n+\n{qual}{qual}\n",
        seq = BASE_SEQ,
        qual = &BASE_QUAL[0..BASE_SEQ.len() - 1]
    )
    .into_bytes();
    assert_err!(
        count_base(data),
        Err(kseq::record::ParseError::InvalidFastq(_))
    );
}
