use std::fs::{DirBuilder, File};
use std::io::{stdin, stdout};
use std::io::{Read, Write};
use std::path::Path;

use anyhow::Result;

use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;

fn is_terminal_io_file(path: &str) -> bool {
    path == "-" || path == "/dev/stdout" || path == "/dev/stdin"
}

pub fn get_writer(path: &str) -> Result<Box<dyn Write>> {
    if is_terminal_io_file(path) {
        Ok(Box::new(stdout()))
    } else {
        if let Some(parent) = Path::new(path).parent() {
            DirBuilder::new().recursive(true).create(parent)?;
        }
        let out = File::create(&path)?;
        if path.ends_with(".gz") {
            Ok(Box::new(GzEncoder::new(out, Compression::default())))
        } else {
            Ok(Box::new(out))
        }
    }
}

pub fn get_reader(path: &str) -> Result<Box<dyn Read>> {
    if is_terminal_io_file(path) {
        Ok(Box::new(stdin()))
    } else {
        let fd = File::open(&path)?;
        if path.ends_with(".gz") {
            Ok(Box::new(GzDecoder::new(fd)))
        } else {
            Ok(Box::new(fd))
        }
    }
}
