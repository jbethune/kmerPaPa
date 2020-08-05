use std::fmt;

#[derive(Debug)]
pub struct MissingCommandLineArgumentError {
    argument: &'static str,
}

impl MissingCommandLineArgumentError {
    pub fn new(argument: &'static str) -> Self {
        Self { argument }
    }
}

impl fmt::Display for MissingCommandLineArgumentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Please provide the command line argument {}",
            self.argument
        )
    }
}

impl std::error::Error for MissingCommandLineArgumentError {}

#[derive(Debug)]
pub struct ParseError {
    message: String,
}

impl ParseError {
    pub fn new(message: String) -> Self {
        Self { message }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for ParseError {}
