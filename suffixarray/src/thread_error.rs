use std::error::Error;
use std::fmt::{Display, Formatter};

/// Custom error type to indicate when one of the threads that computes the results failed
#[derive(Debug)]
pub struct ThreadError {
    message: String,
}

impl Error for ThreadError {}

impl Display for ThreadError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl ThreadError {
    pub fn new(message: &str) -> Self {
        Self {
            message: message.to_string(),
        }
    }
}
