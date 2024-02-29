use std::error::Error;
use std::fmt::{Display, Formatter};

/// Custom error type to indicate that loading an existing SA failed
#[derive(Debug)]
pub struct IndexLoadError {
    message: String,
}

impl Error for IndexLoadError {}

impl Display for IndexLoadError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl IndexLoadError {
    pub fn new(message: &str) -> Self {
        Self {
            message: message.to_string(),
        }
    }
}