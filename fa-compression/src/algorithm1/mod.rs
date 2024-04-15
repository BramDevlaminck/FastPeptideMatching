//! The `fa-compression` crate provides functions to encode and decode annotations following a
//! specific format

use std::ops::BitOr;

mod decode;
mod encode;

pub use decode::decode;
pub use encode::encode;

/// Trait for encoding a value into a character set.
trait Encode {
    /// Encodes the given value into a character set.
    ///
    /// # Arguments
    ///
    /// * `value` - The value to be encoded.
    ///
    /// # Returns
    ///
    /// The encoded character set.
    fn encode(value: u8) -> CharacterSet;
}

/// Trait for decoding a value from a character set.
trait Decode {
    /// Decodes the given value from a character set into a character.
    ///
    /// # Arguments
    ///
    /// * `value` - The value to be decoded.
    ///
    /// # Returns
    ///
    /// The decoded character.
    fn decode(value: u8) -> char;

    /// Decodes a pair of values from a character set into a pair of characters.
    ///
    /// # Arguments
    ///
    /// * `value` - The value to be decoded.
    ///
    /// # Returns
    ///
    /// A tuple containing the decoded characters.
    fn decode_pair(value: u8) -> (char, char) {
        (Self::decode(value >> 4), Self::decode(value & 0b1111))
    }
}

/// Enum representing the set of characters that can be encoded.
#[repr(u8)]
#[cfg_attr(test, derive(Clone, Copy))]
#[derive(PartialEq, Eq, Debug)]
enum CharacterSet {
    /// Empty placeholder character
    Empty,

    /// Numeric characters
    Zero,
    One,
    Two,
    Three,
    Four,
    Five,
    Six,
    Seven,
    Eight,
    Nine,

    /// Special Enzyme Commission characters
    Dash,
    Point,

    /// Different annotation type separator
    Comma,

    /// Annotation separator
    Semicolon
}

impl Encode for CharacterSet {
    /// Encodes the given value into a character set.
    ///
    /// # Arguments
    ///
    /// * `value` - The value to be encoded.
    ///
    /// # Returns
    ///
    /// The encoded character set.
    fn encode(value: u8) -> CharacterSet {
        match value {
            b'$' => CharacterSet::Empty,
            b'0' => CharacterSet::Zero,
            b'1' => CharacterSet::One,
            b'2' => CharacterSet::Two,
            b'3' => CharacterSet::Three,
            b'4' => CharacterSet::Four,
            b'5' => CharacterSet::Five,
            b'6' => CharacterSet::Six,
            b'7' => CharacterSet::Seven,
            b'8' => CharacterSet::Eight,
            b'9' => CharacterSet::Nine,
            b'-' => CharacterSet::Dash,
            b'.' => CharacterSet::Point,
            b',' => CharacterSet::Comma,
            b';' => CharacterSet::Semicolon,
            _ => panic!("Invalid character")
        }
    }
}

impl Decode for CharacterSet {
    /// Decodes the given value from a character set into a character.
    ///
    /// # Arguments
    ///
    /// * `value` - The value to be decoded.
    ///
    /// # Returns
    ///
    /// The decoded character.
    fn decode(value: u8) -> char {
        match value {
            0 => '$',
            1 => '0',
            2 => '1',
            3 => '2',
            4 => '3',
            5 => '4',
            6 => '5',
            7 => '6',
            8 => '7',
            9 => '8',
            10 => '9',
            11 => '-',
            12 => '.',
            13 => ',',
            14 => ';',
            _ => panic!("Invalid character")
        }
    }
}

impl BitOr for CharacterSet {
    type Output = u8;

    /// Performs a bitwise OR operation between two character sets.
    ///
    /// # Arguments
    ///
    /// * `self` - The left-hand side character set.
    /// * `rhs` - The right-hand side character set.
    ///
    /// # Returns
    ///
    /// The result of the bitwise OR operation.
    fn bitor(self, rhs: Self) -> Self::Output {
        ((self as u8) << 4) | rhs as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static CHARACTERS: [u8; 15] =
        [b'$', b'0', b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9', b'-', b'.', b',', b';'];

    static CHARACTER_SETS: [CharacterSet; 15] = [
        CharacterSet::Empty,
        CharacterSet::Zero,
        CharacterSet::One,
        CharacterSet::Two,
        CharacterSet::Three,
        CharacterSet::Four,
        CharacterSet::Five,
        CharacterSet::Six,
        CharacterSet::Seven,
        CharacterSet::Eight,
        CharacterSet::Nine,
        CharacterSet::Dash,
        CharacterSet::Point,
        CharacterSet::Comma,
        CharacterSet::Semicolon
    ];

    #[test]
    fn test_or() {
        for i in 0 .. CHARACTERS.len() {
            for j in 0 .. CHARACTERS.len() {
                assert_eq!(CHARACTER_SETS[i] | CHARACTER_SETS[j], ((i as u8) << 4) | (j as u8));
            }
        }
    }

    #[test]
    fn test_encode() {
        for i in 0 .. CHARACTERS.len() {
            assert_eq!(CHARACTER_SETS[i], CharacterSet::encode(CHARACTERS[i]));
        }
    }

    #[test]
    fn test_decode() {
        for (i, c) in CHARACTERS.iter().enumerate() {
            assert_eq!(CharacterSet::decode(i as u8), *c as char);
        }
    }

    #[test]
    fn test_decode_pair() {
        for (i1, c1) in CHARACTERS.iter().enumerate() {
            for (i2, c2) in CHARACTERS.iter().enumerate() {
                let encoded = CharacterSet::encode(*c1) | CharacterSet::encode(*c2);
                assert_eq!(
                    CharacterSet::decode_pair(encoded),
                    (CharacterSet::decode(i1 as u8), CharacterSet::decode(i2 as u8))
                );
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_encode_invalid() {
        CharacterSet::encode(b'A');
    }

    #[test]
    #[should_panic]
    fn test_decode_invalid() {
        CharacterSet::decode(15);
    }

    #[test]
    #[should_panic]
    fn test_decode_pair_invalid() {
        CharacterSet::decode_pair(0b11111111);
    }

    #[test]
    fn test_clone() {
        let character_set = CharacterSet::Empty;
        let character_set_clone = character_set.clone();
        assert_eq!(character_set, character_set_clone);
    }

    #[test]
    fn test_copy() {
        let character_set = CharacterSet::Empty;
        let character_set_copy = character_set;
        assert_eq!(character_set, character_set_copy);
    }
}
