//! This module contains the function to encode the input string into a compressed byte vector.

use super::{
    CharacterSet,
    Encode
};

/// Encodes the input string into a compressed byte vector.
///
/// The input string is expected to contain annotations separated by semicolons (;).
/// The annotations are categorized into three types: InterPro (IPR), Gene Ontology (GO), and Enzyme
/// Commission (EC). The function splits the input string into these annotation types and encodes
/// them into a compressed byte vector.
///
/// # Arguments
///
/// * `input` - The input string containing annotations.
///
/// # Returns
///
/// A compressed byte vector representing the encoded annotations.
///
/// # Examples
///
/// ```
/// use fa_compression::algorithm1::encode;
///
/// let input = "IPR:IPR016364;EC:1.1.1.-;GO:0009279";
/// let encoded = encode(input);
///
/// assert_eq!(encoded, vec![ 44, 44, 44, 189, 17, 26, 56, 173, 18, 116, 117 ]);
/// ```
pub fn encode(input: &str) -> Vec<u8> {
    if input.is_empty() {
        return Vec::new();
    }

    // ==========================================================================================
    // !!!!! The code between the equal signs can be removed if the input is already sorted !!!!!
    // ==========================================================================================

    // Create vectors to store the different types of annotations
    // Note: We assume an average of 12 characters per annotation
    //       We also know that we have 3 types of annotations
    //       So we can pre-allocate a vector with a capacity of input.len() / 12 / 3
    let mut interpros: Vec<&str> = Vec::with_capacity(input.len() / 36);
    let mut gos: Vec<&str> = Vec::with_capacity(input.len() / 36);
    let mut ecs: Vec<&str> = Vec::with_capacity(input.len() / 36);

    // Read the input and split the annotations into the corresponding vectors
    for annotation in input.split(';') {
        if annotation.starts_with("IPR") {
            interpros.push(&annotation[7 ..]);
        } else if annotation.starts_with("GO") {
            gos.push(&annotation[3 ..]);
        } else if annotation.starts_with("EC") {
            ecs.push(&annotation[3 ..]);
        }
    }

    // ==========================================================================================
    // ==========================================================================================
    // ==========================================================================================

    // Create a string without any unnecessary characters
    let mut result = String::with_capacity(input.len());
    result.push_str(&ecs.join(";"));
    result.push(',');
    result.push_str(&gos.join(";"));
    result.push(',');
    result.push_str(&interpros.join(";"));

    // Take two characters at a time and encode them into a single byte
    let mut encoded: Vec<u8> = Vec::with_capacity(result.len() / 2);
    for bytes in result.as_bytes().chunks(2) {
        if bytes.len() == 2 {
            encoded.push(CharacterSet::encode(bytes[0]) | CharacterSet::encode(bytes[1]));
        } else {
            encoded.push(CharacterSet::encode(bytes[0]) | CharacterSet::Empty);
        }
    }

    encoded
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_empty() {
        assert_eq!(encode(""), vec![])
    }

    #[test]
    fn test_encode_single_ec() {
        assert_eq!(encode("EC:1.1.1.-"), vec![44, 44, 44, 189, 208])
    }

    #[test]
    fn test_encode_single_go() {
        assert_eq!(encode("GO:0009279"), vec![209, 17, 163, 138, 208])
    }

    #[test]
    fn test_encode_single_ipr() {
        assert_eq!(encode("IPR:IPR016364"), vec![221, 18, 116, 117])
    }

    #[test]
    fn test_encode_no_ec() {
        assert_eq!(
            encode("IPR:IPR016364;GO:0009279;IPR:IPR008816"),
            vec![209, 17, 163, 138, 209, 39, 71, 94, 17, 153, 39]
        )
    }

    #[test]
    fn test_encode_no_go() {
        assert_eq!(
            encode("IPR:IPR016364;EC:1.1.1.-;EC:1.2.1.7"),
            vec![44, 44, 44, 190, 44, 60, 44, 141, 209, 39, 71, 80]
        )
    }

    #[test]
    fn test_encode_no_ipr() {
        assert_eq!(
            encode("EC:1.1.1.-;GO:0009279;GO:0009279"),
            vec![44, 44, 44, 189, 17, 26, 56, 174, 17, 26, 56, 173]
        )
    }

    #[test]
    fn test_encode_all() {
        assert_eq!(
            encode("IPR:IPR016364;EC:1.1.1.-;IPR:IPR032635;GO:0009279;IPR:IPR008816"),
            vec![44, 44, 44, 189, 17, 26, 56, 173, 18, 116, 117, 225, 67, 116, 110, 17, 153, 39]
        )
    }
}
