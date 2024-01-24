pub mod taxon_id_calculator;

use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::Path;
use umgap::taxon::TaxonId;

// END_CHARACTER should ALWAYS be lexicographically than SEPARATION_CHARACTER
// otherwise the building of the suffix array will not happen correctly
pub static SEPARATION_CHARACTER: u8 = b'-';
pub static END_CHARACTER: u8 = b'$';


// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub struct Proteins {
    pub input_string: String,
    pub proteins: Vec<Protein>
}

impl Proteins {

    pub fn get_sequence(&self, protein: &Protein) -> &str {
        let (begin, size) = protein.sequence;
        &self.input_string[begin..begin + size as usize]
    }
}


/// The useful information about a protein for our use case
pub struct Protein {
    pub sequence: (usize, u32),
    pub id: TaxonId,
}

/// Parse the given database tsv file into a Vector of Proteins with the data from the tsv file
pub fn get_proteins_from_database_file(database_file: &str) -> Proteins {
    let mut input_string: String = String::new();
    let mut proteins: Vec<Protein> = vec![];
    let mut begin_index: usize = 0;
    if let Ok(lines) = read_lines(database_file) {
        for line in lines.into_iter().flatten() {
            let [_, _, protein_id_str, _, _, protein_sequence]: [&str; 6] = line.splitn(6, '\t').collect::<Vec<&str>>().try_into().unwrap();
            let protein_id_as_taxon_id = protein_id_str.parse::<TaxonId>().expect("Could not parse id of protein to usize!");
            input_string.push_str(&protein_sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER as char);
            proteins.push(
                Protein {
                    sequence: (begin_index, protein_sequence.len() as u32),
                    id: protein_id_as_taxon_id,
                }
            );
            begin_index += protein_sequence.len() + 1;
        }
    } else {
        eprintln!("Database file {} could not be opened!", database_file);
        std::process::exit(1);
    }
    // change the last character from SEPARATION_CHARACTER to END_CHARACTER
    input_string.pop();
    input_string.push(END_CHARACTER as char);
    Proteins {
        input_string,
        proteins
    }
}

#[cfg(test)]
mod tests {
    use crate::{END_CHARACTER, get_proteins_from_database_file, SEPARATION_CHARACTER};

    #[test]
    fn test_protein_read() {
        let proteins = get_proteins_from_database_file("../testfiles/test_database.tsv");
        let sequences_from_file = ["ac", "bracvaa", "ac", "bacrqz"];
        let mut expected_input_string = sequences_from_file
            .map(|s| s.to_uppercase())
            .join(&format!("{}", SEPARATION_CHARACTER as char));
        expected_input_string.push( END_CHARACTER as char);
        assert_eq!(proteins.input_string, expected_input_string);
    }

    #[test]
    fn test_protein_empty() {
        let proteins = get_proteins_from_database_file("../testfiles/empty_database.tsv");
        assert_eq!(proteins.input_string, format!("{}", END_CHARACTER as char));
    }
}