pub mod taxon_id_calculator;

use std::error::Error;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::Path;
use umgap::taxon::TaxonId;
use crate::taxon_id_calculator::{TaxonIdVerifier};

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
    pub input_string: Vec<u8>,
    pub proteins: Vec<Protein>
}

impl Proteins {

    pub fn get_sequence(&self, protein: &Protein) ->&str {
        let (begin, size) = protein.sequence;
        std::str::from_utf8( &self.input_string[begin..begin + size as usize]).unwrap() // should never fail since the input string will always be utf8
    }
}


/// The useful information about a protein for our use case
pub struct Protein {
    pub uniprot_id: String,
    pub sequence: (usize, u32),
    pub id: TaxonId,
    pub ec_numbers: Vec<String>,
    pub go_terms: Vec<String>,
    pub interpro: Vec<String>
}

/// Aid function that splits the input string on ';' as separator and returns a vector with the resulting strings
fn csv_split(data: String) -> Vec<String> {
    // empty string should return empty vector, and not vector with empty string inside it
    if data.is_empty() {
        Vec::new()
    } else {
        data.split(';').map(str::to_string).collect()
    }
}


/// Parse the given database tsv file into a Vector of Proteins with the data from the tsv file
pub fn get_proteins_from_database_file(database_file: &str, taxon_id_calculator: &dyn TaxonIdVerifier) -> Result<Proteins, Box<dyn Error>> {
    let mut input_string: String = "".to_string();
    let mut proteins: Vec<Protein> = vec![];
    let mut begin_index: usize = 0;
    let lines = read_lines(database_file)?;
    for line in lines.into_iter().map_while(Result::ok) {
        let parts: Vec<String> = line.split('\t').map(str::to_string).collect();
        // TODO: description field is wrong, what is called description here is actually the protein name
        let [_unipept_id, uniprot_id, _, protein_id_str, _, _protein_name, protein_sequence, ec_numbers, go_terms, inter_pros]: [String; 10] = parts.try_into().map_err(|e| DatabaseFormatError{ error: e})?;
        let protein_id_as_taxon_id = protein_id_str.parse::<TaxonId>()?;
        // if the taxon ID is not a valid ID in our NCBI taxonomy, skip this protein
        if !taxon_id_calculator.taxon_id_exists(protein_id_as_taxon_id) {
            // eprintln!("Skipped protein with taxon id {}!", protein_id_as_taxon_id);
            continue;
        }

        if begin_index != 0 {
            input_string.push(SEPARATION_CHARACTER as char);
        }
        input_string.push_str(&protein_sequence.to_uppercase());
        proteins.push(
            Protein {
                uniprot_id,
                sequence: (begin_index, protein_sequence.len() as u32),
                id: protein_id_as_taxon_id,
                ec_numbers: csv_split(ec_numbers),
                go_terms: csv_split(go_terms),
                interpro: csv_split(inter_pros)
            }
        );
        begin_index += protein_sequence.len() + 1;
    }
    input_string.push(END_CHARACTER as char);
    Ok(Proteins {
        input_string: input_string.into_bytes(),
        proteins
    })
}

#[derive(Debug)]
struct DatabaseFormatError {
    error: Vec<String>
}

impl std::fmt::Display for DatabaseFormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Expected the protein database file to have the following fields separated by a tab: <unipept id> <Uniprot_id> <TODO> <protein id> <TODO> <protein_name> <sequence> <ec numbers> <go_terms> <inter_pros>.\n The given error is: {:?}", self.error) // TODO: fill in all the fields
    }
}

impl Error for DatabaseFormatError {}

#[cfg(test)]
mod tests {
    use crate::csv_split;

    #[test]
    fn test_csv_split() {
        let data = "a;ab;abc;abcd".to_string();
        assert_eq!(vec!["a", "ab", "abc", "abcd"], csv_split(data));
    }

    #[test]
    fn test_csv_split_empty() {
        let data = String::new();
        assert_eq!(vec![] as Vec<String>, csv_split(data));
    }
}