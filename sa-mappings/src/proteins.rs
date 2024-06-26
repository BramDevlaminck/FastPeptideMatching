//! This module contains the `Protein` and `Proteins` structs, which are used to represent proteins
//! and collections of proteins, respectively.

use std::{
    error::Error,
    fs::File,
    io::BufReader,
    ops::Index,
    str::from_utf8
};

use bytelines::ByteLines;
use fa_compression::algorithm1::decode;
use umgap::taxon::TaxonId;

use crate::taxonomy::TaxonAggregator;

/// The separation character used in the input string
pub static SEPARATION_CHARACTER: u8 = b'-';

/// The termination character used in the input string
/// This character should be smaller than the separation character
pub static TERMINATION_CHARACTER: u8 = b'$';

/// A struct that represents a protein and its linked information
pub struct Protein {
    /// The id of the protein
    pub uniprot_id: String,

    /// the taxon id of the protein
    pub taxon_id: TaxonId,

    /// The encoded functional annotations of the protein
    pub functional_annotations: Vec<u8>
}

/// A struct that represents a collection of proteins
pub struct Proteins {
    /// The input string containing all proteins
    pub input_string: Vec<u8>,

    /// The proteins in the input string
    pub proteins: Vec<Protein>
}

impl Protein {
    /// Returns the decoded functional annotations of the protein
    pub fn get_functional_annotations(&self) -> String {
        decode(&self.functional_annotations)
    }
}

impl Proteins {
    /// Creates a new `Proteins` struct from a database file and a `TaxonAggregator`
    ///
    /// # Arguments
    /// * `file` - The path to the database file
    /// * `taxon_aggregator` - The `TaxonAggregator` to use
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the `Proteins` struct
    ///
    /// # Errors
    ///
    /// Returns a `Box<dyn Error>` if an error occurred while reading the database file
    pub fn try_from_database_file(
        file: &str,
        taxon_aggregator: &TaxonAggregator
    ) -> Result<Self, Box<dyn Error>> {
        let mut input_string: String = String::new();
        let mut proteins: Vec<Protein> = Vec::new();

        let file = File::open(file)?;
        
        // Read the lines as bytes, since the input string is not guaranteed to be utf8
        // because of the encoded functional annotations
        let mut lines = ByteLines::new(BufReader::new(file));

        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');

            // uniprot_id, taxon_id and sequence should always contain valid utf8
            let uniprot_id = from_utf8(fields.next().unwrap())?;
            let taxon_id = from_utf8(fields.next().unwrap())?.parse::<TaxonId>()?;
            let sequence = from_utf8(fields.next().unwrap())?;
            let functional_annotations: Vec<u8> = fields.next().unwrap().to_vec();

            if !taxon_aggregator.taxon_exists(taxon_id) {
                continue;
            }

            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());

            proteins.push(Protein {
                uniprot_id: uniprot_id.to_string(),
                taxon_id,
                functional_annotations
            });

        }

        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());
        input_string.shrink_to_fit();
        proteins.shrink_to_fit();
        Ok(Self {
            input_string: input_string.into_bytes(),
            proteins
        })
    }

    /// Creates a `vec<u8>` which represents all the proteins concatenated from the database file
    ///
    /// # Arguments
    /// * `file` - The path to the database file
    /// * `taxon_aggregator` - The `TaxonAggregator` to use
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the `Vec<u8>`
    ///
    /// # Errors
    ///
    /// Returns a `Box<dyn Error>` if an error occurred while reading the database file
    pub fn try_from_database_file_without_annotations(database_file: &str, taxon_aggregator: &TaxonAggregator) -> Result<Vec<u8>, Box<dyn Error>> {
        let mut input_string: String = String::new();

        let file = File::open(database_file)?;

        // Read the lines as bytes, since the input string is not guaranteed to be utf8
        // because of the encoded functional annotations
        let mut lines = ByteLines::new(BufReader::new(file));

        while let Some(Ok(line)) = lines.next() {
            let mut fields = line.split(|b| *b == b'\t');

            // only get the taxon id and sequence from each line, we don't need the other parts
            fields.next();
            let taxon_id = from_utf8(fields.next().unwrap())?.parse::<TaxonId>()?;
            let sequence = from_utf8(fields.next().unwrap())?;
            fields.next();

            if !taxon_aggregator.taxon_exists(taxon_id) {
                continue;
            }

            input_string.push_str(&sequence.to_uppercase());
            input_string.push(SEPARATION_CHARACTER.into());
        }

        input_string.pop();
        input_string.push(TERMINATION_CHARACTER.into());

        input_string.shrink_to_fit();
        Ok(input_string.into_bytes())
    }
    
}

impl Index<usize> for Proteins {
    type Output = Protein;

    fn index(&self, index: usize) -> &Self::Output {
        &self.proteins[index]
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::Write,
        path::PathBuf
    };

    use fa_compression::algorithm1::decode;
    use tempdir::TempDir;

    use super::*;
    use crate::taxonomy::AggregationMethod;

    fn create_database_file(tmp_dir: &TempDir) -> PathBuf {
        let database_file = tmp_dir.path().join("database.tsv");
        let mut file = File::create(&database_file).unwrap();

        file.write("P12345\t1\tMLPGLALLLLAAWTARALEV\t".as_bytes())
            .unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27])
            .unwrap();
        file.write("\n".as_bytes()).unwrap();
        file.write("P54321\t2\tPTDGNAGLLAEPQIAMFCGRLNMHMNVQNG\t".as_bytes())
            .unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27])
            .unwrap();
        file.write("\n".as_bytes()).unwrap();
        file.write("P67890\t6\tKWDSDPSGTKTCIDT\t".as_bytes())
            .unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27])
            .unwrap();
        file.write("\n".as_bytes()).unwrap();
        file.write("P13579\t17\tKEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH\t".as_bytes())
            .unwrap();
        file.write_all(&[0xD1, 0x11, 0xA3, 0x8A, 0xD1, 0x27, 0x47, 0x5E, 0x11, 0x99, 0x27])
            .unwrap();
        file.write("\n".as_bytes()).unwrap();

        database_file
    }

    fn create_taxonomy_file(tmp_dir: &TempDir) -> PathBuf {
        let taxonomy_file = tmp_dir.path().join("taxonomy.tsv");
        let mut file = File::create(&taxonomy_file).unwrap();

        writeln!(file, "1\troot\tno rank\t1\t\x01").unwrap();
        writeln!(file, "2\tBacteria\tsuperkingdom\t1\t\x01").unwrap();
        writeln!(file, "6\tAzorhizobium\tgenus\t1\t\x01").unwrap();
        writeln!(file, "7\tAzorhizobium caulinodans\tspecies\t6\t\x01").unwrap();
        writeln!(file, "9\tBuchnera aphidicola\tspecies\t6\t\x01").unwrap();
        writeln!(file, "10\tCellvibrio\tgenus\t6\t\x01").unwrap();
        writeln!(file, "11\tCellulomonas gilvus\tspecies\t10\t\x01").unwrap();
        writeln!(file, "13\tDictyoglomus\tgenus\t11\t\x01").unwrap();
        writeln!(file, "14\tDictyoglomus thermophilum\tspecies\t10\t\x01").unwrap();
        writeln!(file, "16\tMethylophilus\tgenus\t14\t\x01").unwrap();
        writeln!(file, "17\tMethylophilus methylotrophus\tspecies\t16\t\x01").unwrap();
        writeln!(file, "18\tPelobacter\tgenus\t17\t\x01").unwrap();
        writeln!(file, "19\tSyntrophotalea carbinolica\tspecies\t17\t\x01").unwrap();
        writeln!(file, "20\tPhenylobacterium\tgenus\t19\t\x01").unwrap();

        taxonomy_file
    }

    #[test]
    fn test_new_protein() {
        let protein = Protein {
            uniprot_id:             "P12345".to_string(),
            taxon_id:               1,
            functional_annotations: vec![0xD1, 0x11]
        };

        assert_eq!(protein.uniprot_id, "P12345");
        assert_eq!(protein.taxon_id, 1);
        assert_eq!(protein.functional_annotations, vec![0xD1, 0x11]);
    }

    #[test]
    fn test_new_proteins() {
        let proteins = Proteins {
            input_string: "MLPGLALLLLAAWTARALEV-PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG"
                .as_bytes()
                .to_vec(),
            proteins:     vec![
                Protein {
                    uniprot_id:             "P12345".to_string(),
                    taxon_id:               1,
                    functional_annotations: vec![0xD1, 0x11]
                },
                Protein {
                    uniprot_id:             "P54321".to_string(),
                    taxon_id:               2,
                    functional_annotations: vec![0xD1, 0x11]
                },
            ]
        };

        assert_eq!(
            proteins.input_string,
            "MLPGLALLLLAAWTARALEV-PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG".as_bytes()
        );
        assert_eq!(proteins.proteins.len(), 2);
        assert_eq!(proteins.proteins[0].uniprot_id, "P12345");
        assert_eq!(proteins.proteins[0].taxon_id, 1);
        assert_eq!(proteins.proteins[0].functional_annotations, vec![0xD1, 0x11]);
        assert_eq!(proteins.proteins[1].uniprot_id, "P54321");
        assert_eq!(proteins.proteins[1].taxon_id, 2);
        assert_eq!(proteins.proteins[1].functional_annotations, vec![0xD1, 0x11]);
    }

    #[test]
    fn test_get_taxon() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_taxon").unwrap();

        let database_file = create_database_file(&tmp_dir);
        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
        .unwrap();
        let proteins =
            Proteins::try_from_database_file(database_file.to_str().unwrap(), &taxon_aggregator)
                .unwrap();

        let taxa = vec![1, 2, 6, 17];
        for (i, protein) in proteins.proteins.iter().enumerate() {
            assert_eq!(protein.taxon_id, taxa[i]);
        }
    }

    #[test]
    fn test_get_functional_annotations() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);
        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
        .unwrap();
        let proteins =
            Proteins::try_from_database_file(database_file.to_str().unwrap(), &taxon_aggregator)
                .unwrap();

        for protein in proteins.proteins.iter() {
            assert_eq!(
                decode(&protein.functional_annotations),
                "GO:0009279;IPR:IPR016364;IPR:IPR008816"
            );
        }
    }

    #[test]
    fn test_get_concatenated_proteins() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_get_fa").unwrap();

        let database_file = create_database_file(&tmp_dir);
        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
            .unwrap();
        let proteins =
            Proteins::try_from_database_file_without_annotations(database_file.to_str().unwrap(), &taxon_aggregator)
                .unwrap();
        
        let sep_char = SEPARATION_CHARACTER as char;
        let end_char = TERMINATION_CHARACTER as char;
        let expected = format!("MLPGLALLLLAAWTARALEV{}PTDGNAGLLAEPQIAMFCGRLNMHMNVQNG{}KWDSDPSGTKTCIDT{}KEGILQYCQEVYPELQITNVVEANQPVTIQNWCKRGRKQCKTHPH{}", sep_char, sep_char, sep_char, end_char);
        assert_eq!(proteins, expected.as_bytes());
    }
}
