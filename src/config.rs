//! Database Configuration 
//! 
//! read configuration from configuration file ChiralDB.toml

#![allow(dead_code)]

use serde_derive::Deserialize;
use crate::types;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub fp_doc: Option<Vec<FingerprintDocumentConfig>>
}

impl Config {}

#[derive(Debug, Deserialize)]
pub struct FingerprintDocumentConfig {
    pub name: String,
    pub kind: types::FingerprintKind, 
    pub nbits: u32,
    pub filepath: String,
    pub source: types::Source 
}

pub fn load_configuration() -> Config {
    let default_filepath = std::path::Path::new("./ChiralDB.toml");
    let contents = std::fs::read_to_string(default_filepath).expect("Error on reading configuration file");
    toml::from_str(contents.as_str()).unwrap()
}


#[cfg(test)]
mod test_config {
    use super::*;

    #[test]
    fn test_fp_doc() {
        let toml_str = r#"
            [[fp_doc]]
            name = "ChEMBL"
            kind = "OpenBabelECFP0"
            nbits = 2048
            filepath = "./chembl.txt"
            source = "Chembl"

            [[fp_doc]]
            name = "ZINC15"
            kind = "OpenBabelECFP4"
            nbits = 1024
            filepath = "./zinc.txt"
            source = "Zinc"
        "#;

        let decoded: Config = toml::from_str(toml_str).unwrap();
        assert!(decoded.fp_doc.is_some());
        let fp_docs = decoded.fp_doc.unwrap();
        assert_eq!(fp_docs[0].name, "ChEMBL");
        assert_eq!(fp_docs[0].kind, types::FingerprintKind::OpenBabelECFP0);
        assert_eq!(fp_docs[0].nbits, 2048);
        assert_eq!(fp_docs[0].filepath, "./chembl.txt");
        assert_eq!(fp_docs[0].source, types::Source::Chembl);
        assert_eq!(fp_docs[1].name, "ZINC15");
        assert_eq!(fp_docs[1].kind, types::FingerprintKind::OpenBabelECFP4);
        assert_eq!(fp_docs[1].nbits, 1024);
        assert_eq!(fp_docs[1].filepath, "./zinc.txt");
        assert_eq!(fp_docs[1].source, types::Source::Zinc);
    }
}

