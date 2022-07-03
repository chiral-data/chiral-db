//! ChiralDB Fingerprint
//! 
//! # FingerprintDocument: fingerprint data for a group of molecules
//! 
//! Create a FingerprintDocument from a list of SMILES and use SMILES as data IDs.
//! ```
//! use chiral_db::types;
//! use chiral_db::fingerprint;
//! 
//! let fpk = types::FingerprintKind::OpenBabelECFP4;
//! let nbits: u32 = 2048;
//! let smiles_1 = String::from("c1ccccc1");
//! let smiles_2 = String::from("CCCCCCN");
//! let smiles_vec = vec![ &smiles_1, &smiles_2 ];
//! let ids: Vec<String> = smiles_vec.iter().map(|&smiles| smiles.clone()).collect();
//! let _fp_doc = fingerprint::FingerprintDocument::new_from_smiles_vec(&fpk, nbits, &smiles_vec, &ids);
//! ```

extern crate openbabel;
extern crate chiral_db_sources;
use crate::types;
use crate::utils;
use crate::config;

pub type FingerprintData = Vec<u32>;
pub type FingerprintDB = std::collections::HashMap<String, std::sync::Arc<FingerprintDocument>>;

///
/// Fingerprint Kind
/// 

#[derive(Clone, Debug)]
enum Kind {
    OpenBabel { kind: openbabel::fingerprint::Kind }
}

impl Kind {
    fn as_string(&self) -> String {
        match self {
            Kind::OpenBabel { kind: ob_fp_kind } => {
                match ob_fp_kind {
                    openbabel::fingerprint::Kind::FP2 { nbits } => format!("OpenBabel FP2 with {} bits", nbits),
                    openbabel::fingerprint::Kind::FP3 { nbits } => format!("OpenBabel FP3 with {} bits", nbits),
                    openbabel::fingerprint::Kind::FP4 { nbits } => format!("OpenBabel FP4 with {} bits", nbits),
                    openbabel::fingerprint::Kind::ECFP0 { nbits } => format!("OpenBabel ECFP0 with {} bits", nbits),
                    openbabel::fingerprint::Kind::ECFP2 { nbits } => format!("OpenBabel ECFP2 with {} bits", nbits),
                    openbabel::fingerprint::Kind::ECFP4 { nbits } => format!("OpenBabel ECFP4 with {} bits", nbits),
                    openbabel::fingerprint::Kind::ECFP6 { nbits } => format!("OpenBabel ECFP6 with {} bits", nbits),
                    openbabel::fingerprint::Kind::ECFP8 { nbits } => format!("OpenBabel ECFP8 with {} bits", nbits),
                    openbabel::fingerprint::Kind::ECFP10 { nbits } => format!("OpenBabel ECFP10 with {} bits", nbits),
                }
            }
        }
    }
}

fn kind_openbabel_ecfp0(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP0 { nbits } } }
fn kind_openbabel_ecfp2(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP2 { nbits } } }
fn kind_openbabel_ecfp4(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits } } }
fn kind_openbabel_ecfp6(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP6 { nbits } } }
fn kind_openbabel_ecfp8(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP8 { nbits } } }
fn kind_openbabel_ecfp10(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP10 { nbits } } }


fn convert_kind(fpk: &types::FingerprintKind, nbits: u32) -> Kind {
    match fpk {
        types::FingerprintKind::OpenBabelECFP0 =>  kind_openbabel_ecfp0(nbits),
        types::FingerprintKind::OpenBabelECFP2 =>  kind_openbabel_ecfp2(nbits),
        types::FingerprintKind::OpenBabelECFP4 =>  kind_openbabel_ecfp4(nbits),
        types::FingerprintKind::OpenBabelECFP6 =>  kind_openbabel_ecfp6(nbits),
        types::FingerprintKind::OpenBabelECFP8 =>  kind_openbabel_ecfp8(nbits),
        types::FingerprintKind::OpenBabelECFP10 =>  kind_openbabel_ecfp10(nbits)
    }
}

///
/// Fingerprint Generation
/// 

fn ob_get_fingerprint_fpg_mol(fpg: &openbabel::fingerprint::FingerprintGenerator, mol: &openbabel::molecule::Molecule) -> Vec<u32> {
    let fpd_cxx = fpg.get_fingerprint(&mol);
    utils::cxx_vector_into_vector(&fpd_cxx)
}

pub fn get_fingerprint_for_smiles(fp_doc: &FingerprintDocument, smiles: &String) -> Vec<u32> {
    match &fp_doc.kind {
        Kind::OpenBabel { kind: ob_fp_kind } => {
            let fpg = openbabel::fingerprint::FingerprintGenerator::new(ob_fp_kind.clone());
            let mol = openbabel::molecule::Molecule::new_from_smiles(smiles.as_str());
            ob_get_fingerprint_fpg_mol(&fpg, &mol)
        }
    }
}


///
/// FingerprintDocument
/// 

#[derive(Debug)]
pub struct FingerprintDocument {
    kind: Kind,
    ids: Vec<String>,
    data: FingerprintData,
    nints: usize
}

impl FingerprintDocument {
    pub fn desc(&self) -> String {
        vec![
            format!("{}", self.len()),
            self.kind.as_string()
        ].join("\t\t")
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    // pub fn nbits(&self) -> usize {
    //     match &self.kind {
    //         Kind::OpenBabel { kind: ob_fp_kind } => ob_fp_kind.get_nbits().clone() as usize
    //     }
    // }

    pub fn new_from_smiles_vec(fpk: &types::FingerprintKind, nbits: u32, smiles_vec: &Vec<&String>, ids_input: &Vec<String>) -> Self {
        let kind = convert_kind(fpk, nbits);
        let ids = ids_input.to_vec();
        match &kind {
            Kind::OpenBabel { kind: ob_fp_kind } => {
                let fpg = openbabel::fingerprint::FingerprintGenerator::new(ob_fp_kind.clone());
                let data: FingerprintData = fpg.get_fingerprint_for_smiles_vec(smiles_vec).iter()
                    .map(|fpd_cxx| utils::cxx_vector_into_vector(&fpd_cxx))
                    .flatten()
                    .collect();
                let nints = (nbits / 32) as usize;

                Self { kind, ids, data, nints }
            }
        }
    }

    pub fn new_from_chembl(fpk: &types::FingerprintKind, nbits: u32, filepath: &std::path::Path) -> Self {
        let source_chembl = chiral_db_sources::chembl::SourceChembl::new(filepath);
        let (smiles_vec, ids) = source_chembl.get_smiles_id_pairs();
        FingerprintDocument::new_from_smiles_vec(fpk, nbits, &smiles_vec, &ids)
    }

    pub fn get_data(&self, idx: usize) -> &[u32] {
        &self.data[(idx * self.nints)..(idx * self.nints + self.nints)]
    }

    pub fn get_id(&self, idx: usize) -> &String {
        &self.ids[idx]
    }
}

///
/// Data Loading
/// 

fn load_fp_chembl(fp_doc: &config::FingerprintDocumentConfig) -> FingerprintDocument {
    let filepath = std::path::Path::new(&fp_doc.filepath);
    FingerprintDocument::new_from_chembl(&fp_doc.kind, fp_doc.nbits, &filepath)
}
    
pub fn load_fingerprint_db(conf: &config::Config) -> FingerprintDB {
    let mut fp_doc_dict = FingerprintDB::new();

    if let Some(fp_docs) = &conf.fp_doc {
        for fp_doc in fp_docs.iter() {
            match &fp_doc.source {
                types::Source::Chembl => {
                    fp_doc_dict.insert(fp_doc.name.clone(), std::sync::Arc::new(load_fp_chembl(fp_doc)));
                }
                types::Source::Zinc => {
                    unimplemented!();
                }
            }
        }
    }

    fp_doc_dict
}

#[cfg(test)]
mod test_fp {
    use super::*;

    #[test]
    fn test_kind() {
        let kind = Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits: 2048} };
        assert_eq!(kind.as_string(), "OpenBabel ECFP4 with 2048 bits");
    }


    #[test]
    fn test_new_from_smiles_vec() {
        let fpk = types::FingerprintKind::OpenBabelECFP4;
        let nbits: u32 = 2048;
        let smiles_1 = String::from("c1ccccc1");
        let smiles_2 = String::from("CCCCCCN");
        let smiles_vec = vec![ &smiles_1, &smiles_2 ];
        let ids: Vec<String> = smiles_vec.iter().map(|&smiles| smiles.clone()).collect();
        let fp_doc = FingerprintDocument::new_from_smiles_vec(&fpk, nbits, &smiles_vec, &ids);
        assert_eq!(smiles_vec[0].to_string(), fp_doc.ids[0]);
        assert_eq!(smiles_vec[1].to_string(), fp_doc.ids[1]);
        assert_eq!(fp_doc.data.len(), fp_doc.nints * 2);
    }

    #[test]
    fn test_new_from_chembl() {
        let fpk = types::FingerprintKind::OpenBabelECFP4;
        let nbits: u32 = 2048;
        let filepath = std::path::Path::new("./data/chembl_30_chemreps_100.txt");
        let fp_doc = FingerprintDocument::new_from_chembl(&fpk, nbits, filepath);
        assert_eq!(fp_doc.data.len(), fp_doc.nints * 100);
        assert_eq!(fp_doc.desc(), "100\t\tOpenBabel ECFP4 with 2048 bits");
    }
}