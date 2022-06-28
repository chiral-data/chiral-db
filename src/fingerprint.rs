//! ChiralDB Fingerprint
//! 
//! # FingerprintDocument: a group of fingerprint data
//! 
//! Create a FingerprintDocument from a list of SMILES and use SMILES as data IDs.
//! ```
//! use chiral_db::fingerprint;
//! 
//! let fpk = fingerprint::Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits: 2048} };
//! let smiles_vec = vec![
//!     String::from("c1ccccc1"),
//!     String::from("CCCCCCN"),
//! ];
//! let _fp_doc = fingerprint::FingerprintDocument::new_from_smiles_vec(&smiles_vec, &fpk, &smiles_vec);
//! ```

extern crate openbabel;
extern crate chiral_db_sources;
use crate::utils;
use crate::common;


fn ob_get_fingerprint_fpg_mol(fpg: &openbabel::fingerprint::FingerprintGenerator, mol: &openbabel::molecule::Molecule) -> Vec<u32> {
    let fpd_cxx = fpg.get_fingerprint(&mol);
    utils::cxx_vector_into_vector(&fpd_cxx)
}

pub fn get_fingerprint_from_smiles(fpk: &Kind, smiles: &String) -> Vec<u32> {
    match fpk {
        Kind::OpenBabel { kind: ob_fp_kind } => {
            let fpg = openbabel::fingerprint::FingerprintGenerator::new(ob_fp_kind.clone());
            let mol = openbabel::molecule::Molecule::new_from_smiles(smiles.as_str());
            ob_get_fingerprint_fpg_mol(&fpg, &mol)
        }
    }
}

#[derive(Clone, Debug)]
pub enum Kind {
    OpenBabel { kind: openbabel::fingerprint::Kind }
}

pub fn kind_openbabel_ecfp0(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP0 { nbits } } }
pub fn kind_openbabel_ecfp2(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP2 { nbits } } }
pub fn kind_openbabel_ecfp4(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits } } }
pub fn kind_openbabel_ecfp6(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP6 { nbits } } }
pub fn kind_openbabel_ecfp8(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP8 { nbits } } }
pub fn kind_openbabel_ecfp10(nbits: u32) -> Kind { Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP10 { nbits } } }

#[derive(Debug)]
pub struct FingerprintDocument {
    kind: Kind,
    ids: Vec<String>,
    data: common::FingerprintData,
    span: usize
}

impl FingerprintDocument {
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn bits_count(&self) -> usize {
        match &self.kind {
            Kind::OpenBabel { kind: ob_fp_kind } => ob_fp_kind.get_nbits().clone() as usize
        }
    }

    pub fn new_from_smiles_vec(smiles_vec: &Vec<String>, fpk: &Kind, ids_input: &Vec<String>) -> Self {
        let mut ids: Vec<String> = vec![];
        let mut data: common::FingerprintData = vec![];
        let kind = fpk.clone();

        match fpk {
            Kind::OpenBabel { kind: ob_fp_kind } => {
                let fpg = openbabel::fingerprint::FingerprintGenerator::new(ob_fp_kind.clone());
                for (smiles, id) in smiles_vec.iter().zip(ids_input.iter()) {
                    let mol = openbabel::molecule::Molecule::new_from_smiles(smiles);
                    let fpd_cxx = fpg.get_fingerprint(&mol);
                    let mut fpd = utils::cxx_vector_into_vector(&fpd_cxx);
                    ids.push(id.clone());
                    data.append(&mut fpd);
                }
                let span: usize = ((*ob_fp_kind.get_nbits()) / 32) as usize;

                Self { kind, ids, data, span }
            }
        }
    }

    pub fn new_from_chembl(fpk: &Kind) -> Self {
        let mut ids: Vec<String> = vec![];
        let mut data: common::FingerprintData = vec![];
        let kind = fpk.clone();
        let mut source_chembl = chiral_db_sources::chembl::SourceChembl::new();
        source_chembl.load();

        match fpk {
            Kind::OpenBabel { kind: ob_fp_kind } => {
                let fpg = openbabel::fingerprint::FingerprintGenerator::new(ob_fp_kind.clone());
                for ec in source_chembl.get_all().values() {
                    let mol = openbabel::molecule::Molecule::new_from_smiles(&ec.smiles.as_str());
                    let mut fpd = ob_get_fingerprint_fpg_mol(&fpg, &mol);
                    ids.push(ec.chembl_id.clone());
                    data.append(&mut fpd);
                }
                let span: usize = ((*ob_fp_kind.get_nbits()) / 32) as usize;

                Self { kind, ids, data, span }
            }
        }
    }

    pub fn get_data(&self, idx: usize) -> &[u32] {
        &self.data[(idx * self.span)..(idx * self.span + self.span)]
    }

    pub fn get_id(&self, idx: usize) -> &String {
        &self.ids[idx]
    }

    pub fn get_kind(&self) -> &Kind {
        &self.kind
    }
}

#[cfg(test)]
mod test_fp {
    use super::*;

    #[test]
    fn test_new_from_smiles_vec() {
        let fpk = Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits: 2048} };
        let smiles_vec = vec![
            String::from("c1ccccc1"),
            String::from("CCCCCCN"),
        ];
        let fp_doc = FingerprintDocument::new_from_smiles_vec(&smiles_vec, &fpk, &smiles_vec);
        assert_eq!(smiles_vec[0], fp_doc.ids[0]);
        assert_eq!(smiles_vec[1], fp_doc.ids[1]);
        assert_eq!(fp_doc.data.len(), fp_doc.bits_count() / 32 * 2);
    }

    #[test]
    fn test_new_from_chembl() {
        let fpk = Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits: 2048} };
        let fp_doc = FingerprintDocument::new_from_chembl(&fpk);
        assert_eq!(fp_doc.data.len(), fp_doc.bits_count() / 32 * 99);
    }
}