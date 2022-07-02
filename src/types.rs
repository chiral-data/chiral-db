//! Shared Types Definition
//! 

use serde_derive::Deserialize;

type ID = String;


#[derive(Debug, Deserialize, PartialEq)]
pub enum Source {
    Chembl,
    Zinc
}


#[derive(Debug, Deserialize, PartialEq)]
pub enum FingerprintKind {
    OpenBabelECFP0,
    OpenBabelECFP2,
    OpenBabelECFP4,
    OpenBabelECFP6,
    OpenBabelECFP8,
    OpenBabelECFP10,
}

type SimilarityScore = f32;
pub type ResultSimilarity = std::collections::HashMap<ID, SimilarityScore>;