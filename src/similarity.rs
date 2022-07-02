//! Molecule Similarity
//! 

use crate::fingerprint;
use crate::types;

/// Tanimoto coefficient
pub fn similarity_tanimoto(fpd_1: &[u32], fpd_2: &[u32]) -> f32 {
    let mut andbits: u32 = 0;
    let mut orbits: u32 = 0;
    for i in 0..fpd_1.len() {
        let andfp: u32 = fpd_1[i] & fpd_2[i];
        let orfp: u32 = fpd_1[i] | fpd_2[i];
        andbits += andfp.count_ones();
        orbits += orfp.count_ones();
    }

    if orbits > 0 {
        andbits as f32 / orbits as f32
    } else {
        0.0
    }
}

pub fn query_similarity(fpd_query: &Vec<u32>, fp_doc: &fingerprint::FingerprintDocument, cut_off: f32) -> Vec<(String, f32)> {
    (0..(fp_doc.len()))
        .map(|idx| (
            fp_doc.get_id(idx).to_string(),
            similarity_tanimoto(fpd_query.as_slice(), fp_doc.get_data(idx))
        ))
        .filter(|(_id, tani_eff)| *tani_eff >= cut_off)
        .collect()
}

pub fn query_similarity_for_smiles(smiles: &String, fp_doc: &fingerprint::FingerprintDocument, cut_off: f32) -> types::ResultSimilarity {
    let fpd_query = fingerprint::get_fingerprint_for_smiles(fp_doc, smiles); 
    query_similarity(&fpd_query, fp_doc, cut_off)
        .into_iter()
        .collect() 
}

#[cfg(test)]
mod test_similarity {
    use super::*;

    #[test]
    fn test_tanimoto() {
        let fpk = types::FingerprintKind::OpenBabelECFP4;
        let nbits: u32 = 4096;
        let smiles_1 = String::from("c1ccccc1");
        let smiles_2 = String::from("O=C(C)Oc1ccccc1C(=O)O");
        let smiles_vec = vec![&smiles_1, &smiles_2];
        let ids: Vec<String> = smiles_vec.iter().map(|&smiles| smiles.clone()).collect();
        let fp_doc = fingerprint::FingerprintDocument::new_from_smiles_vec(&fpk, nbits, &smiles_vec, &ids);
        assert!((similarity_tanimoto(fp_doc.get_data(0), fp_doc.get_data(1)) - 0.0666) < 1e4);
    }

    #[test]
    fn test_query_similarity() {
        let fpk = types::FingerprintKind::OpenBabelECFP4;
        let nbits: u32 = 4096;
        let filepath = std::path::Path::new("./data/chembl_30_chemreps_100.txt");
        let fp_doc = fingerprint::FingerprintDocument::new_from_chembl(&fpk, nbits, &filepath);
        let smiles = String::from("Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1");
        let fpd_query = fingerprint::get_fingerprint_for_smiles(&fp_doc, &smiles);
        let res_fpd = query_similarity(&fpd_query, &fp_doc, 1.0);
        assert_eq!(res_fpd.len(), 1);
        assert_eq!(res_fpd[0].0, "CHEMBL263810");
        let res_smiles = query_similarity_for_smiles(&smiles, &fp_doc, 1.0);
        assert_eq!(res_smiles.keys().len(), 1);
        assert!(res_smiles.contains_key("CHEMBL263810"));
    }
}