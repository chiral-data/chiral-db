//! Molecule Similarity
//! 

use crate::fingerprint;

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

pub fn query_similarity(fpd_query: &Vec<u32>, fp_doc: &fingerprint::FingerprintDocument, cut_off: f32) -> Vec<(f32, String)> {
    (0..(fp_doc.len()))
        .map(|idx| (similarity_tanimoto(fpd_query.as_slice(), fp_doc.get_data(idx)), fp_doc.get_id(idx).to_string()))
        .filter(|(tani_eff, _id)| *tani_eff >= cut_off)
        .collect()
}

#[cfg(test)]
mod test_similarity {
    use super::*;

    #[test]
    fn test_tanimoto() {
        let fpk = fingerprint::Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits: 4096 } };
        let smiles_vec = vec![
            String::from("c1ccccc1"),
            String::from("O=C(C)Oc1ccccc1C(=O)O"),
        ];

        let fp_doc = fingerprint::FingerprintDocument::new_from_smiles_vec(&smiles_vec, &fpk, &smiles_vec);
        assert!((similarity_tanimoto(fp_doc.get_data(0), fp_doc.get_data(1)) - 0.0666) < 1e4);
    }

    #[test]
    fn test_query_similarity() {
        let fpk = fingerprint::Kind::OpenBabel { kind: openbabel::fingerprint::Kind::ECFP4 { nbits: 2048} };
        let mut sc = chiral_db_sources::chembl::SourceChembl::new();
        sc.load();
        let fp_doc = fingerprint::FingerprintDocument::new_from_chembl(&sc, &fpk);
        let fpd_query = fingerprint::get_fingerprint_from_smiles(&fpk, &String::from("Cc1cc(NC(=O)c2cc(Cl)cc(Cl)c2O)ccc1Sc1nc2ccccc2s1"));
        let res = query_similarity(&fpd_query, &fp_doc, 1.0);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].1, "CHEMBL263810");
    }
}