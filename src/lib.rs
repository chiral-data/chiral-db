pub mod types;
mod config;
mod utils;
pub mod fingerprint;
mod similarity;

#[derive(Debug)]
pub struct ChiralDB {
    fp_docs: fingerprint::FingerprintDB 
}

impl ChiralDB {
    pub fn new() -> Self {
        let conf = config::load_configuration();
        let fp_docs = fingerprint::load_fingerprint_db(&conf);

        Self { fp_docs }
    }

    fn desc_fingerprint_db(&self) -> String {
        let mut desc: Vec<String> = vec![vec!["Doc", "Entries", "FP Type"].join("\t\t"), "=".repeat(50)];
        desc.extend(self.fp_docs.keys().map(|key| format!("{}\t\t{}", key, self.fp_docs[key].as_ref().desc())).collect::<Vec<String>>());
        desc.join("\n")
    }

    pub fn desc(&self) -> String {
        vec![
            self.desc_fingerprint_db()
        ].join("\n".repeat(3).as_str())
    }

    pub fn query_similarity_for_smiles(&self, doc_name: &String, smiles: &String, cut_off: f32) -> types::ResultSimilarity {
        if self.fp_docs.contains_key(doc_name) {
            similarity::query_similarity_for_smiles(smiles, &self.fp_docs[doc_name], cut_off)
        } else {
            types::ResultSimilarity::new()
        }
    }
}
