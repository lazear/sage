use sage_core::database::{Builder, IndexedDatabase, PeptideIx};
use sage_core::mass::Tolerance;

const BUCKET_SIZE: usize = 23;

fn mk_database(path: &str) -> IndexedDatabase {
    let builder = Builder {
        bucket_size: Some(BUCKET_SIZE),
        fasta: Some(path.into()),
        generate_decoys: Some(false),
        ..Default::default()
    };

    builder.make_parameters().build().unwrap()
}

#[test]
fn visit_all_ions() {
    let database = mk_database("../../tests/Q99536.fasta");

    // Map PeptideIx -> number of fragments between 500 & 700 m/z
    // We want to make sure that IndexedDatabase::query hits *all* of them
    let mut expected = vec![0usize; database.peptides.len()];

    for (chunk_idx, chunk) in database.fragments.chunks(database.bucket_size).enumerate() {
        // Check for total ordering by PeptideIx within a chunk
        let mut last = PeptideIx(0);
        for frag in chunk {
            assert!(frag.peptide_index >= last);
            assert!(frag.fragment_mz >= database.buckets()[chunk_idx]);
            if chunk_idx + 1 < database.buckets().len() {
                assert!(frag.fragment_mz <= database.buckets()[chunk_idx + 1]);
            }

            if frag.fragment_mz >= 500.0 && frag.fragment_mz <= 700.0 {
                expected[frag.peptide_index.0 as usize] += 1;
            }
            last = frag.peptide_index;
        }
    }

    let mut visited = vec![0usize; database.peptides.len()];

    // Hit all peptides in database, track how many of the 500-700 fragment m/z's
    // are returned to us by searching the database.
    let query = database.query(
        1000.0,
        Tolerance::Da(-5000.0, 5000.0),
        Tolerance::Da(-100.0, 100.0),
        0,
        0,
    );

    for fragment in query.page_search(600.0) {
        visited[fragment.peptide_index.0 as usize] += 1;
    }

    // Make sure we visited every possible fragment
    assert_eq!(expected, visited);
}
