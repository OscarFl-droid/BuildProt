# Validation record

## Local validation performed before packaging

The packaged repository was tested with Node.js 22 and Python 3.13 in the build environment.

`npm test` completed with **14/14 passing tests**. These tests exercise:

- one exactly unique peptide locus;
- identical peptide sequence in multiple proteins;
- two exact loci within one protein while `Protein_match_count = 1`;
- canonical-only versus canonical-plus-reviewed-isoform universes;
- a one-residue local off-target neighbour with correct mismatch position;
- contact-weighted uniqueness search;
- human-versus-mouse fixture differences;
- the exact score gate and headroom formula;
- BLOSUM62 self-normalisation;
- suffix-array nearest-neighbour results versus an independent exhaustive sliding-window Hamming implementation;
- bare, tab-separated, comma-separated, and Excel-column input parsing;
- duplicate input retention/identification;
- invalid amino-acid character rejection;
- the **actual Web Worker module** loading gzip-compressed sequence, metadata, and suffix-array fixture assets and returning the expected exact/homology result.

Python source compilation (`python -m py_compile scripts/*.py`) passed. A small local FASTA was also parsed by the proteome builder and a fixture suffix array was constructed successfully.

## Full UniProt validation

The execution environment used to package this repository does not provide general outbound network access to `rest.uniprot.org`/the UniProt FTP host, so the full current human and mouse proteome bundles could not be downloaded here. The repository therefore does not pretend that a particular UniProt release was locally validated when it was not.

Instead, `.github/workflows/pages.yml` performs the full network-dependent validation on GitHub Actions during deployment:

1. download the current human and mouse UniProt reference proteomes;
2. build the static indexes and record the actual UniProt release/build date;
3. install NCBI BLAST+;
4. run `scripts/validate_blast.py` independently on both canonical FASTA files;
5. deterministically select one unique and one multi-protein exact peptide from each FASTA;
6. compare raw-FASTA expected protein matches with full-length 100%-identity `blastp -task blastp-short` hits;
7. fail deployment if any discrepancy is observed.

BLAST is used as an independent protein-level exact-match check. It is not used to define locus count because BLAST HSP reporting may merge or suppress repeated loci within one protein. Repeated-locus correctness is tested directly in the core automated tests.

## Expected interpretation of discrepancies

If the GitHub Action fails BLAST validation after a future UniProt release, the deployment is intentionally stopped. The appropriate response is to inspect whether the difference reflects BLAST short-peptide heuristics, a source-format change, accession parsing, or an indexing bug. The workflow must not silently deploy a changed result set until the discrepancy is understood.
