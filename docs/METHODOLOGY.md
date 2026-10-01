# Scientific methodology

## Scope

This application evaluates **proteomic sequence uniqueness and local off-target sequence similarity** for short peptide targets. It is intended to rank candidates before peptide-specific synthetic-binder design. It does **not** estimate peptide abundance, MS detectability, digestion behaviour, chemical tractability, Binder Designability Score, structural accessibility, affinity, or experimentally observed cross-reactivity.

Sequence uniqueness is necessary in many target-selection contexts but is not sufficient to establish binder specificity. A binder may cross-react with a non-identical peptide, particularly when substitutions occur outside the structural interaction hotspot. Every result therefore requires downstream structural and experimental specificity testing.

## Reference proteomes and versioning

The build pipeline uses UniProt reference proteomes. Human is `UP000005640` and mouse is `UP000000589`. Canonical/reference sequences are downloaded from UniProt's current reference-proteome FASTA distribution. The extended universe adds manually reviewed UniProtKB isoform sequences obtained during the same build. The application never sends individual peptide queries to UniProt.

Every generated dataset records the UniProt release reported by the REST response headers, the proteome identifier, build date, number of canonical proteins, and number of added reviewed isoforms. These values are shown in the interface and exported results. Because the build follows UniProt's `current_release`, a later rebuild may legitimately change exact-match or nearest-neighbour results.

UniProt documentation describes proteome identifiers as unique identifiers for sequence sets corresponding to individual genome assemblies and provides canonical as well as canonical-plus-isoform sequence downloads. See: UniProt Consortium, *UniProt: the Universal Protein Knowledgebase in 2025*, **Nucleic Acids Research**; and the UniProt Proteomes help pages.

## Input normalisation

Sequences are converted to uppercase and all whitespace is removed. Only the 20 standard amino-acid symbols `ACDEFGHIKLMNPQRSTVWY` are accepted as query input. Ambiguity symbols (`B`, `J`, `X`, `Z`), selenocysteine (`U`), pyrrolysine (`O`), stop symbols, punctuation, and modification notation are rejected rather than silently interpreted. Duplicate input rows are retained in output but their sequence analysis is cached and reused.

UniProt source proteins may themselves contain ambiguity letters such as `X`; these are retained as positional barriers and cannot exactly match a standard-residue query at that position. Isoleucine and leucine are treated as distinct residues. This is deliberate for proteomic sequence definition. Users who wish to model MS indistinguishability of I/L should implement that as a separate analysis universe because it changes the definition of an exact match.

## Exact proteome uniqueness

For a peptide of length `L`, every exact occurrence in every protein sequence of the selected analysis universe is enumerated. A **proteomic locus** is one start position in one protein/isoform sequence. Repeated occurrences within a single protein are therefore separate loci.

`Exact_match_count` is the number of proteomic loci. `Protein_match_count` is the number of distinct accessions containing at least one locus. Gene symbols and all protein positions are retained separately.

`Proteome_unique = TRUE` only when:

\[
Exact\_match\_count = 1.
\]

Accordingly, a peptide occurring twice within one protein, or once in each of two isoforms of the same gene, is **not** proteome-unique.

The status categories are `UNIQUE`, `MULTIPLE EXACT MATCHES`, `NOT FOUND`, and `INVALID INPUT`.

### Relationship among multiple exact matches

If all exact loci map to one gene symbol, the application reports `same gene / isoforms / repeated locus`. If multiple genes share an identical UniProt protein-family annotation, it reports `probable paralogous proteins (shared UniProt family annotation)`. Otherwise it reports `multiple genes — apparently unrelated or family relationship unannotated`. This last category must not be interpreted as proof that the proteins are evolutionarily unrelated; absence of a shared family annotation is not a phylogenetic analysis.

## Browser search index

The build script concatenates all sequences for a species, separated by a non-amino-acid sentinel, and generates a suffix array with `pydivsufsort`. Sequence, metadata, and suffix-array files are gzip-compressed static assets. They are fetched once by a Web Worker, decompressed locally, and never submitted to a server.

Suffix-array binary search provides exact occurrences of arbitrary query substrings without precomputing indexes for every peptide length.

## Primary off-target homology: exact same-length Hamming search

Whole-protein identity is not used. For a query peptide `Q` of length `L`, the relevant search space is every same-length sliding window in every protein of the selected universe.

For a window `W`:

\[
M(Q,W)=\sum_{i=1}^{L} [Q_i \ne W_i]
\]

\[
Identity(Q,W)=100\times\frac{L-M(Q,W)}{L}.
\]

Mismatch positions are reported using 1-based peptide coordinates.

The sole exact target locus is excluded from the off-target search. If the query has multiple exact loci, only the first locus used as the display target is excluded; the other exact loci remain off-targets with 100% identity and zero mismatches. This makes the exact-uniqueness failure visible in both the exact-match and nearest-neighbour results.

### Exactness of the accelerated nearest-neighbour search

The worker does not heuristically sample the proteome. For a mismatch threshold `d`, the query is split into `d+1` non-overlapping segments. By the pigeonhole principle, any same-length window with at most `d` substitutions must contain at least one segment with zero substitutions. The suffix array enumerates every occurrence of those exact segments; candidate windows are reconstructed and verified residue by residue.

Once at least five verified alternatives exist at distance `d`, every window with smaller or equal Hamming distance has been enumerated, so the five closest identity neighbours are known exactly. For pathological cases that cannot be certified by the seeded search, the implementation falls back to exhaustive same-length scanning in the Web Worker.

## Proteome Unique Score

The score is deliberately independent of any Binder Designability Score.

If `Exact_match_count > 1`:

\[
Proteome\ Unique\ Score = 0
\]

with `FAIL — exact sequence not proteome-unique`.

If `Exact_match_count = 0`, the score is `NA — peptide not found` because species, sequence, or proteome-version errors cannot be excluded.

For an exactly unique peptide:

\[
NearestOffTargetIdentity = \max_{W\ne target} Identity(Q,W)
\]

\[
Proteome\ Unique\ Score = 100\times(1-NearestOffTargetIdentity/100).
\]

The descriptive bands are engineering-prioritisation labels, **not experimentally validated specificity thresholds**: 80–100 very high sequence separation; 60–79 high; 40–59 moderate; 20–39 limited; below 20 very close proteomic neighbour. The raw nearest-neighbour identity and sequence are always shown alongside the score.

## ContactWeightedUniqueness

By default all positions have weight 1. Users may nominate 1-based contact positions; those positions receive a configurable multiplier `c >= 1` while all other positions retain weight 1.

For weights `w_i`:

\[
WeightedIdentity(Q,W)=100\times\frac{\sum_i w_i[Q_i=W_i]}{\sum_i w_i}.
\]

The application searches for the non-target window with the **maximum weighted identity** and reports:

\[
ContactWeightedUniqueness=100-MaxWeightedIdentity.
\]

The search remains exact. After all windows with at most `d` substitutions are enumerated, the minimum possible weighted penalty for any unseen window is bounded by the sum of the `d+1` smallest positional weights. Search continues until that bound proves no unseen window can exceed the current weighted-identity maximum. With no nominated contacts, all weights are equal and this metric reduces to the ordinary sequence-separation score.

## Advanced homology: BLOSUM62

Advanced mode retains ordinary identity/Hamming distance as the principal metric and additionally scores ungapped same-length windows with BLOSUM62.

The raw score is:

\[
S(Q,W)=\sum_i BLOSUM62(Q_i,W_i).
\]

A query-specific normalised score is reported on a 0–100 scale:

\[
S_{norm}=100\times\frac{S-S_{min}}{S_{self}-S_{min}},
\]

where `S_self` is the BLOSUM62 self-alignment score of the query and `S_min` is the sum, over positions, of the minimum BLOSUM62 score possible for the query residue at that position. Values are clipped to 0–100.

The BLOSUM-best alternative is searched independently of the identity-best alternative. Search can stop exactly because, for every query position, the difference between the diagonal score and the best possible non-identical residue gives a minimum score loss. After enumerating all windows through Hamming distance `d`, these losses bound the maximum BLOSUM score any unseen window can attain.

No gapped alignment is used. For short diagnostic peptides this keeps the primary question local, position-preserving, and computationally tractable.

## Validation

Automated fixture tests cover: an exactly unique peptide; multiple proteins with the same peptide; two loci within one protein; a one-residue paralog-like neighbour; an absent peptide; duplicate inputs; invalid characters; canonical-versus-isoform universes; and a mouse-versus-human difference. A separate brute-force sliding-window implementation is compared against suffix-array nearest-neighbour output.

The GitHub Pages workflow additionally installs NCBI BLAST+ and runs `scripts/validate_blast.py` against the freshly downloaded human and mouse UniProt FASTA files. The script deterministically selects one unique and one multi-protein exact peptide, establishes expected matches directly from the raw FASTA, and verifies full-length 100%-identity BLASTP-short protein hits. The workflow fails if a discrepancy is found. BLAST protein-hit validation does not replace explicit locus counting because BLAST may merge repeated HSPs within one protein; repeated-locus behaviour is therefore validated by the direct core tests.

## Limitations

The selected UniProt proteome is an analysis universe, not an exhaustive catalogue of every possible biological peptide. Sequence variants, somatic mutations, unreviewed alternative products not included in the selected universe, proteolytic processing, PTMs, non-reference alleles, microbial peptides, contaminants, and experimental sample composition are outside scope. A peptide absent from the selected proteome is therefore not automatically erroneous, but the application intentionally withholds a normal uniqueness score.

Sequence similarity is not a physical binding model. BLOSUM62 is an evolutionary substitution matrix, not a binder-energy function. Contact weighting is a transparent prioritisation device whose biological usefulness depends on whether the nominated positions truly form the binder interface.
