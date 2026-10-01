# Peptide Proteome Uniqueness & Homology

A complete static GitHub Pages application for screening mass-spectrometry-derived peptide sequences before peptide-specific synthetic-binder development. It reports exact proteomic loci, protein/gene mappings, exact proteome uniqueness, the closest same-length proteomic windows, mismatch positions, optional BLOSUM62 similarity, a transparent Proteome Unique Score, and an optional contact-weighted uniqueness metric.

The browser does **not** query a live API for each peptide. Human and mouse UniProt reference-proteome data are downloaded and indexed during the build/deployment workflow, then served as static compressed assets. Peptide analyses run locally in a Web Worker.

## Scientific scope

The tool evaluates **proteomic sequence uniqueness and local off-target sequence similarity**. It does not evaluate peptide abundance, MS quality, tryptic suitability, chemistry, structural accessibility, affinity, or Binder Designability Score. Sequence uniqueness does not guarantee binder specificity; non-identical peptides may cross-react and all candidates require downstream structural design and experimental specificity testing.

Full definitions, formulae, search guarantees, and limitations are in [`docs/METHODOLOGY.md`](docs/METHODOLOGY.md). The packaging-time and deployment validation plan is recorded in [`docs/VALIDATION.md`](docs/VALIDATION.md).

## Supported species

The initial configuration includes:

- Homo sapiens — UniProt reference proteome `UP000005640`
- Mus musculus — UniProt reference proteome `UP000000589`

The species registry is at the top of `scripts/build_proteomes.py`. Additional reference proteomes require only a label, NCBI taxon ID, and UniProt proteome ID.

## Proteome universes

**Canonical/reference proteins only** uses the UniProt reference-proteome canonical FASTA.

**Canonical + reviewed UniProt isoforms** uses the same canonical reference proteome plus manually reviewed UniProtKB isoform sequences fetched during the build. The exact UniProt release, proteome ID, build date, canonical count, and reviewed-isoform count are written to `data/manifest.json`, displayed in the UI, and exported with every result.

The build uses UniProt's authoritative reference-proteome distribution and UniProt REST service. UniProt documents that proteome identifiers uniquely identify the protein set corresponding to a particular genome assembly and that canonical and canonical-plus-isoform sequence sets are available. Reference: UniProt Consortium and UniProt Proteomes documentation.

## Exact uniqueness definition

A peptide is `Proteome_unique = TRUE` **only when its complete sequence occurs at exactly one proteomic locus in the selected analysis universe**. A locus is one start position in one protein/isoform. Two occurrences inside the same protein therefore count as two loci and fail uniqueness even though `Protein_match_count = 1`.

The application reports every accession, gene, protein name, and 1-based position for exact matches. Multiple matches are described as same-gene/isoform/repeated-locus, probable paralogous proteins when a shared UniProt family annotation supports that description, or multiple genes with relationship undetermined/apparently unrelated when no shared family annotation is available.

## Homology algorithm

The principal off-target metric is local, same-length, ungapped peptide identity. Every window of the query length is conceptually considered; whole-protein identity is not used.

A suffix array accelerates the search without making it approximate. At Hamming threshold `d`, the query is split into `d+1` disjoint exact seeds. Every window with at most `d` substitutions must contain at least one intact seed, so suffix-array seed lookup followed by full candidate verification recovers all windows through that distance. Search stops only when the top-five nearest identity neighbours are mathematically certified. A brute-force fallback handles pathological edge cases.

Advanced mode scores local alternatives with BLOSUM62 and independently certifies the highest BLOSUM-similarity alternative using an upper-bound argument described in the methodology.

## Proteome Unique Score

For more than one exact locus, the score is `0` and the result is flagged `FAIL — exact sequence not proteome-unique`. A peptide absent from the selected proteome receives `NA` rather than a normal score.

For one exact locus:

```text
Proteome Unique Score = 100 × (1 - NearestOffTargetIdentity / 100)
```

The displayed bands are engineering-prioritisation descriptors, not experimentally validated specificity thresholds:

```text
80–100  very high sequence separation
60–79   high sequence separation
40–59   moderate sequence separation
20–39   limited separation
<20     very close proteomic neighbour
```

The nearest off-target sequence, identity, mismatches, mismatch positions, accession, gene, and position remain visible; the score never replaces the raw evidence.

## Repository structure

```text
index.html                    application shell
css/app.css                   scientific UI styling
js/app.js                     input, results table, export, UI state
js/parser.js                  peptide/ID parsing and validation
js/search-core.js             suffix-array exact and homology algorithms
js/scoring.js                 BLOSUM62 and uniqueness scoring
js/worker.js                  background proteome loading/search
scripts/build_proteomes.py    UniProt download/update and index builder
scripts/validate_blast.py     independent BLASTP-short spot validation
tests/*.test.mjs              deterministic automated tests
docs/METHODOLOGY.md           scientific methods, formulae, limitations
.github/workflows/pages.yml   tests, proteome update, validation, Pages deploy
data/manifest.json            generated dataset registry/version metadata
```

## Local tests

Node.js 20+ is sufficient for the deterministic unit tests:

```bash
npm test
```

The tests include an independent brute-force same-length search and cover unique, multi-match, repeated-locus, one-residue-neighbour, absent, duplicate, invalid-character, isoform-universe, and human-versus-mouse cases.

Python source can be syntax-checked with:

```bash
python -m py_compile scripts/*.py
```

## Build the full proteome data locally

Use Python 3.11+:

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python scripts/build_proteomes.py --species all --output data --cache .cache/uniprot --refresh
```

`pydivsufsort` constructs the suffix array efficiently. The generated `.gz` files are intentionally git-ignored because the normal GitHub Pages workflow builds them from the current UniProt release during deployment.

Then serve the repository through HTTP; module workers do not work reliably from `file://` URLs:

```bash
python -m http.server 8000
```

Open `http://localhost:8000`.

## Independent BLAST validation

With NCBI BLAST+ installed and the raw FASTA retained in `.cache/uniprot`:

```bash
python scripts/validate_blast.py --fasta .cache/uniprot/human.canonical.fasta.gz
python scripts/validate_blast.py --fasta .cache/uniprot/mouse.canonical.fasta.gz
```

The GitHub deployment workflow runs these automatically after each fresh proteome build and fails on discrepancies.

## Deploy to GitHub Pages — step by step

### 1. Create the repository

On GitHub, create a new public or private repository, for example `peptide-proteome-uniqueness`. Do not initialise it with another README if you are pushing this folder as-is.

### 2. Push the files

From the repository folder:

```bash
git init
git add .
git commit -m "Initial peptide uniqueness application"
git branch -M main
git remote add origin https://github.com/YOUR-USER/peptide-proteome-uniqueness.git
git push -u origin main
```

### 3. Enable GitHub Pages

In the GitHub repository open **Settings → Pages**. Under **Build and deployment**, choose **GitHub Actions** as the source. The included `.github/workflows/pages.yml` workflow will test the code, download the current human and mouse UniProt proteomes, build the static indexes, run independent BLASTP-short validation, and deploy the site.

You can also open **Actions → Build proteomes, validate, and deploy Pages → Run workflow** to trigger a build manually. The workflow is scheduled monthly so the deployed reference proteomes can track UniProt releases. Every deployed result records its actual release and build date; it never silently assumes a version.

### 4. Run/update the proteome build

For GitHub Pages, pushing to `main` or manually running the workflow is sufficient. Locally, use the `build_proteomes.py` command shown above. Add `--refresh` when you explicitly want to redownload source data instead of using `.cache/uniprot`.

### 5. Add another species later

Add an entry to `SPECIES` in `scripts/build_proteomes.py`:

```python
"rat": {
    "label": "Rattus norvegicus",
    "taxid": 10116,
    "proteome_id": "UP000002494"
}
```

Rebuild and redeploy. The command-line species choices are generated from the registry, and no browser-search code changes are required because the UI is populated from `data/manifest.json`.

For a new species, verify the UniProt reference-proteome ID and taxon ID before building. The build script records the selected IDs in the manifest so the analysis universe remains auditable.

## Input formats

All of the following are accepted:

```text
DVFQELIAPK
pep_001<TAB>DVFQELIAPK
pep_002,DVFQELIAPK
```

A single Excel `Stripped.Sequence` column can be pasted directly. Input is uppercased and whitespace is removed. Non-standard amino-acid symbols are flagged rather than silently coerced. Duplicate inputs remain as duplicate output rows but sequence analysis is cached.

## Output and export

The sortable/filterable table contains the requested core fields plus `ContactWeightedUniqueness`. Clicking a row expands all exact loci and at least the five closest off-target windows when five are available. Mismatch positions and a residue-level alignment are shown.

Exports are available as CSV, TSV, UTF-8-BOM Excel-compatible CSV, and clipboard TSV suitable for direct paste back into an Excel peptide workflow.

## Performance architecture

The UI remains responsive because proteome loading and sequence analysis run in a Web Worker. Duplicate peptide sequences are cached. Static proteome bundles are gzip-compressed. Exact matching uses suffix-array binary search; nearest-neighbour calculations use mathematically complete seed-and-verify search rather than repeated DOM work or live remote queries.

Short peptides with extremely weak similarity to the proteome can require wider searches and are intrinsically more expensive. The worker prevents these cases from freezing the page, and an exhaustive fallback preserves correctness.

## Open-source dependencies

The browser application is dependency-free vanilla HTML/CSS/JavaScript. The build pipeline uses `requests`, `numpy`, and `pydivsufsort`; independent validation uses the open-source NCBI BLAST+ command-line tools. GitHub Actions and GitHub Pages provide the build/deployment infrastructure. No paid API, commercial service, or backend is required.

## Data provenance and licensing

Protein data originate from UniProt. Cite UniProt in downstream scientific use and comply with UniProt's current database licence/attribution terms. The application code is provided under the included MIT licence; UniProt-derived sequence data retain the terms specified by UniProt.
