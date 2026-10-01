#!/usr/bin/env python3
"""Build static browser-searchable UniProt reference-proteome bundles.

The deployed application never queries UniProt per peptide. This script is the only
network-dependent stage. It downloads canonical reference-proteome FASTA files,
adds manually reviewed UniProtKB isoforms from the reference-proteome additional FASTA,
constructs one suffix array per species,
and writes compressed static data files plus a versioned manifest.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, os, re, sys, time, zlib
from urllib3.util.retry import Retry
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import numpy as np
import requests

SPECIES = {
    "human": {"label": "Homo sapiens", "taxid": 9606, "proteome_id": "UP000005640"},
    "mouse": {"label": "Mus musculus", "taxid": 10090, "proteome_id": "UP000000589"},
}
PROTEOME_SEQ_RE = re.compile(r"^[A-Z]+$")

@dataclass
class Protein:
    accession: str
    gene: str
    protein: str
    family: str
    sequence: str
    isCanonical: bool
    parentAccession: str
    start: int = 0
    length: int = 0


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": "peptide-proteome-uniqueness/1.1 (+GitHub Pages static build)"})
    # urllib3 covers connection/status failures; get_bytes below additionally retries
    # failures raised when the remote server aborts part-way through a response body.
    retry = Retry(total=2, connect=2, read=2, status=2, backoff_factor=1,
                  status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=frozenset({"GET"}), respect_retry_after_header=True)
    adapter = requests.adapters.HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    return s


def get_bytes(s: requests.Session, url: str, params=None, timeout=(25, 180), attempts=5) -> tuple[bytes, requests.Response]:
    """Retry incomplete bodies as well as request/header failures. No partial cache writes.

    A successful requests.get() does *not* guarantee that .content can finish: UniProt
    /stream previously raised ChunkedEncodingError only during body consumption.
    Loading bytes inside this try-block ensures that case is also retried. For gzipped
    data, the gzip trailer (CRC32 + uncompressed size) is checked before acceptance.
    """
    for attempt in range(1, attempts + 1):
        try:
            with s.get(url, params=params, timeout=timeout, stream=True) as r:
                r.raise_for_status()
                buffer = io.BytesIO()
                for chunk in r.iter_content(chunk_size=256 * 1024):
                    if chunk:
                        buffer.write(chunk)
                raw = buffer.getvalue()
                if not raw:
                    raise ValueError("Received empty response")
                if raw[:2] == b"\x1f\x8b":
                    gzip.decompress(raw)  # Detect an incomplete/corrupt gzip transfer.
                return raw, r
        except (requests.exceptions.RequestException, EOFError, OSError, ValueError, zlib.error) as exc:
            print(f"  Download attempt {attempt}/{attempts} failed: {type(exc).__name__}: {exc}", flush=True)
            if attempt == attempts:
                raise RuntimeError(f"UniProt download failed after {attempts} attempts: {url}") from exc
            time.sleep(min(2 ** attempt, 20))
    raise AssertionError("Unreachable retry path")


def maybe_gunzip(data: bytes) -> bytes:
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


def parse_fasta(raw: bytes, is_canonical: bool, reviewed_only: bool = False) -> list[Protein]:
    text = maybe_gunzip(raw).decode("utf-8")
    records: list[Protein] = []
    header = None
    chunks: list[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
            return
        # Reference-proteome *_additional.fasta.gz contains additional entries for
        # both Swiss-Prot (sp) and TrEMBL (tr). The extended selection intentionally
        # includes *reviewed* additional isoforms only, never all unreviewed variants.
        if reviewed_only and not header.startswith("sp|"):
            return
        seq = "".join(chunks).replace(" ", "").upper()
        if not PROTEOME_SEQ_RE.fullmatch(seq):
            bad = sorted({c for c in seq if not ("A" <= c <= "Z")})
            raise ValueError(f"Non-letter residue(s) {bad} in {header[:100]}")
        parts = header.split("|", 2)
        accession = parts[1] if len(parts) >= 3 else header.split()[0]
        tail = parts[2] if len(parts) >= 3 else header
        desc = tail.split(" ", 1)[1] if " " in tail else tail
        protein_name = desc.split(" OS=", 1)[0].strip()
        gm = re.search(r"\sGN=([^\s]+)", " " + desc)
        gene = gm.group(1) if gm else ""
        records.append(Protein(accession, gene, protein_name, "", seq, is_canonical, accession.split("-",1)[0]))

    for line in text.splitlines():
        if line.startswith(">"):
            flush(); header = line[1:]; chunks = []
        else:
            chunks.append(line.strip())
    flush()
    return records


def fetch_metadata(s: requests.Session, proteome_id: str, cache_path: Path, refresh: bool) -> tuple[dict[str, dict], str|None, str|None]:
    if cache_path.exists() and not refresh:
        obj = json.loads(gzip.decompress(cache_path.read_bytes()).decode())
        if obj.get("metadata"):
            return obj["metadata"], obj.get("release"), obj.get("release_date")

    # Use the resumable-in-practice paginated /search endpoint instead of a single
    # potentially enormous /stream. Any interrupted page is individually retried.
    params = {
        "query": f"proteome:{proteome_id}",
        "format": "tsv",
        "fields": "accession,gene_primary,protein_name,protein_families",
        "size": "500",
    }
    url = "https://rest.uniprot.org/uniprotkb/search"
    meta: dict[str, dict] = {}
    release = release_date = None
    page = 0
    total_expected = None
    seen_urls = set()
    while url:
        if url in seen_urls:
            raise RuntimeError(f"Repeated UniProt pagination URL: {url}")
        seen_urls.add(url)
        raw, resp = get_bytes(s, url, params=params)
        text = maybe_gunzip(raw).decode("utf-8")
        reader = csv.DictReader(io.StringIO(text), delimiter="\t")
        required = {"Entry", "Gene Names (primary)", "Protein names", "Protein families"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise RuntimeError(f"Unexpected UniProt metadata headers: {reader.fieldnames}")
        page_count = 0
        for row in reader:
            acc = (row.get("Entry") or "").strip()
            if not acc:
                continue
            meta[acc] = {
                "gene": (row.get("Gene Names (primary)") or "").strip(),
                "protein": (row.get("Protein names") or "").strip(),
                "family": (row.get("Protein families") or "").strip(),
            }
            page_count += 1
        if page_count == 0:
            raise RuntimeError("Unexpected empty UniProt metadata page")
        this_release = resp.headers.get("x-uniprot-release")
        this_date = resp.headers.get("x-uniprot-release-date")
        if release and this_release and this_release != release:
            raise RuntimeError("UniProt release changed during pagination; rerun the build")
        release = release or this_release
        release_date = release_date or this_date
        if total_expected is None:
            total_header = resp.headers.get("x-total-results")
            total_expected = int(total_header) if total_header else None
        match = re.search(r'<([^>]+)>;\s*rel="next"', resp.headers.get("Link", ""))
        url = match.group(1) if match else None
        params = None  # Next URL already includes full pagination parameters.
        page += 1
        if page % 10 == 0 or not url:
            print(f"  Metadata: {len(meta):,} entries, {page} pages", flush=True)
        if page > 10000:
            raise RuntimeError("Unexpectedly many UniProt metadata pages")
    if total_expected is not None and len(meta) != total_expected:
        raise RuntimeError(f"Metadata incomplete: got {len(meta)}, expected {total_expected}")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    # Write the cache only after the complete paginated dataset has been checked.
    cache_path.write_bytes(gzip.compress(json.dumps({"metadata":meta,"release":release,"release_date":release_date}).encode(), compresslevel=6))
    return meta, release, release_date


def download_source_fastas(s: requests.Session, key: str, cfg: dict, cache_dir: Path, refresh: bool) -> tuple[bytes, bytes, str|None, str|None]:
    pid, taxid = cfg["proteome_id"], cfg["taxid"]
    canonical_cache = cache_dir / f"{key}.canonical.fasta.gz"
    # The official reference-proteome additional FASTA avoids the large REST
    # isoform /stream endpoint (the source of the original ChunkedEncodingError).
    additional_cache = cache_dir / f"{key}.additional.fasta.gz"
    root = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{pid}/{pid}_{taxid}"

    def fetch_fasta(destination: Path, url: str) -> bytes:
        if destination.exists() and not refresh:
            try:
                raw = destination.read_bytes()
                if maybe_gunzip(raw).lstrip().startswith(b">"):
                    return raw
                raise ValueError("Missing FASTA header")
            except (EOFError, OSError, ValueError):
                print(f"  Invalid cached FASTA: {destination}; downloading again", flush=True)
        raw, _ = get_bytes(s, url)
        if not maybe_gunzip(raw).lstrip().startswith(b">"):
            raise RuntimeError(f"UniProt returned non-FASTA content: {url}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(destination.name + ".part")
        temporary.write_bytes(raw)
        temporary.replace(destination)  # Never retain an incomplete source as valid cache.
        print(f"  Downloaded {destination.name}: {len(raw):,} bytes", flush=True)
        return raw

    canonical_raw = fetch_fasta(canonical_cache, root + ".fasta.gz")
    additional_raw = fetch_fasta(additional_cache, root + "_additional.fasta.gz")
    return canonical_raw, additional_raw, None, None


def build_suffix_array(data: bytes) -> np.ndarray:
    try:
        import pydivsufsort
        sa = pydivsufsort.divsufsort(data)
        return np.asarray(sa, dtype="<u4")
    except ImportError:
        if len(data) > 200_000:
            raise RuntimeError("pydivsufsort is required for full proteomes; install requirements.txt")
        order = sorted(range(len(data)), key=lambda i: data[i:])
        return np.asarray(order, dtype="<u4")


def enrich(records: list[Protein], metadata: dict[str, dict]) -> None:
    for p in records:
        m = metadata.get(p.parentAccession) or metadata.get(p.accession)
        if not m: continue
        p.gene = m.get("gene") or p.gene
        p.protein = m.get("protein") or p.protein
        p.family = m.get("family") or p.family


def build_one(key: str, cfg: dict, output: Path, cache: Path, refresh: bool, s: requests.Session) -> dict:
    print(f"[{key}] downloading/loading UniProt source data", flush=True)
    canonical_raw, iso_raw, iso_release, iso_release_date = download_source_fastas(s, key, cfg, cache, refresh)
    metadata, meta_release, meta_release_date = fetch_metadata(s, cfg["proteome_id"], cache / f"{key}.metadata.json.gz", refresh)
    canonical = parse_fasta(canonical_raw, True)
    # The additional FASTA can contain reviewed and unreviewed isoforms/variants.
    # Select sp| entries with isoform accessions only; do not add tr| entries.
    reviewed_all = parse_fasta(iso_raw, False, reviewed_only=True)
    isoforms = [p for p in reviewed_all if "-" in p.accession]
    # Keep only one copy of each accession. Canonical sequence always wins for canonical accession.
    seen = {p.accession for p in canonical}
    isoforms = [p for p in isoforms if p.accession not in seen]
    for p in isoforms: p.isCanonical = False
    records = canonical + isoforms
    enrich(records, metadata)

    pieces: list[str] = []
    cursor = 0
    for p in records:
        if pieces:
            pieces.append("|"); cursor += 1
        p.start = cursor; p.length = len(p.sequence)
        pieces.append(p.sequence); cursor += p.length
    concat = "".join(pieces)
    if len(concat) >= 2**32:
        raise RuntimeError("Proteome exceeds 32-bit suffix-array address space")
    print(f"[{key}] {len(canonical):,} canonical + {len(isoforms):,} reviewed isoforms; {len(concat):,} symbols", flush=True)
    print(f"[{key}] constructing suffix array", flush=True)
    sa = build_suffix_array(concat.encode("ascii"))

    output.mkdir(parents=True, exist_ok=True)
    sequence_name = f"{key}.sequence.txt.gz"
    sa_name = f"{key}.suffix_array.u32.bin.gz"
    metadata_name = f"{key}.metadata.json.gz"
    (output / sequence_name).write_bytes(gzip.compress(concat.encode("ascii"), compresslevel=6))
    (output / sa_name).write_bytes(gzip.compress(sa.tobytes(order="C"), compresslevel=6))
    summary = {
        "species": cfg["label"], "taxid": cfg["taxid"], "proteome_id": cfg["proteome_id"],
        "canonical_proteins": len(canonical), "reviewed_isoforms": len(isoforms), "concat_symbols": len(concat)
    }
    meta_out = {
        "schema_version": 1,
        "summary": summary,
        "proteins": [{k:v for k,v in asdict(p).items() if k != "sequence"} for p in records]
    }
    (output / metadata_name).write_bytes(gzip.compress(json.dumps(meta_out, separators=(",",":")).encode(), compresslevel=6))
    release = meta_release or iso_release or "unknown"
    release_date = meta_release_date or iso_release_date
    build_date = datetime.now(timezone.utc).date().isoformat()
    return {
        "id": key,
        "label": cfg["label"],
        "taxid": cfg["taxid"],
        "proteome_id": cfg["proteome_id"],
        "uniprot_release": release,
        "uniprot_release_date": release_date,
        "build_date": build_date,
        "canonical_proteins": len(canonical),
        "reviewed_isoforms": len(isoforms),
        "universe_definition": {
            "canonical": "UniProt reference-proteome canonical FASTA",
            "extended": "Canonical reference proteome plus manually reviewed UniProtKB isoform sequences"
        },
        "source_sha256": {
            "canonical_fasta_decompressed": hashlib.sha256(maybe_gunzip(canonical_raw)).hexdigest(),
            "additional_reference_fasta_decompressed": hashlib.sha256(maybe_gunzip(iso_raw)).hexdigest()
        },
        "files": {
            "sequence": f"./data/{sequence_name}",
            "suffix_array": f"./data/{sa_name}",
            "metadata": f"./data/{metadata_name}"
        }
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", choices=[*SPECIES.keys(), "all"], default="all")
    ap.add_argument("--output", type=Path, default=Path("data"))
    ap.add_argument("--cache", type=Path, default=Path(".cache/uniprot"))
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    keys = list(SPECIES) if args.species == "all" else [args.species]
    s = session()
    datasets = [build_one(k, SPECIES[k], args.output, args.cache, args.refresh, s) for k in keys]
    manifest = {
        "schema_version": 1,
        "generated": True,
        "build_date": datetime.now(timezone.utc).isoformat(),
        "uniprot_release": "; ".join(sorted(set(d["uniprot_release"] for d in datasets))),
        "datasets": datasets
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {args.output/'manifest.json'}", flush=True)

if __name__ == "__main__": main()
