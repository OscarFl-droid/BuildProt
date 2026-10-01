#!/usr/bin/env python3
"""Build static browser-searchable UniProt reference-proteome bundles.

The deployed application never queries UniProt per peptide. This script is the only
network-dependent stage. It downloads canonical reference-proteome FASTA files,
adds manually reviewed UniProtKB isoforms, constructs one suffix array per species,
and writes compressed static data files plus a versioned manifest.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, os, re, sys, time
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
    s.headers.update({"User-Agent": "peptide-proteome-uniqueness/1.0 (+GitHub Pages static build)"})
    adapter = requests.adapters.HTTPAdapter(max_retries=4)
    s.mount("https://", adapter)
    return s


def get_bytes(s: requests.Session, url: str, params=None, timeout=180) -> tuple[bytes, requests.Response]:
    r = s.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.content, r


def maybe_gunzip(data: bytes) -> bytes:
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


def parse_fasta(raw: bytes, is_canonical: bool) -> list[Protein]:
    text = maybe_gunzip(raw).decode("utf-8")
    records: list[Protein] = []
    header = None
    chunks: list[str] = []

    def flush():
        nonlocal header, chunks
        if header is None:
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
        return obj["metadata"], obj.get("release"), obj.get("release_date")
    params = {
        "query": f"proteome:{proteome_id}",
        "format": "tsv",
        "fields": "accession,gene_primary,protein_name,protein_families",
        "compressed": "true",
    }
    raw, resp = get_bytes(s, "https://rest.uniprot.org/uniprotkb/stream", params=params)
    text = maybe_gunzip(raw).decode("utf-8")
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    meta = {}
    for row in reader:
        acc = row.get("Entry", "").strip()
        if not acc: continue
        meta[acc] = {
            "gene": row.get("Gene Names (primary)", "").strip(),
            "protein": row.get("Protein names", "").strip(),
            "family": row.get("Protein families", "").strip(),
        }
    release = resp.headers.get("x-uniprot-release")
    release_date = resp.headers.get("x-uniprot-release-date")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_bytes(gzip.compress(json.dumps({"metadata":meta,"release":release,"release_date":release_date}).encode(), compresslevel=6))
    return meta, release, release_date


def download_source_fastas(s: requests.Session, key: str, cfg: dict, cache_dir: Path, refresh: bool) -> tuple[bytes, bytes, str|None, str|None]:
    pid, taxid = cfg["proteome_id"], cfg["taxid"]
    canonical_cache = cache_dir / f"{key}.canonical.fasta.gz"
    iso_cache = cache_dir / f"{key}.reviewed_isoforms.fasta.gz"
    base = f"https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{pid}/{pid}_{taxid}.fasta.gz"
    if refresh or not canonical_cache.exists():
        raw, _ = get_bytes(s, base)
        canonical_cache.write_bytes(raw)
    canonical_raw = canonical_cache.read_bytes()

    params = {
        "query": f"(proteome:{pid} AND reviewed:true)",
        "format": "fasta",
        "includeIsoform": "true",
        "compressed": "true",
    }
    if refresh or not iso_cache.exists():
        raw, iso_resp = get_bytes(s, "https://rest.uniprot.org/uniprotkb/stream", params=params)
        iso_cache.write_bytes(raw)
        release = iso_resp.headers.get("x-uniprot-release")
        release_date = iso_resp.headers.get("x-uniprot-release-date")
    else:
        release = release_date = None
    return canonical_raw, iso_cache.read_bytes(), release, release_date


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
    reviewed_all = parse_fasta(iso_raw, True)
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
            "reviewed_isoform_fasta_decompressed": hashlib.sha256(maybe_gunzip(iso_raw)).hexdigest()
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
