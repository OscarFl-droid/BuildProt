#!/usr/bin/env python3
"""Independent exact-match spot validation against BLASTP-short.

Requires NCBI BLAST+ (`makeblastdb`, `blastp`). It chooses deterministic peptides
from a UniProt FASTA, establishes locus/protein counts by direct raw-FASTA scan,
then asks BLASTP-short for independent 100%-identity full-length protein hits.
BLAST validates protein-level exact mapping; repeated loci within one protein are
validated by the browser/core fixture tests because BLAST HSP reporting can merge hits.
"""
from __future__ import annotations
import argparse, gzip, shutil, subprocess, tempfile
from pathlib import Path
from build_proteomes import parse_fasta


def protein_hits(records, pep):
    return sorted(p.accession for p in records if pep in p.sequence)

def choose(records):
    candidates=[]
    for p in records[:500]:
        if len(p.sequence) >= 24:
            mid=max(0,len(p.sequence)//2-6); candidates.append(p.sequence[mid:mid+12])
    unique=next((x for x in candidates if len(protein_hits(records,x))==1),None)
    seen={}; multi=None
    for p in records[:1500]:
        for i in range(0,max(0,len(p.sequence)-7+1),3):
            k=p.sequence[i:i+7]
            if k in seen and seen[k]!=p.accession: multi=k; break
            seen[k]=p.accession
        if multi: break
    if not unique or not multi: raise RuntimeError("Could not construct validation peptides")
    return {"unique":unique,"multiple":multi}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--fasta",type=Path,required=True); args=ap.parse_args()
    if not shutil.which("blastp") or not shutil.which("makeblastdb"):
        raise SystemExit("BLAST+ not installed")
    raw=args.fasta.read_bytes(); records=parse_fasta(raw,True); tests=choose(records)
    with tempfile.TemporaryDirectory() as td:
        td=Path(td); fa=td/'db.fasta'
        text=gzip.decompress(raw).decode() if raw[:2]==b'\x1f\x8b' else raw.decode(); fa.write_text(text)
        subprocess.run(["makeblastdb","-in",str(fa),"-dbtype","prot","-parse_seqids"],check=True,stdout=subprocess.DEVNULL)
        q=td/'q.fasta'; q.write_text(''.join(f'>{k}\n{v}\n' for k,v in tests.items()))
        out=td/'out.tsv'
        subprocess.run(["blastp","-task","blastp-short","-query",str(q),"-db",str(fa),"-seg","no","-evalue","1000000","-max_target_seqs","100000","-outfmt","6 qseqid sseqid pident length qlen","-out",str(out)],check=True)
        blast={k:set() for k in tests}
        for line in out.read_text().splitlines():
            qid,sid,pident,length,qlen=line.split('\t')
            if float(pident)==100 and int(length)==int(qlen): blast[qid].add(sid.split('|')[1] if '|' in sid else sid)
        bad=[]
        for k,pep in tests.items():
            expected=set(protein_hits(records,pep)); got=blast[k]
            print(k,pep,"raw_fasta",len(expected),"blast",len(got))
            if expected!=got: bad.append((k,sorted(expected-got),sorted(got-expected)))
        if bad: raise SystemExit(f"BLAST validation discrepancy: {bad}")
        print("BLASTP-short validation: no discrepancies for selected full-length exact matches.")
if __name__=='__main__': main()
