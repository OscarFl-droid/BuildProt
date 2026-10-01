import test from 'node:test';
import assert from 'node:assert/strict';
import { SuffixProteome, mismatchInfo } from '../js/search-core.js';
import { uniquenessScore, blosumScore } from '../js/scoring.js';

function fixture(records) {
  let sequence=''; const proteins=[];
  for (const r of records) {
    if (sequence) sequence+='|';
    const start=sequence.length; sequence+=r.sequence;
    proteins.push({accession:r.accession,gene:r.gene,protein:r.protein||r.accession,family:r.family||'',isCanonical:r.isCanonical!==false,start,length:r.sequence.length});
  }
  const sa=Uint32Array.from(Array.from({length:sequence.length},(_,i)=>i).sort((a,b)=>sequence.slice(a)<sequence.slice(b)?-1:sequence.slice(a)>sequence.slice(b)?1:0));
  return new SuffixProteome(sequence,sa,proteins);
}

const human=fixture([
  {accession:'H1',gene:'GENE1',sequence:'MSTDVFQELIAPKAAAQQQ'},
  {accession:'H2',gene:'PARA1',family:'Kinase family',sequence:'GGGACDEFGHIKTTT'},
  {accession:'H3',gene:'PARA2',family:'Kinase family',sequence:'GGGACDEYGHIKTTT'},
  {accession:'H4',gene:'REP1',sequence:'QQQPEPTIDERAAA'},
  {accession:'H5',gene:'REP2',sequence:'TTTPEPTIDERCCC'},
  {accession:'H6',gene:'INTRA',sequence:'AAAAKQQQAAAAK'},
  {accession:'H2-2',gene:'PARA1',family:'Kinase family',sequence:'GGGACDEFGHIKSSS',isCanonical:false}
]);
const mouse=fixture([
  {accession:'M1',gene:'Gene1',sequence:'MSTDVFQELVAPKAAAQQQ'},
  {accession:'M2',gene:'Rep1',sequence:'QQQPEPTIDERAAA'}
]);

test('exactly unique peptide maps to one proteomic locus',()=>{
  const hits=human.exactMatches('DVFQELIAPK','canonical');
  assert.equal(hits.length,1); assert.equal(hits[0].accession,'H1'); assert.equal(hits[0].position,'4-13');
});

test('multiple proteins and repeated loci are counted separately',()=>{
  assert.equal(human.exactMatches('PEPTIDER','canonical').length,2);
  const intra=human.exactMatches('AAAAK','canonical');
  assert.equal(intra.length,2); assert.equal(new Set(intra.map(x=>x.accession)).size,1);
});

test('canonical universe excludes reviewed isoforms but extended includes them',()=>{
  assert.equal(human.exactMatches('ACDEFGHIK','canonical').length,1);
  assert.equal(human.exactMatches('ACDEFGHIK','extended').length,2);
});

test('nearest same-length off-target recovers one-residue paralog window',()=>{
  const exact=human.exactMatches('ACDEFGHIK','canonical');
  const h=human.nearestOfftargets('ACDEFGHIK',exact,{universe:'canonical',topK:5});
  assert.equal(h.offTargets[0].sequence,'ACDEYGHIK');
  assert.equal(h.offTargets[0].mismatches,1);
  assert.deepEqual(h.offTargets[0].mismatchPositions,[5]);
  assert.ok(Math.abs(h.offTargets[0].identity-88.8888889)<1e-5);
});

test('contact weighting finds global worst contact-preserving neighbour',()=>{
  const exact=human.exactMatches('ACDEFGHIK','canonical');
  const h=human.nearestOfftargets('ACDEFGHIK',exact,{universe:'canonical',contactPositions:[5],contactMultiplier:5});
  assert.ok(h.weightedBest);
  assert.ok(h.weightedBest.weightedIdentity>=0 && h.weightedBest.weightedIdentity<=100);
});

test('mouse versus human difference',()=>{
  assert.equal(human.exactMatches('DVFQELIAPK','canonical').length,1);
  assert.equal(mouse.exactMatches('DVFQELIAPK','canonical').length,0);
});

test('score gate and headroom are exact',()=>{
  assert.equal(uniquenessScore(2,100).score,0);
  assert.equal(uniquenessScore(0,null).score,null);
  assert.ok(Math.abs(uniquenessScore(1,80).score-20)<1e-12);
});

test('BLOSUM62 normalization gives 100 for self comparison',()=>{
  const s=blosumScore('ACDEFGHIK','ACDEFGHIK');
  assert.equal(s.normalized,100);
});

test('suffix-array nearest results agree with independent brute-force Hamming search',()=>{
  const q='DVFQELIAPK'; const exact=human.exactMatches(q,'canonical');
  const got=human.nearestOfftargets(q,exact,{universe:'canonical',topK:5}).offTargets.map(x=>[x.sequence,x.mismatches,x.accession,x.position]);
  const all=[];
  for(let pi=0;pi<human.proteins.length;pi++){
    const p=human.proteins[pi]; if(!p.isCanonical||p.length<q.length)continue;
    for(let rel=0;rel<=p.length-q.length;rel++){
      const g=p.start+rel;if(g===exact[0].globalStart)continue;
      const s=human.sequence.slice(g,g+q.length);const m=mismatchInfo(q,s).mismatches;
      all.push([s,m,p.accession,`${rel+1}-${rel+q.length}`]);
    }
  }
  all.sort((a,b)=>a[1]-b[1]||a[2].localeCompare(b[2])||Number(a[3].split('-')[0])-Number(b[3].split('-')[0]));
  assert.deepEqual(got,all.slice(0,5));
});
