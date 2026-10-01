import { SuffixProteome, classifyExactRelationship } from './search-core.js';
import { uniquenessScore, separationBand } from './scoring.js';

let index = null;
let dataset = null;
const analysisCache = new Map();

async function ungzipResponse(url, kind = 'arrayBuffer') {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`Failed to load ${url}: HTTP ${r.status}`);
  if (!url.endsWith('.gz')) return kind === 'text' ? r.text() : r.arrayBuffer();
  if (typeof DecompressionStream === 'undefined') throw new Error('This browser lacks DecompressionStream support required for compressed proteome bundles.');
  const stream = r.body.pipeThrough(new DecompressionStream('gzip'));
  const out = await new Response(stream)[kind]();
  return out;
}

async function loadDataset(config) {
  dataset = config;
  postMessage({ type: 'load-progress', message: 'Loading protein metadata…', fraction: 0.1 });
  const metadataText = await ungzipResponse(config.files.metadata, 'text');
  const metadata = JSON.parse(metadataText);
  postMessage({ type: 'load-progress', message: 'Loading concatenated proteome sequence…', fraction: 0.45 });
  const sequence = await ungzipResponse(config.files.sequence, 'text');
  postMessage({ type: 'load-progress', message: 'Loading suffix-array index…', fraction: 0.75 });
  const saBuffer = await ungzipResponse(config.files.suffix_array, 'arrayBuffer');
  const sa = new Uint32Array(saBuffer);
  index = new SuffixProteome(sequence, sa, metadata.proteins);
  analysisCache.clear();
  postMessage({ type: 'ready', proteinCount: metadata.proteins.length, residues: sequence.length, metadata: metadata.summary });
}

function uniqueJoin(items) {
  return [...new Set(items.filter(Boolean))].join('; ');
}

function targetFields(matches) {
  if (!matches.length) return { gene: '', accession: '', protein: '', position: '' };
  if (matches.length === 1) {
    const m = matches[0];
    return { gene: m.gene, accession: m.accession, protein: m.protein, position: m.position };
  }
  return {
    gene: uniqueJoin(matches.map(m => m.gene)),
    accession: uniqueJoin(matches.map(m => m.accession)),
    protein: uniqueJoin(matches.map(m => m.protein)),
    position: matches.map(m => `${m.accession}:${m.position}`).join('; ')
  };
}

function analyzePeptide(peptide, options) {
  const cacheKey = JSON.stringify([peptide, options.universe, options.advanced, options.contactPositions, options.contactMultiplier]);
  if (analysisCache.has(cacheKey)) return structuredClone(analysisCache.get(cacheKey));

  const exact = index.exactMatches(peptide, options.universe);
  const exactCount = exact.length;
  const proteinCount = new Set(exact.map(x => x.accession)).size;
  const geneCount = new Set(exact.map(x => x.gene).filter(Boolean)).size;
  const target = targetFields(exact);
  const relationship = classifyExactRelationship(exact);
  let homology = { offTargets: [], weightedBest: null, blosumBest: null };
  if (exactCount > 0) homology = index.nearestOfftargets(peptide, exact, options);
  const nearest = homology.offTargets[0] ?? null;
  const baseScore = uniquenessScore(exactCount, nearest?.identity ?? 0);
  let contactScore = null;
  if (exactCount > 1) contactScore = 0;
  else if (exactCount === 1) contactScore = homology.weightedBest ? 100 - homology.weightedBest.weightedIdentity : 100;

  let status = 'NOT FOUND';
  if (exactCount === 1) status = 'UNIQUE';
  else if (exactCount > 1) status = 'MULTIPLE EXACT MATCHES';

  const notes = [];
  if (exactCount === 1) notes.push(`Engineering prioritisation band: ${separationBand(baseScore.score)}.`);
  if (exactCount > 1) notes.push(`Exact-match relationship: ${relationship}.`);
  if (options.contactPositions?.length) notes.push(`Contact-weighted score uses positions ${options.contactPositions.join(',')} at ×${options.contactMultiplier} weight.`);
  if (exactCount === 0) notes.push('No score assigned: verify sequence, species, and selected proteome universe.');

  const result = {
    peptide,
    length: peptide.length,
    foundInProteome: exactCount > 0,
    exactMatchCount: exactCount,
    proteinMatchCount: proteinCount,
    geneMatchCount: geneCount,
    proteomeUnique: exactCount === 1,
    exactRelationship: relationship,
    targetGene: target.gene,
    targetAccession: target.accession,
    targetProtein: target.protein,
    targetPosition: target.position,
    nearestOfftargetSequence: nearest?.sequence ?? '',
    nearestOfftargetIdentity: nearest?.identity ?? null,
    nearestOfftargetMismatches: nearest?.mismatches ?? null,
    nearestOfftargetGene: nearest?.gene ?? '',
    nearestOfftargetAccession: nearest?.accession ?? '',
    proteomeUniqueScore: baseScore.score,
    contactWeightedUniqueness: contactScore,
    advancedSimilarityScore: options.advanced ? (homology.blosumBest?.blosumNormalized ?? null) : null,
    advancedSimilarityRaw: options.advanced ? (homology.blosumBest?.blosumRaw ?? null) : null,
    status,
    notes: notes.join(' '),
    exactMatches: exact,
    offTargets: homology.offTargets,
    weightedBest: homology.weightedBest,
    blosumBest: homology.blosumBest
  };
  analysisCache.set(cacheKey, structuredClone(result));
  return result;
}

self.onmessage = async (event) => {
  const msg = event.data;
  try {
    if (msg.type === 'init') {
      await loadDataset(msg.dataset);
      return;
    }
    if (msg.type === 'analyze') {
      if (!index) throw new Error('Proteome index not loaded.');
      const { peptides, options } = msg;
      const results = [];
      for (let i = 0; i < peptides.length; i++) {
        results.push(analyzePeptide(peptides[i], options));
        if (i === peptides.length - 1 || i % 5 === 0) {
          postMessage({ type: 'analysis-progress', completed: i + 1, total: peptides.length });
        }
      }
      postMessage({ type: 'results', results });
    }
  } catch (error) {
    postMessage({ type: 'error', message: error?.message || String(error), stack: error?.stack || '' });
  }
};
