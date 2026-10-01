import { blosumScore, blosumMismatchPenalties, contactWeights } from './scoring.js';

function compareSuffixToPattern(sequence, pos, pattern) {
  const n = pattern.length;
  for (let i = 0; i < n; i++) {
    const s = sequence.charCodeAt(pos + i) || -1;
    const p = pattern.charCodeAt(i);
    if (s < p) return -1;
    if (s > p) return 1;
  }
  return 0;
}

function lowerBound(sequence, sa, pattern) {
  let lo = 0, hi = sa.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (compareSuffixToPattern(sequence, sa[mid], pattern) < 0) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

function upperBound(sequence, sa, pattern) {
  let lo = 0, hi = sa.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (compareSuffixToPattern(sequence, sa[mid], pattern) <= 0) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

export function splitSegments(length, mismatches) {
  const parts = Math.min(length, mismatches + 1);
  const base = Math.floor(length / parts);
  const rem = length % parts;
  const out = [];
  let offset = 0;
  for (let i = 0; i < parts; i++) {
    const size = base + (i < rem ? 1 : 0);
    out.push({ offset, length: size });
    offset += size;
  }
  return out;
}

export function mismatchInfo(query, subject, weights = null) {
  const positions = [];
  let weightedMismatch = 0;
  let totalWeight = 0;
  for (let i = 0; i < query.length; i++) {
    const w = weights ? weights[i] : 1;
    totalWeight += w;
    if (query[i] !== subject[i]) {
      positions.push(i + 1);
      weightedMismatch += w;
    }
  }
  const mismatches = positions.length;
  return {
    mismatches,
    mismatchPositions: positions,
    identity: 100 * (query.length - mismatches) / query.length,
    weightedIdentity: 100 * (1 - weightedMismatch / totalWeight),
    weightedMismatch,
    totalWeight
  };
}

export class SuffixProteome {
  constructor(sequence, suffixArray, proteins) {
    this.sequence = sequence;
    this.sa = suffixArray;
    this.proteins = proteins;
    this.starts = proteins.map(p => p.start);
  }

  allowed(proteinIndex, universe) {
    return universe === 'extended' || this.proteins[proteinIndex].isCanonical;
  }

  proteinIndexAt(globalPos) {
    let lo = 0, hi = this.starts.length - 1, best = -1;
    while (lo <= hi) {
      const mid = (lo + hi) >>> 1;
      if (this.starts[mid] <= globalPos) { best = mid; lo = mid + 1; }
      else hi = mid - 1;
    }
    if (best < 0) return -1;
    const p = this.proteins[best];
    return globalPos < p.start + p.length ? best : -1;
  }

  occurrenceRange(pattern) {
    return [lowerBound(this.sequence, this.sa, pattern), upperBound(this.sequence, this.sa, pattern)];
  }

  exactMatches(query, universe = 'canonical') {
    const [lo, hi] = this.occurrenceRange(query);
    const hits = [];
    for (let i = lo; i < hi; i++) {
      const start = this.sa[i];
      const pi = this.proteinIndexAt(start);
      if (pi < 0 || !this.allowed(pi, universe)) continue;
      const p = this.proteins[pi];
      if (start + query.length > p.start + p.length) continue;
      hits.push(this.formatHit(start, pi, query.length));
    }
    hits.sort((a, b) => a.accession.localeCompare(b.accession) || a.start - b.start);
    return hits;
  }

  formatHit(globalStart, proteinIndex, length) {
    const p = this.proteins[proteinIndex];
    const rel = globalStart - p.start;
    return {
      globalStart,
      proteinIndex,
      accession: p.accession,
      gene: p.gene || '',
      protein: p.protein || '',
      family: p.family || '',
      isCanonical: !!p.isCanonical,
      positionStart: rel + 1,
      positionEnd: rel + length,
      position: `${rel + 1}-${rel + length}`
    };
  }

  _verifyCandidate(globalStart, query, universe, threshold, weights, advanced) {
    const pi = this.proteinIndexAt(globalStart);
    if (pi < 0 || !this.allowed(pi, universe)) return null;
    const p = this.proteins[pi];
    if (globalStart + query.length > p.start + p.length) return null;
    const subject = this.sequence.slice(globalStart, globalStart + query.length);
    let mismatches = 0;
    const mismatchPositions = [];
    let weightedMismatch = 0, totalWeight = 0;
    for (let i = 0; i < query.length; i++) {
      const w = weights[i]; totalWeight += w;
      if (query[i] !== subject[i]) {
        mismatches += 1; mismatchPositions.push(i + 1); weightedMismatch += w;
        if (threshold != null && mismatches > threshold) return null;
      }
    }
    const identity = 100 * (query.length - mismatches) / query.length;
    const out = {
      ...this.formatHit(globalStart, pi, query.length),
      sequence: subject,
      mismatches,
      mismatchPositions,
      identity,
      weightedIdentity: 100 * (1 - weightedMismatch / totalWeight)
    };
    if (advanced) Object.assign(out, {
      blosumRaw: blosumScore(query, subject).raw,
      blosumNormalized: blosumScore(query, subject).normalized
    });
    return out;
  }

  _allWindows(query, universe, excludedStart, weights, advanced, candidates) {
    for (let pi = 0; pi < this.proteins.length; pi++) {
      if (!this.allowed(pi, universe)) continue;
      const p = this.proteins[pi];
      if (p.length < query.length) continue;
      for (let rel = 0; rel <= p.length - query.length; rel++) {
        const g = p.start + rel;
        if (g === excludedStart || candidates.has(g)) continue;
        const v = this._verifyCandidate(g, query, universe, null, weights, advanced);
        if (v) candidates.set(g, v);
      }
    }
  }

  nearestOfftargets(query, exactMatches, {
    universe = 'canonical',
    topK = 5,
    advanced = false,
    contactPositions = [],
    contactMultiplier = 2
  } = {}) {
    if (!exactMatches.length) return { offTargets: [], weightedBest: null, blosumBest: null };
    const excludedStart = exactMatches[0].globalStart;
    const weights = contactWeights(query.length, contactPositions, contactMultiplier);
    const candidates = new Map();
    const sortedWeights = [...weights].sort((a, b) => a - b);
    const blosumPenalties = blosumMismatchPenalties(query).sort((a, b) => a - b);
    const selfRaw = advanced ? blosumScore(query, query).raw : null;
    let primaryDone = false, weightedDone = false, advancedDone = !advanced;
    let weightedBest = null, blosumBest = null;

    const refreshBest = () => {
      for (const c of candidates.values()) {
        if (!weightedBest || c.weightedIdentity > weightedBest.weightedIdentity ||
          (c.weightedIdentity === weightedBest.weightedIdentity && c.identity > weightedBest.identity)) weightedBest = c;
        if (advanced && (!blosumBest || c.blosumRaw > blosumBest.blosumRaw ||
          (c.blosumRaw === blosumBest.blosumRaw && c.identity > blosumBest.identity))) blosumBest = c;
      }
    };

    for (let d = 0; d <= query.length - 1; d++) {
      for (const seg of splitSegments(query.length, d)) {
        const pattern = query.slice(seg.offset, seg.offset + seg.length);
        const [lo, hi] = this.occurrenceRange(pattern);
        for (let i = lo; i < hi; i++) {
          const g = this.sa[i] - seg.offset;
          if (g < 0 || g === excludedStart || candidates.has(g)) continue;
          const v = this._verifyCandidate(g, query, universe, d, weights, advanced);
          if (v) candidates.set(g, v);
        }
      }
      refreshBest();
      const withinD = [...candidates.values()].filter(c => c.mismatches <= d);
      primaryDone = withinD.length >= topK;

      const unseenMismatchCount = d + 1;
      const minWeightedLoss = sortedWeights.slice(0, unseenMismatchCount).reduce((a, b) => a + b, 0);
      const totalWeight = weights.reduce((a, b) => a + b, 0);
      const unseenWeightedUpper = 100 * (1 - minWeightedLoss / totalWeight);
      weightedDone = !!weightedBest && unseenWeightedUpper <= weightedBest.weightedIdentity;

      if (advanced && blosumBest) {
        const minBlosumLoss = blosumPenalties.slice(0, unseenMismatchCount).reduce((a, b) => a + b, 0);
        const unseenRawUpper = selfRaw - minBlosumLoss;
        advancedDone = unseenRawUpper <= blosumBest.blosumRaw;
      }
      if (primaryDone && weightedDone && advancedDone) break;
    }

    if (!primaryDone || !weightedDone || !advancedDone) {
      this._allWindows(query, universe, excludedStart, weights, advanced, candidates);
      refreshBest();
    }

    const offTargets = [...candidates.values()]
      .sort((a, b) => a.mismatches - b.mismatches || b.identity - a.identity || a.accession.localeCompare(b.accession) || a.positionStart - b.positionStart)
      .slice(0, topK);
    return { offTargets, weightedBest, blosumBest };
  }
}

export function classifyExactRelationship(matches) {
  if (matches.length <= 1) return 'single proteomic locus';
  const genes = new Set(matches.map(m => m.gene).filter(Boolean));
  if (genes.size <= 1) return 'same gene / isoforms / repeated locus';
  const familySets = matches.map(m => new Set(String(m.family || '').split(/[;,]/).map(x => x.trim()).filter(Boolean)));
  let shared = familySets.length ? new Set(familySets[0]) : new Set();
  for (const s of familySets.slice(1)) shared = new Set([...shared].filter(x => s.has(x)));
  if (shared.size) return 'probable paralogous proteins (shared UniProt family annotation)';
  return 'multiple genes — apparently unrelated or family relationship unannotated';
}
