const AA = 'ARNDCQEGHILKMFPSTWYV';
const ROWS = [
  [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
  [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
  [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
  [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
  [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
  [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
  [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
  [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
  [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
  [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
  [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
  [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
  [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
  [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
  [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
  [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
  [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
  [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
  [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
  [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4]
];
const IDX = Object.fromEntries([...AA].map((a, i) => [a, i]));

export function blosum62(a, b) {
  return ROWS[IDX[a]][IDX[b]];
}

export function blosumScore(query, subject) {
  let raw = 0, self = 0, minPossible = 0;
  for (let i = 0; i < query.length; i++) {
    const qi = query[i];
    raw += blosum62(qi, subject[i]);
    self += blosum62(qi, qi);
    minPossible += Math.min(...ROWS[IDX[qi]]);
  }
  const denom = self - minPossible;
  const normalized = denom === 0 ? 100 : Math.max(0, Math.min(100, 100 * (raw - minPossible) / denom));
  return { raw, normalized, self, minPossible };
}

export function blosumMismatchPenalties(query) {
  return [...query].map(q => {
    const i = IDX[q];
    const self = ROWS[i][i];
    let bestOther = -Infinity;
    for (let j = 0; j < AA.length; j++) if (j !== i) bestOther = Math.max(bestOther, ROWS[i][j]);
    return self - bestOther;
  });
}

export function separationBand(score) {
  if (score == null || Number.isNaN(score)) return 'NA';
  if (score >= 80) return 'very high sequence separation';
  if (score >= 60) return 'high sequence separation';
  if (score >= 40) return 'moderate sequence separation';
  if (score >= 20) return 'limited separation';
  return 'very close proteomic neighbour';
}

export function uniquenessScore(exactCount, nearestIdentity) {
  if (exactCount === 0) return { score: null, status: 'NA — peptide not found' };
  if (exactCount > 1) return { score: 0, status: 'FAIL — exact sequence not proteome-unique' };
  const identity = nearestIdentity ?? 0;
  const score = 100 * (1 - identity / 100);
  return { score, status: separationBand(score) };
}

export function contactWeights(length, contactPositions = [], multiplier = 2) {
  const contacts = new Set(contactPositions);
  return Array.from({ length }, (_, i) => contacts.has(i + 1) ? multiplier : 1);
}
