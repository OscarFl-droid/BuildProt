export const STANDARD_AA = /^[ACDEFGHIKLMNPQRSTVWY]+$/;

function cleanToken(token) {
  return token.replace(/^['"]|['"]$/g, '').replace(/\s+/g, '').toUpperCase();
}

function looksSequenceLike(token) {
  const s = cleanToken(token);
  return s.length > 0 && /^[A-Z*.-]+$/.test(s);
}

export function parseContactPositions(text, length) {
  if (!text || !text.trim()) return [];
  const out = new Set();
  for (const raw of text.split(',')) {
    const part = raw.trim();
    if (!part) continue;
    const m = part.match(/^(\d+)\s*-\s*(\d+)$/);
    if (m) {
      let a = Number(m[1]), b = Number(m[2]);
      if (a > b) [a, b] = [b, a];
      for (let i = a; i <= b; i++) if (i >= 1 && i <= length) out.add(i);
      continue;
    }
    const n = Number(part);
    if (Number.isInteger(n) && n >= 1 && n <= length) out.add(n);
  }
  return [...out].sort((a, b) => a - b);
}

export function parsePeptideInput(text) {
  const rows = [];
  const firstSeen = new Map();
  const lines = String(text ?? '').split(/\r?\n/);
  let ordinal = 0;

  for (const originalLine of lines) {
    if (!originalLine.trim()) continue;
    ordinal += 1;
    let fields;
    if (originalLine.includes('\t')) fields = originalLine.split('\t');
    else if (originalLine.includes(',')) fields = originalLine.split(',');
    else fields = [originalLine];
    fields = fields.map(x => x.trim()).filter((x, i) => x.length || i === 0);

    let id = `PEP_${String(ordinal).padStart(4, '0')}`;
    let seqField = fields[0] ?? '';
    if (fields.length > 1) {
      id = fields[0] || id;
      const candidates = fields.slice(1).filter(looksSequenceLike);
      seqField = candidates.length ? candidates[candidates.length - 1] : fields[fields.length - 1];
    }

    const peptide = cleanToken(seqField);
    const invalidChars = [...new Set([...peptide].filter(c => !'ACDEFGHIKLMNPQRSTVWY'.includes(c)))];
    const valid = peptide.length > 0 && invalidChars.length === 0;
    const duplicateOf = valid && firstSeen.has(peptide) ? firstSeen.get(peptide) : null;
    if (valid && !firstSeen.has(peptide)) firstSeen.set(peptide, id);

    rows.push({
      inputId: id,
      peptide,
      originalLine,
      valid,
      invalidChars,
      duplicateOf,
      length: peptide.length
    });
  }
  return rows;
}
