import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';

// Exercise the actual browser Worker module with compressed static fixture assets.
globalThis.self = globalThis;
const waiters = [];
globalThis.postMessage = msg => {
  for (let i = waiters.length - 1; i >= 0; i--) {
    if (waiters[i].type === msg.type) {
      const w = waiters.splice(i, 1)[0]; w.resolve(msg);
    }
  }
};
function waitFor(type) { return new Promise(resolve => waiters.push({type, resolve})); }

globalThis.fetch = async url => {
  const data = await fs.readFile(String(url));
  return new Response(data, {status: 200});
};

await import('../js/worker.js');

test('actual worker loads gzip bundle and returns expected homology result', async () => {
  const root = new URL('./browser-data/', import.meta.url).pathname;
  const ready = waitFor('ready');
  await self.onmessage({data:{type:'init',dataset:{files:{
    sequence: root+'fixture.sequence.txt.gz',
    suffix_array: root+'fixture.sa.bin.gz',
    metadata: root+'fixture.metadata.json.gz'
  }}}});
  await ready;
  const done = waitFor('results');
  await self.onmessage({data:{type:'analyze',peptides:['ACDEFGHIK'],options:{universe:'canonical',advanced:true,contactPositions:[5],contactMultiplier:2,topK:5}}});
  const msg = await done;
  const r = msg.results[0];
  assert.equal(r.exactMatchCount,1);
  assert.equal(r.proteomeUnique,true);
  assert.equal(r.nearestOfftargetSequence,'ACDEYGHIK');
  assert.equal(r.nearestOfftargetMismatches,1);
  assert.equal(r.nearestOfftargetGene,'PARA2');
  assert.equal(typeof r.advancedSimilarityScore,'number');
});
