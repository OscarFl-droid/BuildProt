import test from 'node:test';
import assert from 'node:assert/strict';
import { parsePeptideInput, parseContactPositions } from '../js/parser.js';

test('parses bare, tab-separated, comma-separated and Excel-style sequences',()=>{
  const rows=parsePeptideInput('dvfqeliapk\nid2\t ACDEFGHIK \nid3,PEPTIDER');
  assert.deepEqual(rows.map(x=>x.peptide),['DVFQELIAPK','ACDEFGHIK','PEPTIDER']);
  assert.deepEqual(rows.map(x=>x.inputId),['PEP_0001','id2','id3']);
});

test('retains duplicate inputs and identifies them',()=>{
  const rows=parsePeptideInput('x1\tACDEFGHIK\nx2\tacdefghik');
  assert.equal(rows.length,2); assert.equal(rows[1].duplicateOf,'x1');
});

test('flags non-standard residues',()=>{
  const rows=parsePeptideInput('bad\tACDXFGHIK');
  assert.equal(rows[0].valid,false); assert.deepEqual(rows[0].invalidChars,['X']);
});

test('contact position parser accepts ranges and clips by peptide length',()=>{
  assert.deepEqual(parseContactPositions('2,5-7,99',9),[2,5,6,7]);
});
