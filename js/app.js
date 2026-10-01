import { parsePeptideInput, parseContactPositions } from './parser.js';

const $ = s => document.querySelector(s);
const input = $('#peptideInput');
const parseSummary = $('#parseSummary');
const speciesSelect = $('#speciesSelect');
const universeSelect = $('#universeSelect');
const advancedToggle = $('#advancedToggle');
const contactPositionsInput = $('#contactPositions');
const contactMultiplierInput = $('#contactMultiplier');
const datasetInfo = $('#datasetInfo');
const analyzeBtn = $('#analyzeBtn');
const cancelBtn = $('#cancelBtn');
const progressWrap = $('#progressWrap');
const progressBar = $('#progressBar');
const progressText = $('#progressText');
const table = $('#resultsTable');
const filterInput = $('#filterInput');
const resultSummary = $('#resultSummary');
const exportBtns = [$('#copyBtn'), $('#csvBtn'), $('#tsvBtn'), $('#excelBtn')];

const COLUMNS = [
  ['Input_ID','inputId'],['Peptide','peptide'],['Length','length'],['Species','species'],['Proteome_version','proteomeVersion'],
  ['Found_in_proteome','foundInProteome'],['Exact_match_count','exactMatchCount'],['Protein_match_count','proteinMatchCount'],
  ['Proteome_unique','proteomeUnique'],['Target_gene','targetGene'],['Target_accession','targetAccession'],['Target_protein','targetProtein'],['Target_position','targetPosition'],
  ['Nearest_offtarget_sequence','nearestOfftargetSequence'],['Nearest_offtarget_identity_%','nearestOfftargetIdentity'],['Nearest_offtarget_mismatches','nearestOfftargetMismatches'],
  ['Nearest_offtarget_gene','nearestOfftargetGene'],['Nearest_offtarget_accession','nearestOfftargetAccession'],['Proteome_Unique_Score','proteomeUniqueScore'],
  ['ContactWeightedUniqueness','contactWeightedUniqueness'],['Advanced_similarity_score','advancedSimilarityScore'],['Status','status'],['Notes','notes']
];

let manifest = null;
let datasets = [];
let worker = null;
let loadedDatasetId = null;
let workerReady = false;
let pendingReadyResolve = null;
let currentRows = [];
let displayRows = [];
let sortState = { key: null, dir: 1 };
let parsedRows = [];
let activeRun = null;

function formatNumber(v) {
  return typeof v === 'number' && Number.isFinite(v) ? Number(v.toFixed(2)) : (v ?? '');
}
function boolText(v) { return v === true ? 'TRUE' : v === false ? 'FALSE' : ''; }
function cellValue(row, key) {
  if (['foundInProteome','proteomeUnique'].includes(key)) return boolText(row[key]);
  if (['nearestOfftargetIdentity','proteomeUniqueScore','contactWeightedUniqueness','advancedSimilarityScore'].includes(key)) return formatNumber(row[key]);
  return row[key] ?? '';
}

function updateParseSummary() {
  parsedRows = parsePeptideInput(input.value);
  const valid = parsedRows.filter(r => r.valid).length;
  const invalid = parsedRows.length - valid;
  const dup = parsedRows.filter(r => r.duplicateOf).length;
  parseSummary.textContent = parsedRows.length ? `${parsedRows.length} input rows · ${valid} valid · ${dup} duplicate input${dup===1?'':'s'} · ${invalid} invalid` : 'No peptides parsed.';
}
input.addEventListener('input', updateParseSummary);

function selectedDataset() { return datasets.find(d => d.id === speciesSelect.value); }
function updateDatasetInfo() {
  const d = selectedDataset();
  if (!d) { datasetInfo.textContent = manifest?.message || 'No proteome bundle available.'; return; }
  datasetInfo.innerHTML = `<strong>${d.label}</strong> · ${d.proteome_id} · UniProt ${d.uniprot_release || 'release not reported'} · built ${d.build_date || 'unknown date'} · ${Number(d.canonical_proteins||0).toLocaleString()} canonical proteins${d.reviewed_isoforms ? ` + ${Number(d.reviewed_isoforms).toLocaleString()} reviewed isoforms` : ''}.`;
}
speciesSelect.addEventListener('change', () => { updateDatasetInfo(); loadedDatasetId = null; workerReady = false; });

async function initManifest() {
  const r = await fetch('./data/manifest.json', { cache: 'no-store' });
  manifest = await r.json();
  datasets = manifest.datasets || [];
  speciesSelect.innerHTML = '';
  for (const d of datasets) {
    const o = document.createElement('option'); o.value = d.id; o.textContent = d.label; speciesSelect.appendChild(o);
  }
  if (!datasets.length) {
    const o = document.createElement('option'); o.textContent = 'Proteome data not built'; o.value = ''; speciesSelect.appendChild(o);
    analyzeBtn.disabled = true;
  }
  updateDatasetInfo();
}

function makeWorker() {
  if (worker) worker.terminate();
  workerReady = false;
  loadedDatasetId = null;
  worker = new Worker('./js/worker.js', { type: 'module' });
  worker.onmessage = e => {
    const m = e.data;
    if (m.type === 'load-progress') setProgress(Math.round(m.fraction * 100), m.message);
    else if (m.type === 'ready') {
      workerReady = true;
      setProgress(100, `Proteome loaded: ${Number(m.proteinCount).toLocaleString()} proteins/isoforms.`);
      if (pendingReadyResolve) { pendingReadyResolve(); pendingReadyResolve = null; }
    } else if (m.type === 'analysis-progress') {
      const p = Math.round(100 * m.completed / m.total);
      setProgress(p, `Analysing unique peptide ${m.completed.toLocaleString()} of ${m.total.toLocaleString()}…`);
    } else if (m.type === 'results') finishAnalysis(m.results);
    else if (m.type === 'error') failAnalysis(m.message);
  };
}
makeWorker();

function setProgress(percent, text) {
  progressWrap.classList.remove('hidden'); progressBar.style.width = `${Math.max(0,Math.min(100,percent))}%`; progressText.textContent = text;
}
async function ensureDatasetLoaded() {
  const d = selectedDataset();
  if (!d) throw new Error('No built proteome dataset is available.');
  if (loadedDatasetId === d.id && workerReady) return;
  loadedDatasetId = d.id; workerReady = false;
  setProgress(2, 'Starting proteome loader…');
  const ready = new Promise(resolve => pendingReadyResolve = resolve);
  const workerDataset = structuredClone(d);
  workerDataset.files = Object.fromEntries(Object.entries(d.files).map(([k, v]) => [k, new URL(v, document.baseURI).href]));
  worker.postMessage({ type: 'init', dataset: workerDataset });
  await ready;
}

analyzeBtn.addEventListener('click', async () => {
  updateParseSummary();
  if (!parsedRows.length) { parseSummary.textContent = 'Paste at least one peptide first.'; return; }
  const valid = parsedRows.filter(r => r.valid);
  const uniquePeptides = [...new Set(valid.map(r => r.peptide))];
  analyzeBtn.disabled = true; cancelBtn.disabled = false;
  try {
    await ensureDatasetLoaded();
    const contactPositions = parseContactPositions(contactPositionsInput.value, 10000);
    const contactMultiplier = Math.max(1, Number(contactMultiplierInput.value) || 2);
    activeRun = { dataset: selectedDataset(), universe: universeSelect.value, advanced: advancedToggle.checked, contactPositions, contactMultiplier };
    setProgress(0, `Preparing ${uniquePeptides.length.toLocaleString()} unique peptide sequence${uniquePeptides.length===1?'':'s'}…`);
    worker.postMessage({ type:'analyze', peptides: uniquePeptides, options: { universe: activeRun.universe, advanced: activeRun.advanced, contactPositions, contactMultiplier, topK: 5 }});
  } catch (e) { failAnalysis(e.message); }
});

cancelBtn.addEventListener('click', () => {
  makeWorker(); analyzeBtn.disabled = datasets.length === 0; cancelBtn.disabled = true; setProgress(0, 'Worker reset.');
});

function invalidResult(r, d) {
  return {
    inputId:r.inputId, peptide:r.peptide, length:r.length, species:d?.label||'', proteomeVersion:versionText(d, activeRun?.universe || universeSelect.value),
    foundInProteome:null, exactMatchCount:null, proteinMatchCount:null, proteomeUnique:null,
    targetGene:'',targetAccession:'',targetProtein:'',targetPosition:'',nearestOfftargetSequence:'',nearestOfftargetIdentity:null,nearestOfftargetMismatches:null,
    nearestOfftargetGene:'',nearestOfftargetAccession:'',proteomeUniqueScore:null,contactWeightedUniqueness:null,advancedSimilarityScore:null,
    status:'INVALID INPUT',notes:`Rejected non-standard amino-acid character(s): ${r.invalidChars.join(', ') || 'empty sequence'}.`,exactMatches:[],offTargets:[]
  };
}
function versionText(d, universe = universeSelect.value) { return d ? `${d.proteome_id}; UniProt ${d.uniprot_release || 'unknown'}; ${universe === 'canonical' ? 'canonical/reference' : 'canonical + reviewed isoforms'}; build ${d.build_date || 'unknown'}` : ''; }

function finishAnalysis(uniqueResults) {
  const byPep = new Map(uniqueResults.map(x => [x.peptide, x]));
  const d = activeRun?.dataset || selectedDataset();
  currentRows = parsedRows.map(r => {
    if (!r.valid) return invalidResult(r,d);
    const base = structuredClone(byPep.get(r.peptide));
    const duplicateNote = r.duplicateOf ? ` Duplicate input of ${r.duplicateOf}; analysis reused from cache.` : '';
    return { ...base, inputId:r.inputId, species:d.label, proteomeVersion:versionText(d, activeRun?.universe || universeSelect.value), notes:(base.notes||'')+duplicateNote, duplicateOf:r.duplicateOf };
  });
  displayRows = [...currentRows];
  renderTable();
  const nUnique = currentRows.filter(r => r.status==='UNIQUE').length;
  const nMulti = currentRows.filter(r => r.status==='MULTIPLE EXACT MATCHES').length;
  const nAbsent = currentRows.filter(r => r.status==='NOT FOUND').length;
  const nInvalid = currentRows.filter(r => r.status==='INVALID INPUT').length;
  resultSummary.textContent = `${currentRows.length.toLocaleString()} rows · ${nUnique} unique · ${nMulti} multiple exact matches · ${nAbsent} not found · ${nInvalid} invalid.`;
  exportBtns.forEach(b => b.disabled = false);
  analyzeBtn.disabled = false; cancelBtn.disabled = true;
  setProgress(100, 'Analysis complete.');
}
function failAnalysis(message) {
  analyzeBtn.disabled = datasets.length === 0; cancelBtn.disabled = true; setProgress(0, `Error: ${message}`); progressBar.style.width='0%';
}

function renderTable() {
  const thead = table.querySelector('thead'); const tbody = table.querySelector('tbody');
  thead.innerHTML = `<tr>${COLUMNS.map(([label,key]) => `<th data-key="${key}">${label}${sortState.key===key?(sortState.dir>0?' ▲':' ▼'):''}</th>`).join('')}</tr>`;
  tbody.innerHTML = '';
  const frag = document.createDocumentFragment();
  for (let i=0;i<displayRows.length;i++) {
    const row = displayRows[i]; const tr=document.createElement('tr'); tr.className='data-row'; tr.dataset.idx=String(i);
    tr.innerHTML=COLUMNS.map(([label,key])=>{const v=cellValue(row,key);const cls=key==='status'?(row.status==='UNIQUE'?'status-unique':row.status==='MULTIPLE EXACT MATCHES'||row.status==='INVALID INPUT'?'status-fail':''):'';return `<td class="${cls}" title="${escapeHtml(String(v))}">${escapeHtml(String(v))}</td>`}).join('');
    frag.appendChild(tr);
  }
  tbody.appendChild(frag);
  thead.querySelectorAll('th').forEach(th=>th.addEventListener('click',()=>sortBy(th.dataset.key)));
  tbody.querySelectorAll('tr.data-row').forEach(tr=>tr.addEventListener('click',()=>toggleDetails(tr)));
}
function escapeHtml(s){return s.replace(/[&<>"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));}

function sortBy(key){ if(sortState.key===key)sortState.dir*=-1;else sortState={key,dir:1}; const dir=sortState.dir;displayRows.sort((a,b)=>{const av=a[key],bv=b[key];if(av==null&&bv==null)return 0;if(av==null)return 1;if(bv==null)return -1;if(typeof av==='number'&&typeof bv==='number')return dir*(av-bv);return dir*String(av).localeCompare(String(bv),undefined,{numeric:true,sensitivity:'base'});});renderTable();}

function alignment(query, off){
  const marks=[...query].map((c,i)=>c===off[i]?' ':'^').join('');
  return `Query:     ${query}\nOfftarget: ${off}\n           ${marks}`;
}
function hitLine(h){return `${escapeHtml(h.accession)} · ${escapeHtml(h.gene||'—')} · ${escapeHtml(h.protein||'—')} · aa ${escapeHtml(h.position)}`;}
function toggleDetails(tr){
  const next=tr.nextElementSibling;if(next?.classList.contains('detail-row')){next.remove();return;}
  const row=displayRows[Number(tr.dataset.idx)];
  const d=document.createElement('tr');d.className='detail-row';const td=document.createElement('td');td.colSpan=COLUMNS.length;
  const exact=row.exactMatches?.length?row.exactMatches.map(h=>`<div class="hit">${hitLine(h)}${h.isCanonical?' <span class="pill">canonical</span>':' <span class="pill">isoform</span>'}</div>`).join(''):'<p>No exact match.</p>';
  const offs=row.offTargets?.length?row.offTargets.map((h,i)=>`<div class="hit"><strong>#${i+1}</strong> ${hitLine(h)} · identity ${formatNumber(h.identity)}% · ${h.mismatches} mismatch${h.mismatches===1?'':'es'} at ${h.mismatchPositions.join(',')||'none'}${h.blosumRaw!=null?` · BLOSUM62 raw ${h.blosumRaw}, normalized ${formatNumber(h.blosumNormalized)}`:''}<div class="alignment">${escapeHtml(alignment(row.peptide,h.sequence))}</div></div>`).join(''):'<p>No alternative same-length window available.</p>';
  const weighted=row.weightedBest?`<p><strong>Contact-weighted worst neighbour:</strong> ${escapeHtml(row.weightedBest.sequence)} (${formatNumber(row.weightedBest.weightedIdentity)}% weighted identity; ${escapeHtml(row.weightedBest.accession)}:${row.weightedBest.position}).</p>`:'';
  const blosum=row.blosumBest?`<p><strong>BLOSUM62-best alternative:</strong> ${escapeHtml(row.blosumBest.sequence)} · raw ${row.blosumBest.blosumRaw}, normalized ${formatNumber(row.blosumBest.blosumNormalized)} · ${escapeHtml(row.blosumBest.accession)}:${row.blosumBest.position}.</p>`:'';
  td.innerHTML=`<div class="detail-grid"><div class="detail-block"><h3>Exact proteomic loci</h3><p><strong>Relationship:</strong> ${escapeHtml(row.exactRelationship||'—')}</p>${exact}</div><div class="detail-block"><h3>Closest same-length off-targets</h3>${weighted}${blosum}${offs}</div></div>`;
  d.appendChild(td);tr.after(d);
}

filterInput.addEventListener('input',()=>{const q=filterInput.value.trim().toLowerCase();displayRows=!q?[...currentRows]:currentRows.filter(r=>COLUMNS.some(([,k])=>String(cellValue(r,k)).toLowerCase().includes(q)));renderTable();});

function exportText(delim=',', excel=false){
  const quote=v=>{let s=String(v??'');if(delim===','&&(s.includes(',')||s.includes('"')||s.includes('\n')))s='"'+s.replaceAll('"','""')+'"';return s;};
  const lines=[COLUMNS.map(([h])=>quote(h)).join(delim),...currentRows.map(r=>COLUMNS.map(([,k])=>quote(cellValue(r,k))).join(delim))];
  return (excel?'\ufeff':'')+lines.join('\r\n');
}
function download(name,text,type){const blob=new Blob([text],{type});const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);}
$('#csvBtn').addEventListener('click',()=>download('peptide_proteome_uniqueness.csv',exportText(','), 'text/csv;charset=utf-8'));
$('#excelBtn').addEventListener('click',()=>download('peptide_proteome_uniqueness_excel.csv',exportText(',',true),'text/csv;charset=utf-8'));
$('#tsvBtn').addEventListener('click',()=>download('peptide_proteome_uniqueness.tsv',exportText('\t'),'text/tab-separated-values;charset=utf-8'));
$('#copyBtn').addEventListener('click',async()=>{await navigator.clipboard.writeText(exportText('\t'));$('#copyBtn').textContent='Copied';setTimeout(()=>$('#copyBtn').textContent='Copy table',1200);});

$('#exampleBtn').addEventListener('click',()=>{input.value='pep_unique\tDVFQELIAPK\npep_repeat,DVFQELIAPK\nBAD\tACDXFGHIK\nGILFVGSGVSGK';updateParseSummary();});

await initManifest();
updateParseSummary();
