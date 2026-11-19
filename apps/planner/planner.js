// The DFT budget planning rule in JavaScript, for the planner page (docs/planner/). It is a port of tsso/plan.py and
// must give the same plan: tests/test_planner_page.py runs both on the same inputs and compares them.
// No DOM in this file, so it also runs under node.
(function (root) {
  'use strict';
  const ACTIONS = ['run_dft', 'trust', 'hold', 'regenerate', 'not_checked'];
  const ACTION_LABEL = {run_dft: 'Run DFT', trust: 'Trust', hold: 'Hold (not verified)', regenerate: 'Re-generate', not_checked: 'Fix the input'};
  const VERIFY = 'After each TS optimisation, check for exactly one imaginary frequency and then follow the reaction path ' +
    '(IRC) to confirm it joins this reactant to this product. One imaginary frequency alone is not enough.';

  function num(x) {
    if (x === null || x === undefined || x === '') return null;
    const v = typeof x === 'number' ? x : Number(String(x).trim());
    return Number.isFinite(v) ? v : null;
  }
  function truthy(x) { return x === true || ['true', '1', 'yes'].includes(String(x).trim().toLowerCase()); }
  function cmp(a, b) { return a < b ? -1 : a > b ? 1 : 0; }
  // Python's f'{x:.3f}': round half to even when x sits exactly half way (x.toFixed rounds such ties up).
  function f3(x) {
    const t = x * 16;                 // a double is exactly half way between two 3-decimal numbers only if it is an odd number of sixteenths
    if (Number.isInteger(t) && Math.abs(t) % 2 === 1) { const lo = Math.floor(x * 1000); return ((lo % 2 === 0 ? lo : lo + 1) / 1000).toFixed(3); }
    return x.toFixed(3);
  }

  // ---- SHA-256 (synchronous, so the audit draw matches Python's hashlib on every page, with or without https)
  const K = new Uint32Array([0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]);
  function sha256hex(text) {
    const bytes = Array.from(new TextEncoder().encode(text));
    const bits = bytes.length * 8;
    bytes.push(0x80);
    while (bytes.length % 64 !== 56) bytes.push(0);
    for (let i = 7; i >= 0; i--) bytes.push(i >= 4 ? Math.floor(bits / 2 ** 32) >>> ((i - 4) * 8) & 255 : bits >>> (i * 8) & 255);
    const H = new Uint32Array([0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]);
    const w = new Uint32Array(64), r = (x, n) => (x >>> n) | (x << (32 - n));
    for (let o = 0; o < bytes.length; o += 64) {
      for (let i = 0; i < 16; i++) w[i] = (bytes[o + 4 * i] << 24) | (bytes[o + 4 * i + 1] << 16) | (bytes[o + 4 * i + 2] << 8) | bytes[o + 4 * i + 3];
      for (let i = 16; i < 64; i++) {
        const s0 = r(w[i - 15], 7) ^ r(w[i - 15], 18) ^ (w[i - 15] >>> 3), s1 = r(w[i - 2], 17) ^ r(w[i - 2], 19) ^ (w[i - 2] >>> 10);
        w[i] = (w[i - 16] + s0 + w[i - 7] + s1) | 0;
      }
      let [a, b, c, d, e, f, g, h] = H;
      for (let i = 0; i < 64; i++) {
        const t1 = (h + (r(e, 6) ^ r(e, 11) ^ r(e, 25)) + ((e & f) ^ (~e & g)) + K[i] + w[i]) | 0;
        const t2 = ((r(a, 2) ^ r(a, 13) ^ r(a, 22)) + ((a & b) ^ (a & c) ^ (b & c))) | 0;
        h = g; g = f; f = e; e = (d + t1) | 0; d = c; c = b; b = a; a = (t1 + t2) | 0;
      }
      H[0] += a; H[1] += b; H[2] += c; H[3] += d; H[4] += e; H[5] += f; H[6] += g; H[7] += h;
    }
    return Array.from(H, x => x.toString(16).padStart(8, '0')).join('');
  }
  function auditOrder(ids, seed) {
    return ids.map(i => [sha256hex(seed + ':' + i), i]).sort((p, q) => cmp(p[0], q[0])).map(p => p[1]);
  }

  // ---- the exact audit bound (hypergeometric), with whole-number arithmetic
  function comb(n, k) {
    if (k < 0 || k > n) return 0n;
    k = Math.min(k, n - k);
    let c = 1n;
    for (let i = 1; i <= k; i++) c = c * BigInt(n - k + i) / BigInt(i);
    return c;
  }
  function upperBoundWrong(population, sampled, wrongFound) {      // 95% confidence
    const N = population, k = sampled, x = wrongFound;
    if (!(0 <= x && x <= k && k <= N)) throw new Error('need 0 <= wrong_found <= sampled <= population');
    const total = comb(N, k);
    let best = x;
    for (let M = x; M <= N - (k - x); M++) {
      let tail = 0n;
      for (let i = 0; i <= x; i++) tail += comb(M, i) * comb(N - M, k - i);
      if (tail * 20n > total) best = M; else break;                  // tail / total > 0.05
    }
    return best - x;
  }
  function auditsNeeded(population, target) {
    for (let k = 1; k < population; k++) if (upperBoundWrong(population, k, 0) <= target * (population - k)) return k;
    return null;
  }

  // ---- the plan
  function plan(rows, budget, opts) {
    const o = Object.assign({goal: 'find_errors', audit: 0, seed: 0}, opts || {});
    if (!Number.isInteger(budget) || budget < 0) return {error: 'budget must be a whole number of DFT runs, 0 or more'};
    if (!Number.isInteger(o.audit) || o.audit < 0 || o.audit > budget) return {error: 'audit must be a whole number between 0 and the budget'};
    if (!['find_errors', 'finish_ts'].includes(o.goal)) return {error: "goal must be 'find_errors' or 'finish_ts'"};
    const ids = rows.map(r => String(r.id));
    if (new Set(ids).size !== ids.length) return {error: 'reaction ids must be unique'};
    const items = [], ranked = [], single = [];
    for (const r of rows) {
      const v = String(r.verdict || '').trim().toLowerCase();
      const it = {id: String(r.id), verdict: v || null, disagreement_A: num(r.disagreement_A), endpoint_distance_A: num(r.endpoint_distance_A),
                  mirror_repaired: truthy(r.mirror_repaired), action: null, dft_order: null, why: ''};
      if (!['green', 'amber', 'red'].includes(v)) { it.action = 'not_checked'; it.why = String(r.error || 'no verdict for this reaction').trim(); }
      else if (v === 'red') {
        const broke = String(r.rules_broken || '').trim();
        it.action = 'regenerate'; it.why = 'still breaks a rule after repair' + (broke ? ` (${broke.slice(0, 120)})` : '');
      } else if (it.disagreement_A === null) {
        it.action = 'trust'; it.why = 'passes every rule, but only one guess was given, so model disagreement is unknown';
        if (it.endpoint_distance_A !== null) single.push(it);
      } else ranked.push(it);
      items.push(it);
    }
    const targeted = budget - o.audit;
    let queue;
    if (o.goal === 'find_errors') {
      ranked.sort((a, b) => cmp(-a.disagreement_A, -b.disagreement_A) || cmp(a.id, b.id));
      single.sort((a, b) => cmp(-a.endpoint_distance_A, -b.endpoint_distance_A) || cmp(a.id, b.id));
      queue = ranked.concat(single);
    } else queue = ranked.filter(it => it.verdict === 'amber').sort((a, b) => cmp(a.disagreement_A, b.disagreement_A) || cmp(a.id, b.id));
    const picked = queue.slice(0, targeted);
    picked.forEach((it, k) => {
      it.action = 'run_dft'; it.dft_order = k + 1;
      if (it.disagreement_A === null) it.why = `number ${k + 1}: a single guess, ranked by its distance from the reactant and the product (${f3(it.endpoint_distance_A)} A), a weaker sign than model disagreement`;
      else if (o.goal === 'finish_ts') it.why = `number ${k + 1}: the least disputed of the "check with DFT" reactions (${f3(it.disagreement_A)} A), so the run is likely to end in a usable TS`;
      else it.why = `number ${k + 1} by model disagreement (${f3(it.disagreement_A)} A)` + (it.verdict === 'amber' ? '' :
        '; a spot check: every "check with DFT" reaction is already covered, so this run is optional');
    });
    for (const it of ranked) if (it.action === null) {
      const d = f3(it.disagreement_A) + ' A';
      if (it.verdict === 'amber') { it.action = 'hold'; it.why = `the models disagree (${d}) and the budget did not reach it; do not use it unverified`; }
      else { it.action = 'trust'; it.why = `passes every rule and the models agree (${d})`; }
    }
    const trusted = items.filter(it => it.action === 'trust');
    const drawn = new Set(auditOrder(trusted.map(it => it.id), o.seed).slice(0, o.audit));
    let nRun = picked.length;
    for (const it of trusted.slice().sort((a, b) => cmp(a.id, b.id))) if (drawn.has(it.id)) {
      nRun += 1; it.action = 'run_dft'; it.dft_order = nRun; it.audit = true;
      it.why = 'random audit of a trusted reaction: drawn by lot, not because anything looks wrong, so its outcome says something about all the trusted reactions';
    }
    for (const it of items) if (it.mirror_repaired && ['run_dft', 'trust', 'hold'].includes(it.action)) it.why += '; wrong-handed guess repaired by mirroring';
    const order = Object.fromEntries(ACTIONS.map((a, k) => [a, k]));
    items.sort((a, b) => cmp(order[a.action], order[b.action]) || cmp(a.dft_order || 0, b.dft_order || 0) ||
      cmp(-(a.disagreement_A !== null ? a.disagreement_A : -1), -(b.disagreement_A !== null ? b.disagreement_A : -1)) || cmp(a.id, b.id));
    const counts = Object.fromEntries(ACTIONS.map(a => [a, items.filter(it => it.action === a).length]));
    const nAmber = items.filter(it => it.verdict === 'amber').length;
    const trustPop = counts.trust + drawn.size;
    return {budget, goal: o.goal, reactions: items.length, counts, dft_runs_planned: counts.run_dft, of_which_random_audits: drawn.size,
      of_which_optional_spot_checks: items.filter(it => it.action === 'run_dft' && it.verdict !== 'amber' && !it.audit).length,
      budget_unused: budget - counts.run_dft, audit_seed: o.seed, amber: nAmber, trusted_population: trustPop,
      bound_if_all_audits_right: drawn.size ? upperBoundWrong(trustPop, drawn.size, 0) : null,
      run_dft: items.filter(it => it.action === 'run_dft').map(it => it.id), audited: items.filter(it => it.audit).map(it => it.id), plan: items};
  }

  function explain(p) {
    if (p.error) return `Could not make a plan: ${p.error}.`;
    const c = p.counts, n = p.reactions, b = p.budget;
    const runs = p.plan.filter(it => it.action === 'run_dft' && !it.audit), audits = p.plan.filter(it => it.audit);
    const lines = [`You have ${n} reaction${n !== 1 ? 's' : ''} and ${b} DFT run${b !== 1 ? 's' : ''}.`];
    if (runs.length) {
      const vals = runs.filter(it => it.disagreement_A !== null).map(it => it.disagreement_A);
      const spread = !vals.length ? '' : vals.length > 1 ? ` (${f3(vals[0])} A to ${f3(vals[vals.length - 1])} A)` : ` (${f3(vals[0])} A)`;
      const how = p.goal === 'find_errors' ? 'in order of how much the guesses disagree' :
        'the least disputed of the reactions that need a check, so the runs are likely to end in a usable TS';
      lines.push(`Run DFT on ${runs.length}: ${runs.map(it => it.id).join(', ')}, ${how}${spread}.`);
      const spot = p.of_which_optional_spot_checks;
      if (spot) lines.push((spot === 1 ? 'One of these is an optional spot check on a reaction' : `${spot} of these are optional spot checks on reactions`) +
        ' that would otherwise be trusted: every reaction marked "check with DFT" is already covered.');
    } else if (b && !audits.length) lines.push('No reaction could be ranked, so no run is planned.');
    if (audits.length) lines.push(`Also run DFT on ${audits.length} trusted reaction${audits.length !== 1 ? 's' : ''} drawn at random (${audits.map(it => it.id).join(', ')}): an audit, which is the only way to see an error that all the models share.`);
    if (p.budget_unused > 0 && (runs.length || audits.length)) lines.push(`${p.budget_unused} run(s) are left over: nothing else needs one.`);
    if (c.trust) lines.push(`Trust ${c.trust}: they pass every rule and their guesses agree.`);
    if (c.hold) lines.push(`Hold ${c.hold}: marked "check with DFT" but beyond this budget, so do not use them unverified.`);
    if (c.regenerate) lines.push(`Re-generate ${c.regenerate}: they still break a rule, and a DFT run from such a guess is unlikely to reach the right TS.`);
    if (c.not_checked) lines.push(`${c.not_checked} could not be checked; fix their input files first.`);
    lines.push(VERIFY);
    return lines.join(' ');
  }

  // outcomes: {id: 'right' | 'wrong'} for reactions the plan sent to DFT
  function review(p, outcomes) {
    const sent = p.plan.filter(it => it.action === 'run_dft');
    const targeted = sent.filter(it => !it.audit).map(it => it.id), audits = sent.filter(it => it.audit).map(it => it.id);
    const doneT = targeted.filter(i => outcomes[i]), doneA = audits.filter(i => outcomes[i]);
    const out = {runs_reported: doneT.length + doneA.length,
      targeted: {reported: doneT.length, right: doneT.filter(i => outcomes[i] === 'right').length, wrong: doneT.filter(i => outcomes[i] === 'wrong').length,
                 to_regenerate: doneT.filter(i => outcomes[i] === 'wrong')},
      audits: audits.length, audits_reported: doneA.length, audits_wrong: doneA.filter(i => outcomes[i] === 'wrong').length, not_run: p.counts.trust, at_most_wrong: null};
    if (audits.length && doneA.length === audits.length) out.at_most_wrong = upperBoundWrong(p.counts.trust + audits.length, audits.length, out.audits_wrong);
    return out;
  }

  function parseCsv(text) {                                   // small CSV reader: quoted fields, commas and new lines inside quotes
    const rows = []; let row = [], cell = '', q = false;
    text = text.replace(/^﻿/, '');
    for (let i = 0; i < text.length; i++) {
      const ch = text[i];
      if (q) { if (ch === '"') { if (text[i + 1] === '"') { cell += '"'; i++; } else q = false; } else cell += ch; }
      else if (ch === '"') q = true;
      else if (ch === ',') { row.push(cell); cell = ''; }
      else if (ch === '\n' || ch === '\r') { if (ch === '\r' && text[i + 1] === '\n') i++; row.push(cell); cell = ''; if (row.some(c => c !== '')) rows.push(row); row = []; }
      else cell += ch;
    }
    row.push(cell); if (row.some(c => c !== '')) rows.push(row);
    if (!rows.length) return [];
    const head = rows[0].map(h => h.trim().toLowerCase());
    return rows.slice(1).map((r, k) => {
      const g = name => { const j = head.indexOf(name); return j < 0 ? '' : (r[j] || '').trim(); };
      return {id: g('id') || g('reaction_id') || `row_${k + 1}`, verdict: g('verdict'), disagreement_A: num(g('disagreement_a')),
              endpoint_distance_A: num(g('endpoint_distance_a')), mirror_repaired: truthy(g('mirror_repaired') || g('mirrored')),
              rules_broken: g('rules_broken') || g('failed_rules'), error: g('error')};
    });
  }

  // ---- the agent's answers: values are copied in from tool results by reference, as in tsso/agent.py
  const REF = /\{\{\s*(\d+)\s*:\s*([^}|\s]+)\s*(\|\s*%)?\s*\}\}/g;
  function show(v, percent) {
    if (Array.isArray(v)) return v.length ? v.map(x => show(x, percent)).join(', ') : 'none';
    if (v === true) return 'yes'; if (v === false) return 'no'; if (v === null || v === undefined) return 'not available';
    if (typeof v === 'number') return percent ? (100 * v).toFixed(1) + '%' : Number.isInteger(v) ? String(v) : String(Number(v.toPrecision(4)));
    if (typeof v === 'object') throw new Error('group');
    return String(v);
  }
  function render(answer, outputs) {
    const filled = [], unresolved = [];
    const text = (answer || '').replace(REF, (ref, n, path, pct) => {
      try {
        if (Number(n) < 1 || Number(n) > outputs.length) throw new Error('call');
        let v = outputs[Number(n) - 1];
        for (const part of path.split('.')) { if (v === null || v === undefined || !(part in Object(v))) throw new Error('path'); v = v[part]; }
        const s = show(v, Boolean(pct)); filled.push({reference: ref, value: s}); return s;
      } catch (e) { unresolved.push(ref); return `[no such value: ${n}:${path}]`; }
    });
    return {text, filled, unresolved};
  }
  const NUM = /(?<![\w.])(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?![\w]|\.\d)/g;
  function typedNumbers(answer, userText) {
    const typed = (answer || '').replace(REF, ' ').replace(/^\s*\d+[.)]\s+/gm, '');
    const said = new Set((userText || '').match(NUM) || []);
    return Array.from(new Set(typed.match(NUM) || [])).filter(t => !said.has(t));
  }

  const api = {ACTIONS, ACTION_LABEL, VERIFY, plan, explain, review, upperBoundWrong, auditsNeeded, sha256hex, auditOrder, parseCsv, render, typedNumbers, f3};
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.Planner = api;
})(typeof self !== 'undefined' ? self : this);
