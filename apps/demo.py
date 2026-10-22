# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo>=0.9", "numpy>=1.26", "rdkit>=2024.3; sys_platform != 'emscripten'"]
# ///
"""TS Second Opinion app. Run inside the repository (marimo run apps/demo.py), or on its own anywhere, e.g. molab:
it then downloads the code and data from github.com/tammy357/ts-second-opinion on first start (about 20 MB).
It also runs entirely in a web browser (WebAssembly): https://tammy357.github.io/ts-second-opinion/app/"""
import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", app_title="TS Second Opinion")


@app.cell
async def _():
    import json, shutil, sys
    from pathlib import Path
    import numpy as np
    import marimo as mo

    import os
    IN_BROWSER = sys.platform == 'emscripten' or bool(os.environ.get('TSSO_SIMULATE_BROWSER'))   # the env var is for tests

    async def _browser_root():
        """In the browser (Pyodide), fetch the code and data bundle published with the GitHub Pages site."""
        import io, zipfile
        from pyodide.http import pyfetch
        dest = Path('/tmp/ts-second-opinion')
        if (dest / 'tsso' / 'checks').is_dir():
            return dest
        urls = []
        try:
            _loc = str(mo.notebook_location())
            if _loc.startswith('http'):
                urls.append(_loc.rstrip('/') + '/bundle.zip')
        except Exception:
            pass
        urls += ['https://tammy357.github.io/ts-second-opinion/app/bundle.zip',
                 'https://raw.githubusercontent.com/tammy357/ts-second-opinion/main/docs/app/bundle.zip']
        with mo.status.spinner(title='Loading the rules and the example data (about 2 MB)...'):
            for _u in urls:
                try:
                    _r = await pyfetch(_u)
                    if _r.ok:
                        zipfile.ZipFile(io.BytesIO(await _r.bytes())).extractall(dest)
                        return dest
                except Exception:
                    continue
        raise RuntimeError('Could not download the app bundle from ' + ', '.join(urls))

    def _find_root():
        try:
            here = Path(__file__).resolve().parent.parent
        except NameError:
            here = Path.cwd()
        if (here / 'tsso' / 'checks').is_dir() and (here / 'data').is_dir():
            return here
        import tarfile, tempfile, urllib.request
        cache = Path(tempfile.gettempdir()) / 'ts-second-opinion-main'
        if not (cache / 'tsso' / 'checks').is_dir():
            url = 'https://codeload.github.com/tammy357/ts-second-opinion/tar.gz/refs/heads/main'
            with mo.status.spinner(title='First start: downloading the code and data from GitHub (about 20 MB)...'):
                tmp = Path(tempfile.mkdtemp())
                archive = tmp / 'repo.tar.gz'
                urllib.request.urlretrieve(url, archive)
                with tarfile.open(archive) as tf:
                    try:
                        tf.extractall(tmp, filter='data')
                    except TypeError:          # Python without extraction filters
                        tf.extractall(tmp)
                top = next(p for p in tmp.iterdir() if p.is_dir())
                if cache.exists():
                    shutil.rmtree(cache)
                shutil.move(str(top), str(cache))
        return cache

    ROOT = (await _browser_root()) if IN_BROWSER else _find_root()
    sys.path.insert(0, str(ROOT))
    from tsso.data import load_reactions
    from tsso.triage import load_checks, triage, ORCA_TEMPLATE
    from tsso.stereo import stereocentre_status
    from tsso.viewer import compare_html
    from tsso.report import report_html
    from tsso.service import build_reaction
    from tsso.chem import charge_and_multiplicity
    from tsso.geometry import dmae
    from tsso.metrics import handedness_correct
    from tsso.sandbox import run_check
    from tsso.llm import available_backends
    return (IN_BROWSER, ORCA_TEMPLATE, ROOT, available_backends, build_reaction, charge_and_multiplicity, compare_html,
            dmae, handedness_correct, json, load_checks, load_reactions, mo, np, report_html, run_check, shutil,
            stereocentre_status, triage)


@app.cell
def _(ROOT, json, load_checks, load_reactions):
    _slim = ROOT / 'data' / 'app_examples.json.gz'       # the browser bundle carries only the worked examples
    by_id = {r.id: r for r in (load_reactions(_slim) if _slim.exists() else load_reactions())}
    checks = load_checks()
    metrics = json.load(open(ROOT / 'results' / 'metrics.json'))
    tri = json.load(open(ROOT / 'results' / 'triage_summary.json'))
    cert = json.load(open(ROOT / 'results' / 'certify.json'))
    explanations = json.load(open(ROOT / 'results' / 'explanations.json'))['examples']
    authoring = [json.loads(l) for f in ('authoring_log_run1.jsonl', 'authoring_log.jsonl')
                 for l in open(ROOT / 'results' / f)]
    EXAMPLES = {
        3509: {'title': 'Every model got the handedness wrong',
               'note': 'All 8 GoFlow guesses put C2 the wrong way round. Their interatomic distances are almost '
                       'perfect, so D-MAE and chirality-blind RMSD call them accurate, but each one is the TS of the '
                       'mirror-image reaction.'},
        5877: {'title': 'Only part of the molecule is flipped',
               'note': 'One of three stereocentres is inverted: the guess is the TS of a diastereomer, whose energy '
                       'really differs. Mirroring would invert the two correct centres, so it cannot be repaired.'},
        1292: {'title': 'The models disagree',
               'note': 'Two guesses pass every rule, but the models disagree a lot about the geometry. This is '
                       'where a DFT check pays off.'},
        2577: {'title': 'A guess you can trust',
               'note': 'Every guess passes every rule, both stereocentres are right and the models agree closely.'},
        1093: {'title': 'A problem in the dataset',
               'note': 'Every guess breaks the reacting-bond rule. Looking closer, the atom-mapped SMILES says the '
                       'N5–O6 bond breaks, but the reactant, TS and product geometries all keep it near 1.4 Å, so this '
                       'RDB7 entry looks inconsistent rather than the model being wrong.'},
    }
    return EXAMPLES, authoring, by_id, cert, checks, explanations, metrics, tri


@app.cell
def _(cert, metrics, mo):
    _fa = metrics['final_answer']
    _c5 = next(d for d in cert['geometry']['conformal'] if d['target'] == 0.05)

    def _tile(value, label):
        return (f"<div style='flex:1;min-width:200px;padding:14px 16px;border:1px solid #e3e2de;border-radius:12px;background:#fff'>"
                f"<div style='font-size:28px;font-weight:750;line-height:1.1'>{value}</div>"
                f"<div style='color:#55534e;font-size:13px;margin-top:6px'>{label}</div></div>")

    mo.Html(
        "<div style='padding:4px 0 6px'>"
        "<div style='font-size:13px;letter-spacing:.05em;text-transform:uppercase;color:#55534e'>Open Scientific Intelligence Hackathon 2026</div>"
        "<h1 style='margin:6px 0 8px;font-size:34px'>TS Second Opinion</h1>"
        "<p style='font-size:17px;margin:0 0 14px;max-width:760px'><b>Check, repair and triage machine-learned transition-state "
        "guesses before you spend DFT.</b> Try an example or upload your own reaction: you get a verdict, the reasons, "
        "a repaired guess you can see in 3D, and a ready ORCA input.</p>"
        "<div style='display:flex;gap:12px;flex-wrap:wrap'>"
        + _tile(f"{100 * metrics['stereo_checks_recall_of_stereo_errors']:.0f}%",
                "of wrong-handed GoFlow guesses caught by two LLM-written rules. Most pass D-MAE, RMSD and even a DFT frequency check.")
        + _tile(f"{100 * _fa['single_model_avg']['stereo_bad_rate']:.0f}% → {100 * _fa['check_aware_plus_mirror']['stereo_bad_rate']:.0f}%",
                "final answers with the wrong handedness, before → after checking and mirror repair, with no DFT")
        + _tile(f"{100 * _c5['mean_share_trusted']:.0f}%",
                f"of reactions trusted at a guaranteed 5% average geometry-error rate (actual {100 * _c5['mean_error_among_trusted']:.1f}%)")
        + "</div></div>")
    return


@app.cell
def _(EXAMPLES, ROOT, mo):
    example = mo.ui.radio(options={v['title']: str(k) for k, v in EXAMPLES.items()},
                          value=next(iter(EXAMPLES.values()))['title'])
    smiles_in = mo.ui.text_area(placeholder='[C:1]([H:4])([H:5])[H:6]...>>...', full_width=True,
                                label='Atom-mapped reaction SMILES (every atom mapped, hydrogens too)')
    reactant_in = mo.ui.file(filetypes=['.xyz'], label='Reactant .xyz')
    product_in = mo.ui.file(filetypes=['.xyz'], label='Product .xyz')
    guesses_in = mo.ui.file(filetypes=['.xyz'], multiple=True, label='TS guess(es) .xyz')
    order_in = mo.ui.dropdown({'by atom-map number': 'map', 'in the order of the SMILES': 'smiles'},
                              value='by atom-map number', label='Atoms in my files are listed')
    check_btn = mo.ui.run_button(label='Check my TS guesses', kind='success')
    _d = ROOT / 'examples' / 'rdb7_3509'
    _files = mo.hstack([mo.download(data=(_d / f).read_bytes(), filename=f, label=f)
                        for f in ('reaction.smi', 'reactant.xyz', 'product.xyz', 'guesses.xyz')], justify='start')
    mode = mo.ui.tabs({
        'Try an example': mo.vstack([
            mo.md('Each case is a real GoFlow prediction (8 models) for a reaction in the RDB7 test set.'), example]),
        'Check your own reaction': mo.vstack([
            mo.md('Works with any TS generator and any elements. Give one or more guesses: several (models, seeds or '
                  'samples) add the model-disagreement signal. Charge and multiplicity are read from the SMILES.'),
            smiles_in, mo.hstack([reactant_in, product_in, guesses_in], justify='start'), order_in, check_btn,
            mo.accordion({'No files to hand? Download a ready example and upload it here': _files})]),
    })
    mo.vstack([mo.md('## 1 · Check a TS guess'), mode])
    return check_btn, example, guesses_in, mode, order_in, product_in, reactant_in, smiles_in


@app.cell
def _(build_reaction, by_id, charge_and_multiplicity, check_btn, checks, example, guesses_in, mo, mode, order_in,
      product_in, reactant_in, smiles_in, tri, triage):
    own = mode.value == 'Check your own reaction'
    if own:
        mo.stop(not check_btn.value, mo.callout(mo.md('Add the SMILES and the `.xyz` files, then press '
                                                      '**Check my TS guesses**.'), kind='neutral'))
        if not (smiles_in.value.strip() and reactant_in.value and product_in.value and guesses_in.value):
            _rxn, _err = None, 'Please add the SMILES, the reactant, the product and at least one TS guess.'
        else:
            _rxn, _err = build_reaction(smiles_in.value.strip(), reactant_in.contents().decode(),
                                        product_in.contents().decode(),
                                        [f.contents.decode() for f in guesses_in.value], order_in.value)
        if _err is not None:
            import html as _html
            _msg = _html.escape(_err[0].upper() + _err[1:])
            mo.stop(True, mo.callout(mo.Html(f'<b>Could not check this reaction.</b> {_msg}'), kind='danger'))
        rxn, ex, title = _rxn, None, 'Your reaction'
        charge, mult = charge_and_multiplicity(rxn.smiles)
    else:
        ex = int(example.value)
        rxn, title, charge, mult = by_id[ex], f'RDB7 reaction {ex}', 0, 1
    t = triage(rxn, checks, tri['amber_cutoff_A'])
    return charge, ex, mult, rxn, t, title


@app.cell
def _(EXAMPLES, charge, ex, mo, mult, np, rxn, t, title):
    _col = {'green': ('#137333', '#e6f4ea', '✓', 'Trust', 'Passes every rule and the models agree.'),
            'amber': ('#9a5b00', '#fdf1dc', '!', 'Check with DFT', 'Passes the rules, but the models disagree a lot.'),
            'red': ('#b3261e', '#fbe9e7', '✕', 'Re-generate', 'Still breaks a rule after repair.')}[t.verdict]
    _n = len(rxn.preds)
    _facts = [f'{_n - len(t.rejected_models)} of {_n} raw guesses pass every rule',
              'mirror-image repair applied' if t.mirrored else ('no repair needed' if t.verdict != 'red' else 'cannot be repaired'),
              'one guess: no disagreement signal' if np.isnan(t.disagreement) else f'model disagreement {t.disagreement:.3f} Å',
              f'charge {charge}, multiplicity {mult}']
    _note = (f"<p style='margin:10px 0 0;color:#33322e'>{EXAMPLES[ex]['note']}</p>" if ex else '')
    mo.Html(
        f"<div style='border:2px solid {_col[0]};border-radius:14px;padding:16px 18px;background:#fff;margin-top:6px'>"
        f"<div style='display:flex;gap:16px;align-items:center'>"
        f"<div style='flex:none;width:64px;height:64px;border-radius:14px;background:{_col[1]};color:{_col[0]};"
        f"display:grid;place-items:center;font-size:32px;font-weight:800'>{_col[2]}</div>"
        f"<div><div style='color:#55534e;font-size:13px'>{title}</div>"
        f"<div style='font-size:26px;font-weight:750;color:{_col[0]}'>{_col[3]}</div>"
        f"<div style='color:#55534e'>{_col[4]}</div></div></div>"
        f"<div style='display:flex;flex-wrap:wrap;gap:8px;margin-top:12px'>"
        + ''.join(f"<span style='background:#f1f0ec;border-radius:999px;padding:3px 10px;font-size:13px;color:#55534e'>{f}</span>" for f in _facts)
        + f"</div>{_note}</div>")
    return


@app.cell
def _(mo):
    reveal = mo.ui.switch(label='Overlay the reference DFT TS (examples only; never used for the verdict)', value=False)
    return (reveal,)


@app.cell
def _(compare_html, ex, mo, reveal, rxn, stereocentre_status, t):
    _raw = rxn.preds[t.chosen_model]
    _html = compare_html(rxn, _raw, t.xyz, stereocentre_status(rxn, _raw), stereocentre_status(rxn, t.xyz),
                         labels=(f'Best raw guess (guess {t.chosen_model + 1})',
                                 'Recommended (mirror-repaired)' if t.mirrored else 'Recommended'),
                         reference=rxn.ts_ref_xyz if (ex and reveal.value) else None)
    mo.vstack([mo.md('### Why'), mo.md('\n'.join(f'- {r}' for r in t.reasons)),
               mo.md('### Raw guess and recommended guess'), mo.iframe(_html, height='400px'),
               reveal if ex else mo.md('')])
    return


@app.cell
def _(checks, mo, run_check, rxn, stereocentre_status, t):
    _mark = {True: "<span style='color:#137333;font-weight:600'>✓ right</span>",
             False: "<span style='color:#b3261e;font-weight:600'>✕ inverted</span>",
             None: "<span style='color:#55534e'>too flat to judge</span>"}
    _st_raw, _st_rec = stereocentre_status(rxn, rxn.preds[t.chosen_model]), stereocentre_status(rxn, t.xyz)
    _srows = ''.join(f"<tr><td>{a['label']}</td><td>{a['compared_with']}</td><td>{_mark[a['ok']]}</td><td>{_mark[b['ok']]}</td></tr>"
                     for a, b in zip(_st_raw, _st_rec)) or "<tr><td colspan=4 style='color:#55534e'>No stereocentres in this reaction.</td></tr>"
    _rrows = ''
    for _name, _fn in checks.items():
        _ok, _detail = run_check(_fn, rxn, t.xyz)
        _res = {True: "<span style='color:#137333;font-weight:600'>✓ passes</span>",
                False: "<span style='color:#b3261e;font-weight:600'>✕ fails</span>",
                None: "<span style='color:#55534e'>could not run</span>"}[_ok]
        if _ok and str(_detail).startswith('rule not applicable'):
            _res, _detail = "<span style='color:#55534e'>– n/a</span>", 'nothing in this reaction for this rule to check'
        _rrows += f"<tr><td><code>{_name}</code></td><td>{_res}</td><td>{_detail}</td></tr>"
    _tab = "style='width:100%;border-collapse:collapse;font-size:14px'"
    _th = lambda *c: '<tr>' + ''.join(f"<th style='text-align:left;color:#55534e;font-weight:600;padding:6px 8px;border-bottom:1px solid #e3e2de'>{x}</th>" for x in c) + '</tr>'
    mo.accordion({
        'Stereocentres, raw vs recommended': mo.Html(f"<table {_tab}>{_th('Atom', 'Compared with', 'Raw guess', 'Recommended')}{_srows}</table>"),
        'Every rule on the recommended guess': mo.Html(f"<table {_tab}>{_th('Rule', 'Result', 'Detail')}{_rrows}</table>"),
    })
    return


@app.cell
def _(ORCA_TEMPLATE, charge, checks, dmae, ex, handedness_correct, mo, mult, report_html, rxn, t, title):
    _score = ({'dmae_A': dmae(t.xyz, rxn.ts_ref_xyz), 'stereo_correct': handedness_correct(rxn, t.xyz, rxn.ts_ref_xyz)}
              if ex else None)
    _stem = f'rxn{ex}_ts' if ex else 'my_reaction_ts'
    _xyz = '\n'.join([str(rxn.n), f'TS Second Opinion recommended guess ({title})'] +
                     [f'{s} {a:.6f} {b:.6f} {c:.6f}' for s, (a, b, c) in zip(rxn.symbols, t.xyz)]) + '\n'
    _src = f"guess {t.chosen_model + 1}{' (mirror image)' if t.mirrored else ''}; verdict {t.verdict}"
    mo.vstack([mo.md('### Take it with you'), mo.hstack([
        mo.download(data=report_html(rxn, t, checks, charge, mult, title, _score).encode(),
                    filename=f'{_stem}_report.html', label='Report card (.html)'),
        mo.download(data=_xyz.encode(), filename=f'{_stem}.xyz', label='Recommended TS (.xyz)'),
        mo.download(data=ORCA_TEMPLATE.format(rid=title, source=_src, xyzname=f'{_stem}.xyz', charge=charge,
                                              mult=mult).encode(), filename=f'{_stem}.inp',
                    label='ORCA OptTS + Freq input'),
    ], justify='start')])
    return


@app.cell
def _(mo, shutil):
    explain_btn = mo.ui.run_button(label='Ask Claude to explain this verdict', disabled=shutil.which('claude') is None,
                                   tooltip='Needs the Claude Code CLI on the machine running this app')
    return (explain_btn,)


@app.cell
def _(ex, explain_btn, explanations, mo, rxn, t):
    if ex and str(ex) in explanations:
        _out = mo.callout(mo.md("**Claude's explanation** (generated in advance with the Claude Code CLI, "
                                "from the facts above):\n\n" + explanations[str(ex)]['explanation']), kind='info')
    elif explain_btn.value:
        import subprocess as _sp
        _prompt = ('You are helping a computational chemist decide whether to trust a machine-learned transition-state guess. '
                   'In 3 short sentences, plain English, no headings, explain the verdict and what to do next. Use only '
                   'these facts and do not suggest a repair they do not support.\n'
                   f'Reaction (atom-mapped SMILES): {rxn.smiles}\nVerdict: {t.verdict}\n' + '\n'.join(t.reasons))
        with mo.status.spinner(title='Asking Claude...'):
            _res = _sp.run(['claude', '-p'], input=_prompt, capture_output=True, text=True, timeout=180)
        _out = mo.callout(mo.md(_res.stdout.strip() or _res.stderr.strip()), kind='info')
    else:
        _out = explain_btn
    _out
    return


@app.cell
def _(ROOT, mo):
    import re as _re
    _rows = []
    for _p in sorted((ROOT / 'tsso' / 'checks').glob('*.py')):
        _src = _p.read_text()
        _rule = _re.search(r'^RULE = "(.*)"$', _src, _re.M)
        _gate = _re.search(r'passes ([\d.]+)%; rejects (\d+)% of mirror images and (\d+)% of noisy', _src)
        _rows.append({'rule': _rule.group(1) if _rule else _p.stem, 'real TSs accepted': f'{_gate.group(1)}%' if _gate else '',
                      'mirror images caught': f'{_gate.group(2)}%' if _gate else '', 'noisy TSs caught': f'{_gate.group(3)}%' if _gate else ''})
    mo.vstack([
        mo.md('## 2 · Rule Studio: plain English in, tested code out'),
        mo.md('A chemist writes a rule as a sentence. Claude writes the Python, and the **gate** keeps it only if it '
              'accepts at least 99% of 9,540 real DFT transition states, never crashes, and catches at least 1% of '
              'mirrored or distorted ones. Otherwise Claude gets the failures back and rewrites it. The six rules in use:'),
        mo.ui.table(_rows, selection=None, page_size=10)])
    return


@app.cell
def _(authoring, mo):
    _runs = {f"{e['name']}: {len(e['attempts'])} attempts": i for i, e in enumerate(authoring) if len(e['attempts']) > 1}
    replay = mo.ui.dropdown(_runs, value='h_transfer: 2 attempts' if 'h_transfer: 2 attempts' in _runs else next(iter(_runs)),
                            label='Replay a recorded run where the first attempt failed')
    replay
    return (replay,)


@app.cell
def _(authoring, mo, replay):
    _e = authoring[replay.value]
    _cards = []
    for _a in _e['attempts']:
        if _a['accepted']:
            _verdict, _col = 'Accepted: kept', '#137333'
        elif _a['crash_rate'] > 0:
            _verdict, _col = f"Rejected: crashed on {100 * _a['crash_rate']:.2f}% of structures", '#b3261e'
        elif _a['pass_rate'] < 0.99:
            _verdict, _col = f"Rejected: too strict, it called {100 * (1 - _a['pass_rate']):.1f}% of real TSs wrong", '#b3261e'
        else:
            _verdict, _col = 'Rejected: it almost never fires, so it would add nothing', '#b3261e'
        _cards.append(mo.vstack([
            mo.Html(f"<div style='border-left:4px solid {_col};padding:6px 12px'><b>Attempt {_a['attempt']}</b> · accepts "
                    f"{100 * _a['pass_rate']:.1f}% of real TSs (needs ≥ 99%) · catches {100 * _a['rejects_mirror_images']:.1f}% "
                    f"of mirror images and {100 * _a['rejects_noisy_geometries']:.1f}% of distorted TSs<br>"
                    f"<span style='color:{_col};font-weight:600'>{_verdict}</span></div>"),
            mo.accordion({'Code Claude wrote': mo.md(f"```python\n{_a['code']}\n```")})]))
    mo.vstack([mo.md(f"**Rule:** “{_e['rule']}”")] + _cards)
    return


@app.cell
def _(IN_BROWSER, available_backends, mo):
    new_rule = mo.ui.text_area(placeholder='e.g. Hydrogens that move between atoms sit roughly between them at the transition state.',
                               label='Try your own rule', full_width=True)
    api_key = mo.ui.text(kind='password', placeholder='sk-ant-...', label='Your Anthropic API key')
    write_btn = mo.ui.run_button(label='Write and test it with Claude')
    _b = available_backends()
    if IN_BROWSER:
        _note = ('Uses your own Anthropic API key. It stays in this browser tab and is sent only to api.anthropic.com; '
                 'nothing is stored. Each rule takes up to 3 Claude calls (a few cents). The 9,540 test structures (6 MB) '
                 'download the first time. Takes about 30–90 s.')
        _ui = mo.vstack([new_rule, api_key, write_btn, mo.md(f"<span style='color:#55534e;font-size:12px'>{_note}</span>")])
    else:
        _note = (f"LLM route: {_b[0]}{' (no Claude CLI or API key here: run the repository locally to write new rules)' if _b[0] == 'cache' else ''}"
                 '. Takes about 15–60 s.')
        _ui = mo.vstack([new_rule, mo.hstack([write_btn, mo.md(f"<span style='color:#55534e;font-size:12px'>{_note}</span>")],
                                             justify='start')])
    _ui
    return api_key, new_rule, write_btn


@app.cell
async def _(IN_BROWSER, ROOT, api_key, mo, new_rule, write_btn):
    mo.stop(not write_btn.value or not new_rule.value.strip())
    import html as _h
    import re as _re2
    from tsso.author import author, author_async
    _name = 'user_' + _re2.sub(r'[^a-z0-9]+', '_', new_rule.value.lower()).strip('_')[:30]
    _save = ROOT / 'tsso' / 'checks_user'

    def _fail(msg):
        mo.stop(True, mo.callout(mo.Html(f'<b>Could not write the rule.</b> {_h.escape(msg)}'), kind='danger'))

    if IN_BROWSER:
        if not api_key.value.strip():
            _fail('Paste your Anthropic API key first (console.anthropic.com → API keys).')
        from pyodide.http import pyfetch
        from tsso.llm import ask_async
        _train = ROOT / 'data' / 'rdb7_train_reference.json.gz'
        if not _train.exists():
            with mo.status.spinner(title='Downloading the 9,540 reference transition states (6 MB, first time only)...'):
                _r = await pyfetch('https://raw.githubusercontent.com/tammy357/ts-second-opinion/main/data/rdb7_train_reference.json.gz')
                if not _r.ok:
                    _fail('Could not download the reference transition states from GitHub.')
                _train.parent.mkdir(parents=True, exist_ok=True)
                _train.write_bytes(await _r.bytes())
        _key = api_key.value.strip()
        try:
            with mo.status.spinner(title='Claude is writing the check; testing it on 9,540 real transition states...') as _sp:
                _out = await author_async(
                    new_rule.value.strip(), lambda _p: ask_async(_p, _key, pyfetch), name=_name, save_dir=_save,
                    on_attempt=lambda _a: _sp.update(subtitle=f"Attempt {_a['attempt']}: accepts {100 * _a['pass_rate']:.1f}% "
                                                              f"of real TSs{'; kept' if _a['accepted'] else '; asking Claude to revise'}"))
        except Exception as _e:
            _fail(str(_e))
    else:
        with mo.status.spinner(title='Claude is writing the check; testing it on 9,540 real transition states...'):
            _out = author(new_rule.value.strip(), name=_name, verbose=False, save_dir=_save)
    _s = _out['stats']
    _steps = [mo.Html(f"<div style='border-left:4px solid {'#137333' if _a['accepted'] else '#b3261e'};padding:4px 12px'>"
                      f"<b>Attempt {_a['attempt']}</b> · accepts {100 * _a['pass_rate']:.1f}% of real TSs (needs ≥ 99%) · catches "
                      f"{100 * _a['rejects_mirror_images']:.1f}% of mirror images and {100 * _a['rejects_noisy_geometries']:.1f}% of "
                      f"distorted TSs · <b>{'kept' if _a['accepted'] else 'rejected'}</b></div>")
              for _a in _out.get('attempt_log', [])]
    mo.vstack([
        mo.md(f"**{'✓ Accepted' if _out['accepted'] else '✕ Rejected'}** after {_out['attempts']} attempt(s): accepts "
              f"{100 * _s['pass_rate']:.1f}% of real TSs, catches {100 * _s['rejects_mirror_images']:.1f}% of mirror images "
              f"and {100 * _s['rejects_noisy_geometries']:.1f}% of distorted ones."
              + ('' if _out['accepted'] else ' The gate did its job: this rule is not used.'))]
        + _steps + [mo.accordion({'Code Claude wrote': mo.md(f"```python\n{_out['source']}\n```")})])
    return


@app.cell
def _(mo):
    target = mo.ui.slider(steps=[2, 3, 4, 5, 6, 8, 10, 12, 15, 20], value=5, label='Target error rate (%)', show_value=True)
    count = mo.ui.radio({'geometry errors (what matters for the barrier)': 'geometry',
                         'geometry or handedness': 'geometry_or_handedness'},
                        value='geometry errors (what matters for the barrier)', label='Count as wrong')
    mo.vstack([mo.md('## 3 · How much can you trust?'),
               mo.md('Conformal selection turns the triage into a guarantee: choose an error rate, and it picks which '
                     'reactions to trust so that, on average, at most that share of them are wrong. Calibrated on half of '
                     'the 1,194 test reactions, tested on the other half, 2,000 times.'),
               mo.hstack([target, count], justify='start', gap=2)])
    return count, target


@app.cell
def _(cert, count, mo, target):
    _d = next(d for d in cert[count.value]['conformal'] if abs(d['target'] - target.value / 100) < 1e-9)
    mo.Html(f"<div style='border:1px solid #e3e2de;border-radius:12px;padding:16px 18px;background:#fff'>"
            f"<div style='font-size:15px'>At a <b>{target.value}%</b> target, you can trust</div>"
            f"<div style='font-size:40px;font-weight:750;line-height:1.1;margin:4px 0'>{100 * _d['mean_share_trusted']:.0f}% of reactions</div>"
            f"<div style='color:#55534e'>and on average {100 * _d['mean_error_among_trusted']:.1f}% of those were wrong "
            f"(in 95% of splits at most {100 * _d['error_among_trusted_95th_percentile']:.1f}%). The guarantee assumes new "
            f"reactions look like the calibration ones; recalibrate for a new generator or new chemistry "
            f"(<code>tsso certify</code>, or <code>tsso.certify</code> in Python).</div></div>")
    return


@app.cell
def _(ROOT, json, mo):
    _rows = json.load(open(ROOT / 'results' / 'leaderboard.json'))['rows']
    _pc = lambda v: f'{100 * v:.1f}%'
    _td = "style='padding:6px 8px;border-bottom:1px solid #e3e2de;text-align:{}'"
    _body = ''
    for _r in _rows:
        _ours = 'TS Second Opinion' in _r['method']
        _bg = 'background:#f3f8f4;' if _ours else ''
        _sd = _r.get('sd')
        _cell = lambda k: _pc(_r[k]) + (f" ± {_pc(_sd[k])}" if _sd else '')
        _body += (f"<tr style='{_bg}'><td {_td.format('left')}>{'<b>' if _ours else ''}{_r['method']}{'</b>' if _ours else ''}</td>"
                  f"<td {_td.format('left')}>{_r['test_set']}</td><td {_td.format('right')}>{_cell('success_dmae_only')}</td>"
                  f"<td {_td.format('right')}>{_cell('handedness_correct')}</td>"
                  f"<td {_td.format('right')}><b>{_cell('success_dmae_and_handedness')}</b></td></tr>")
    _th = "style='padding:6px 8px;border-bottom:2px solid #e3e2de;color:#55534e;text-align:{}'"
    mo.vstack([mo.md('## 4 · Stereo-aware leaderboard'),
               mo.md('A guess usually counts as a success if its D-MAE is at most 0.2 Å, which cannot see handedness. '
                     'The last column counts only guesses that are close **and** correctly handed (`tsso evaluate`).'),
               mo.Html(f"<table style='width:100%;border-collapse:collapse;font-size:14px'><tr><th {_th.format('left')}>Method</th>"
                       f"<th {_th.format('left')}>Test set</th><th {_th.format('right')}>Success (D-MAE ≤ 0.2 Å)</th>"
                       f"<th {_th.format('right')}>Correct handedness</th><th {_th.format('right')}>Success with correct handedness</th></tr>"
                       f"{_body}</table>"),
               mo.md('Handedness stays at 82–85% correct for every model size, so a bigger model does not fix it. '
                     'Add your model with `tsso evaluate` (see CONTRIBUTING.md).')])
    return


@app.cell
def _(ROOT, mo):
    mo.vstack([
        mo.md('## 5 · Results'),
        mo.hstack([mo.image((ROOT / 'results' / 'fig_final_answer.png').read_bytes(), width=470),
                   mo.image((ROOT / 'results' / 'fig_certify.png').read_bytes(), width=560)], widths='equal'),
        mo.image((ROOT / 'results' / 'fig_risk_coverage.png').read_bytes(), width=760),
        mo.md("<span style='color:#55534e;font-size:12px'>Retrospective: GoFlow's released predictions for the 1,194 RDB7 "
              "test reactions; rules were gated only on the 9,540 training reactions. The GoFlow paper notes that D-MAE "
              "misses chirality errors, and FK-Flow (Mark, Galustian, Kovar & Heid, 2025) measured them and fixed them by "
              "steering generation. TS Second Opinion works after generation, on any model's output, and adds repair, a "
              "statistical guarantee and agent tools. Code: github.com/tammy357/ts-second-opinion (MIT).</span>")])
    return


if __name__ == "__main__":
    app.run()
