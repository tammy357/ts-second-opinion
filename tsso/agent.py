"""A small agent that plans a DFT budget in conversation: Claude chooses which tool to call and explains the result;
the tools make every decision and produce every number.

Three things keep it honest, and each is code rather than a promise:

  1. Tools decide. Claude can only call the functions in a tool registry (by default the planner and a few read-only
     tools from tsso.service). The plan is a fixed function of its inputs (tsso/plan.py).
  2. No structures leave. Tool results are passed through `redact` before Claude sees them: coordinates, SMILES and
     ORCA inputs are removed. Verdicts, rule names, atom labels such as "C3-H7" and disagreement values are sent.
  3. Numbers are traced. `check_numbers` lists every number in Claude's answer that does not appear in a tool result
     or in the person's own messages. It is a flag for the reader, not proof: a traced number can still sit in a
     wrong sentence, and a number Claude worked out (a sum, a difference) is flagged even when right.

The loop is written without a network client so it runs anywhere: `request` builds the next API call and `absorb`
takes the reply, runs any tools and says whether the answer is finished. `run` drives it with any function that
posts JSON to the Anthropic Messages API (`http_sender` uses only the standard library); `run_async` does the same
with an async sender, which is how the in-browser app calls it with pyodide's fetch.

    from tsso.agent import run, http_sender
    out = run('I have 40 reactions and 5 DFT runs. Where should they go?', http_sender(my_key))
    print(out['answer'], out['numbers_not_from_tools'])
"""
import json
import re

from .llm import API_MODEL, API_URL, request_headers

MAX_STEPS = 8
SYSTEM = (
    'You help a computational chemist decide where to spend a fixed number of DFT runs on machine-learned '
    'transition-state (TS) guesses, using the TS Second Opinion tools.\n'
    'Rules you must follow:\n'
    '- Call a tool for every decision and every number. Never choose which reactions get DFT yourself, never estimate '
    'an error rate, and do not do arithmetic on tool output: quote numbers exactly as the tools return them.\n'
    '- If the person has not said how many DFT runs they can afford, ask before planning.\n'
    '- Say what to do next in plain words: which reactions to run, which to trust, which are on hold, which to '
    're-generate. Keep it under 150 words unless asked for more. No headings.\n'
    '- Always pass on the caveats that apply: confirm each optimised TS with a reaction-path (IRC) calculation; model '
    'disagreement cannot see an error all models share; benchmark numbers come from GoFlow guesses on RDB7 and are '
    'not a prediction for other chemistry; a certified error rate is an average that assumes the new reactions '
    'resemble the calibration set.\n'
    '- You cannot see structures (coordinates and SMILES are withheld), so do not comment on a mechanism.\n'
    '- If a tool returns an error, say so plainly and stop; do not guess.')

_STRUCTURE_KEYS = ('recommended_ts_xyz', 'orca_input', 'smiles', 'mapped_reaction_smiles', 'bond_orders_assigned',
                   'reactant_xyz', 'product_xyz', 'guesses_xyz', 'xyz', 'feedback_url', 'inference_notes')


def redact(x):
    """A copy of a tool result without structures: no coordinates, SMILES or ORCA inputs."""
    if isinstance(x, dict):
        return {k: redact(v) for k, v in x.items() if k not in _STRUCTURE_KEYS}
    if isinstance(x, (list, tuple)):
        return [redact(v) for v in x]
    return x


# ------------------------------------------------------------------ tools
def _tool(name, description, properties, required=()):
    return {'name': name, 'description': description,
            'input_schema': {'type': 'object', 'properties': properties, 'required': list(required)}}


def default_tools():
    """{name: (schema for the API, function)}: the planner plus read-only tools from tsso.service."""
    from . import service
    return {
        'plan_dft_budget': (_tool(
            'plan_dft_budget',
            'Plan where a fixed number of DFT runs goes. Give path (a verdicts.csv from `tsso batch`, the folder that '
            'holds it, or a folder or CSV of reactions to triage first), or leave path empty to use the RDB7 benchmark: '
            'then give reaction_ids, or sample (how many reactions to draw at random) and seed. evaluate=true (benchmark '
            'only) adds what the reference TS says about the picks; that is for evaluation, never for deciding. '
            'target_error_rate adds certified trust (benchmark only, calibrated on the reactions not in the batch).',
            {'budget': {'type': 'integer', 'minimum': 0, 'description': 'DFT runs the person can afford'},
             'path': {'type': 'string'}, 'reaction_ids': {'type': 'array', 'items': {'type': 'integer'}},
             'sample': {'type': 'integer', 'minimum': 2}, 'seed': {'type': 'integer'},
             'evaluate': {'type': 'boolean'}, 'target_error_rate': {'type': 'number'}}, ['budget']),
            lambda **a: service.plan_dft_budget(**a)),
        'plan_evidence': (_tool(
            'plan_evidence',
            'How well the planning rule did when tested against reference TSs on RDB7 (GoFlow guesses): wrong guesses '
            'among the picks against picking at random, for several budgets and batch sizes, a second ensemble, what is '
            'left unchecked, what TS optimisations returned, and the limits of the test.', {}),
            lambda: service.plan_evidence()),
        'triage_reaction': (_tool(
            'triage_reaction', 'Verdict and reasons for one RDB7 test reaction (benchmark only).',
            {'reaction_id': {'type': 'integer'}}, ['reaction_id']),
            lambda reaction_id: service.triage_reaction(reaction_id)),
        'benchmark_summary': (_tool(
            'benchmark_summary', 'Headline benchmark results: geometry and handedness error rates before and after, '
            'certified-trust numbers.', {}),
            lambda: service.benchmark_summary()['headline']),
    }


def batch_tools(rows, evidence=None):
    """Tools tied to one batch of verdict rows (what the app uses): Claude can re-plan this batch for another budget
    and read the benchmark evidence, and nothing else. The rows never include reference answers."""
    from .plan import plan_dft_budget, explain_plan
    keep = ('id', 'verdict', 'disagreement_A', 'mirror_repaired', 'rules_broken', 'error')
    rows = [{k: r.get(k) for k in keep} for r in rows]

    def plan(budget):
        p = plan_dft_budget(rows, int(budget) if isinstance(budget, float) and budget.is_integer() else budget)
        if 'error' not in p:
            p['in_plain_words'] = explain_plan(p)
            p['plan'] = [it for it in p['plan'] if it['action'] != 'trust']      # the trusted ones are counted, not listed
        return p

    def evid():
        if not evidence:
            return {'error': 'the benchmark evidence file is not available here'}
        return {k: evidence[k] for k in ('question', 'wrong_means', 'B_batches', 'D_what_is_left', 'E_what_the_runs_return',
                                         'limits') if k in evidence} | {
            'all_1194_reactions_by_budget': evidence.get('A_eight_seed_models_geometry', {}).get('by_budget'),
            'second_ensemble_by_budget': evidence.get('C_five_model_sizes_geometry', {}).get('by_budget')}

    return {
        'plan_dft_budget': (_tool('plan_dft_budget', f'Plan where the DFT runs go for the {len(rows)} reactions in front '
                                  'of the person. A fixed rule: re-generate what breaks a rule, run DFT on the most '
                                  'disputed, hold the other "check with DFT" reactions, trust the rest.',
                                  {'budget': {'type': 'integer', 'minimum': 0, 'description': 'DFT runs the person can afford'}},
                                  ['budget']), plan),
        'plan_evidence': (_tool('plan_evidence', 'How the planning rule did against reference TSs on the RDB7 benchmark '
                                '(GoFlow guesses), with the limits of that test. Benchmark results, not predictions.', {}), evid),
    }


# ------------------------------------------------------------------ the loop, one step at a time
def new_conversation(question, history=None):
    """State for one exchange. history: earlier [{'role', 'content'}] turns with plain-text content."""
    msgs = [{'role': m['role'], 'content': str(m['content'])} for m in (history or [])]
    msgs.append({'role': 'user', 'content': question})
    return {'messages': msgs, 'steps': [], 'answer': None, 'model': None, 'stopped': None}


def request(conv, tools, model=None, max_tokens=1500):
    """The JSON body of the next Messages API call."""
    return {'model': model or API_MODEL, 'max_tokens': max_tokens, 'system': SYSTEM,
            'tools': [schema for schema, _ in tools.values()], 'messages': conv['messages']}


def absorb(conv, reply, tools):
    """Take the API's reply; run any tools it asked for. True when the answer is finished."""
    if reply.get('type') == 'error' or 'content' not in reply:
        raise RuntimeError('the Anthropic API said: ' + str((reply.get('error') or {}).get('message', reply))[:300])
    conv['model'] = reply.get('model', conv['model'])
    content = reply['content']
    conv['messages'].append({'role': 'assistant', 'content': content})
    calls = [b for b in content if b.get('type') == 'tool_use']
    if reply.get('stop_reason') != 'tool_use' or not calls:
        conv['answer'] = ''.join(b.get('text', '') for b in content if b.get('type') == 'text').strip()
        conv['stopped'] = reply.get('stop_reason')
        return True
    results = []
    for b in calls:
        name, args = b.get('name'), b.get('input') or {}
        if name not in tools:
            out = {'error': f'there is no tool called {name}'}
        else:
            try:
                out = redact(tools[name][1](**args))
            except Exception as e:                       # a bad argument or a failed tool: tell the model, keep going
                out = {'error': f'{type(e).__name__}: {e}'}
        conv['steps'].append({'tool': name, 'input': args, 'output': out})
        results.append({'type': 'tool_result', 'tool_use_id': b.get('id'), 'content': json.dumps(out, default=str),
                        'is_error': isinstance(out, dict) and 'error' in out})
    conv['messages'].append({'role': 'user', 'content': results})
    return False


def finish(conv):
    """What a caller shows: the answer, the tool calls behind it, and the numbers that could not be traced."""
    user_text = ' '.join(m['content'] for m in conv['messages'] if m['role'] == 'user' and isinstance(m['content'], str))
    return {'answer': conv['answer'], 'steps': conv['steps'], 'model': conv['model'], 'stopped': conv['stopped'],
            'numbers_not_from_tools': check_numbers(conv['answer'] or '', [s['output'] for s in conv['steps']], user_text)}


def run(question, send, tools=None, history=None, model=None, max_steps=MAX_STEPS):
    """Answer one question. send(body: dict) -> reply dict posts to the Messages API (see http_sender)."""
    tools = default_tools() if tools is None else tools
    conv = new_conversation(question, history)
    for _ in range(max_steps):
        if absorb(conv, send(request(conv, tools, model)), tools):
            return finish(conv)
    conv['answer'], conv['stopped'] = (conv['answer'] or 'Stopped: too many tool calls without an answer.'), 'max_steps'
    return finish(conv)


async def run_async(question, send, tools=None, history=None, model=None, max_steps=MAX_STEPS):
    """Like run, with an async send (for example fetch_sender in the browser)."""
    tools = default_tools() if tools is None else tools
    conv = new_conversation(question, history)
    for _ in range(max_steps):
        if absorb(conv, await send(request(conv, tools, model)), tools):
            return finish(conv)
    conv['answer'], conv['stopped'] = (conv['answer'] or 'Stopped: too many tool calls without an answer.'), 'max_steps'
    return finish(conv)


def http_sender(key, timeout=120):
    """send() for run(): posts to the Anthropic API with your key, using only the standard library."""
    import urllib.error, urllib.request

    def send(body):
        req = urllib.request.Request(API_URL, data=json.dumps(body).encode(), headers=request_headers(key), method='POST')
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:              # the API explains errors in a JSON body
            try:
                return json.loads(e.read().decode())
            except Exception:
                raise RuntimeError(f'the Anthropic API returned HTTP {e.code}') from None
    return send


def fetch_sender(key, fetch):
    """send() for run_async() in the browser: fetch is pyodide.http.pyfetch. The key stays in the visitor's tab."""
    async def send(body):
        resp = await fetch(API_URL, method='POST', headers=request_headers(key), body=json.dumps(body))
        return json.loads(await resp.string())
    return send


# ------------------------------------------------------------------ tracing numbers
_NUM = re.compile(r'(?<![\w.])(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?![\w]|\.\d)')


def _numbers_in_text(text):
    out = []
    for m in _NUM.finditer(text):
        tok = m.group(0)
        try:
            out.append((tok, float(tok.replace(',', ''))))
        except ValueError:
            continue
    return out


def _pool(x, acc):
    if isinstance(x, bool) or x is None:
        return
    if isinstance(x, (int, float)):
        acc.add(float(x))
    elif isinstance(x, str):
        acc.update(v for _, v in _numbers_in_text(x))
        acc.update(float(d) for d in re.findall(r'\d+', x) if len(d) < 12)     # ids such as rxn_1093
    elif isinstance(x, dict):
        for k, v in x.items():
            _pool(str(k), acc); _pool(v, acc)
    elif isinstance(x, (list, tuple)):
        for v in x:
            _pool(v, acc)


def check_numbers(answer, tool_outputs, user_text=''):
    """Numbers in `answer` that appear in no tool output and in none of the person's messages.

    A whole number is traced only by the same whole number in a source (a list of 5 ids also traces "5"). A decimal is
    traced by a source value that rounds to it at the digits shown. A percentage is traced by a source value that
    rounds to it as written or as a fraction (0.117 traces "11.7%" and "12%"). Markdown list numbers are ignored.

    This is a coarse net, not a proof. After the agent has read the benchmark evidence, about a third of arbitrary
    whole numbers below 100 happen to match something (tests/test_plan.py measures it); invented decimals and
    percentages are caught far more often. A traced number can still sit in a wrong sentence.
    """
    text = re.sub(r'(?m)^\s*\d+[.)]\s+', '', answer)
    pool = set()
    for o in tool_outputs:
        _pool(o, pool)
        _lengths(o, pool)
    _pool(user_text, pool)
    missing = []
    for m in _NUM.finditer(text):
        tok = m.group(0)
        value = float(tok.replace(',', ''))
        digits = len(tok.split('.')[1]) if '.' in tok else 0
        percent = text[m.end():m.end() + 2].lstrip().startswith('%')
        if percent:
            ok = any(round(p, digits) == value or round(100 * p, digits) == value for p in pool)
        elif digits == 0:
            ok = value in pool
        else:
            ok = any(round(p, digits) == value for p in pool)
        if not ok and tok not in missing:
            missing.append(tok)
    return missing


def _lengths(x, acc):
    if isinstance(x, dict):
        for v in x.values():
            _lengths(v, acc)
    elif isinstance(x, (list, tuple)):
        acc.add(float(len(x)))
        for v in x:
            _lengths(v, acc)
