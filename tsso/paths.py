"""Where TS Second Opinion finds its bundled data and writes its files.

In a clone of the repository (pip install -e .) everything is local and nothing changes. From an installed package
(pip install git+https://github.com/tammy357/ts-second-opinion, or uvx), checking your own reactions needs no extra
data; the commands that use the RDB7 benchmark (triage of test reactions, the gate, the audit, the mirror test) fetch
the files they need from GitHub on first use (about 11 MB in total) and keep them in a cache folder:
~/.cache/ts-second-opinion, or the folder named by the TSSO_CACHE environment variable. Files that TS Second Opinion
writes go to results/ in a clone and to ./tsso_output otherwise.
"""
import os, shutil, sys, tempfile, urllib.request
from pathlib import Path

from . import __version__

REPO = 'tammy357/ts-second-opinion'
ROOT = Path(__file__).resolve().parent.parent
IN_CLONE = (ROOT / 'pyproject.toml').is_file() and (ROOT / 'tsso' / 'checks').is_dir() and (ROOT / 'data').is_dir()


def cache_dir():
    base = os.environ.get('TSSO_CACHE')
    return Path(base) if base else Path(os.environ.get('XDG_CACHE_HOME') or Path.home() / '.cache') / 'ts-second-opinion'


def fetch(relative_path):
    """Local path to a repository file such as 'data/rdb7_train_reference.json.gz', downloading it once if needed."""
    local = ROOT / relative_path
    if local.exists():
        return local
    target = cache_dir() / __version__ / relative_path
    if target.exists():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    last = None
    print(f'TS Second Opinion: downloading {relative_path} (first use only) ...', file=sys.stderr)
    for ref in (f'v{__version__}', 'main'):       # the files that match this version, else the latest
        url = f'https://raw.githubusercontent.com/{REPO}/{ref}/{relative_path}'
        try:
            with urllib.request.urlopen(url, timeout=60) as resp, \
                    tempfile.NamedTemporaryFile(dir=target.parent, delete=False) as tmp:
                shutil.copyfileobj(resp, tmp)
            Path(tmp.name).replace(target)
            return target
        except Exception as e:                   # try the next ref, then explain
            last = e
    raise FileNotFoundError(f'{relative_path} is not available here and could not be downloaded ({last}). '
                            f'Clone the repository (git clone https://github.com/{REPO}) or put the file in '
                            f'{cache_dir() / __version__}.')


def data_file(name):
    return fetch(f'data/{name}')


def results_file(name):
    """A results file (for example 'triage.csv'), or None if it is neither local nor downloadable."""
    try:
        return fetch(f'results/{name}')
    except FileNotFoundError:
        return None


def output_dir():
    """Folder for files TS Second Opinion writes (ORCA inputs, logs)."""
    if IN_CLONE and os.access(ROOT / 'results', os.W_OK):
        return ROOT / 'results'
    return Path.cwd() / 'tsso_output'


def shown(path):
    """A path as people expect to see it: relative to the repository or the working folder when possible."""
    path = Path(path)
    for base in (ROOT, Path.cwd()):
        try:
            return str(path.resolve().relative_to(base.resolve()))
        except ValueError:
            continue
    return str(path)
