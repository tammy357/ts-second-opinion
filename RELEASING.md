# Releasing

A release puts the package on PyPI (so `pip install ts-second-opinion` and `uvx ts-second-opinion` work) and lists the
MCP server in the official MCP Registry as `io.github.tammy357/ts-second-opinion`, so agents can find it. Both are done
by `.github/workflows/publish.yml`, which runs when a GitHub Release is published. No passwords or tokens are stored in
the repository: PyPI and the MCP Registry trust this repository's GitHub Actions directly.

## One-time setup on PyPI

1. Sign in at [pypi.org](https://pypi.org) (an account needs two-factor authentication). The PyPI page lists the
   maintainer's PyPI user name publicly.
2. Open **Your account → Publishing → Add a new pending publisher → GitHub** and enter:

   | Field | Value |
   |---|---|
   | PyPI project name | `ts-second-opinion` |
   | Owner | `tammy357` |
   | Repository name | `ts-second-opinion` |
   | Workflow name | `publish.yml` |
   | Environment name | `pypi` |

The MCP Registry needs no setup: it accepts names under `io.github.tammy357/` from this repository's workflow, and checks
that the package's README on PyPI contains `mcp-name: io.github.tammy357/ts-second-opinion` (it does, as a hidden comment
near the top of README.md).

## Each release

1. Set the same version in `pyproject.toml`, `tsso/__init__.py`, both `version` fields of `server.json` and
   `CITATION.cff`, and add a section to `CHANGELOG.md`. The workflow stops if the numbers disagree.
2. Merge to `main`, tag it (`git tag v0.3.0 && git push origin v0.3.0`), then on GitHub open **Releases → Draft a new
   release**, choose the tag and press **Publish release**.
3. Watch the run under **Actions → publish**. When it is green:
   - `pip install ts-second-opinion` installs it, and `uvx ts-second-opinion` starts the MCP server;
   - the server is listed at `https://registry.modelcontextprotocol.io/v0/servers?search=ts-second-opinion`.

If the PyPI step fails with an error about trusted publishing, check that the five values above match exactly. A version
already on PyPI cannot be uploaded again, so fix the problem and release the next version number.
