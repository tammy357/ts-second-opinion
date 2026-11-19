"""Helpers that connect TS Second Opinion to AI agents: install the agent skill and print MCP client settings."""
import json, shutil
from pathlib import Path

SKILL = Path(__file__).resolve().parent / 'skill' / 'SKILL.md'
GIT = 'git+https://github.com/tammy357/ts-second-opinion'


def install_skill(where='user'):
    """Copy the skill to ~/.claude/skills (where='user'), ./.claude/skills (where='project') or a given folder."""
    base = {'user': Path.home() / '.claude' / 'skills', 'project': Path.cwd() / '.claude' / 'skills'}.get(where, Path(where))
    target = base / 'ts-second-opinion' / 'SKILL.md'
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(SKILL, target)
    return target


def mcp_settings(client='claude-desktop'):
    """Settings that start the MCP server with uvx, so nothing needs installing by hand."""
    args = ['--from', GIT, 'tsso', 'mcp']
    if client == 'claude-code':
        return 'claude mcp add ts-second-opinion -- uvx ' + ' '.join(args)
    if client in ('claude-desktop', 'cursor'):
        return json.dumps({'mcpServers': {'ts-second-opinion': {'command': 'uvx', 'args': args}}}, indent=2)
    raise ValueError("client must be 'claude-desktop', 'cursor' or 'claude-code'")
