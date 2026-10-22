"""Render the hero figure (RDB7 3509, GoFlow seed-4 guess) from results/hero_mirror.html.

The HTML embeds the reference TS and the guess (aligned by proper rotation) with the vendored 3Dmol.js.
Screenshot it with any browser, or: python scripts/make_hero.py  (needs: pip install playwright; playwright install chromium)
"""
import asyncio
from pathlib import Path
from playwright.async_api import async_playwright

ROOT = Path(__file__).resolve().parent.parent


async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await b.new_page(viewport={'width': 1200, 'height': 640}, device_scale_factor=2)
        await pg.goto((ROOT / 'results' / 'hero_mirror.html').as_uri()); await pg.wait_for_timeout(3000)
        await (await pg.query_selector('.wrap')).screenshot(path=str(ROOT / 'results' / 'fig_hero_mirror.png'))
        await b.close()

asyncio.run(main())
