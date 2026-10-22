The rule passed the gate on attempt 2 and is saved as `h_collinear`. It never fires on the GoFlow guesses, so it adds nothing to triage.

**Final code** (saved to `tsso/checks_user/h_collinear.py`, `use_in_triage=False`):

```python
def check(rxn, xyz):
    """A transferred H is roughly collinear with donor and acceptor, unless a small ring TS forces it to bend."""
    nb_r = rxn.neighbours('reactant')
    nb_p = rxn.neighbours('product')
    checked = 0
    for h in range(rxn.n):
        if rxn.symbols[h] != 'H':
            continue
        heavy_r = [a for a in nb_r.get(h, []) if rxn.symbols[a] != 'H']
        heavy_p = [a for a in nb_p.get(h, []) if rxn.symbols[a] != 'H']
        if len(heavy_r) != 1 or len(heavy_p) != 1 or heavy_r[0] == heavy_p[0]:
            continue
        d, a = heavy_r[0], heavy_p[0]
        # 1,2- and 1,3-shifts (donor and acceptor bonded or sharing a neighbour) go through 3-/4-membered rings and must bend
        close = (a in nb_r.get(d, [])) or (a in nb_p.get(d, [])) \
            or (set(nb_r.get(d, [])) & set(nb_r.get(a, []))) or (set(nb_p.get(d, [])) & set(nb_p.get(a, [])))
        if close:
            continue
        v1 = xyz[d] - xyz[h]
        v2 = xyz[a] - xyz[h]
        n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        checked += 1
        ang = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2) / (n1 * n2))))))
        if ang < 90.0:
            return False, (f"H{h + 1} moves from {rxn.symbols[d]}{d + 1} to {rxn.symbols[a]}{a + 1} "
                           f"but the angle {rxn.symbols[d]}{d + 1}-H{h + 1}-{rxn.symbols[a]}{a + 1} "
                           f"is only {ang:.0f} degrees at the TS (expected >= 90)")
    if checked == 0:
        return True, "rule not applicable (no transferred hydrogen between distant heavy atoms)"
    return True, f"{checked} transferred hydrogen(s) lie roughly between donor and acceptor"
```

**Attempts:** 2 of 4.
- **Attempt 1:** I required a donor–H–acceptor angle of at least 110° for every transferred H. It wrongly rejected 36% of real TSs, which is far over the 1% allowed. Most of the rejected examples were 1,2- and 1,3-H shifts, which pass through 3- or 4-membered rings and are bent by geometry (angles of 55–105°).
- **Attempt 2:** I skipped hydrogens whose donor and acceptor are bonded or share a neighbour, and relaxed the threshold to 90°. This was accepted.

**Gate numbers:**

| Measure | Result |
|---|---|
| Real TSs accepted | 99.91% (needs ≥99%) |
| Mirrored TSs caught | 0.09% |
| Noisy TSs caught | 1.35% (just above the 1% minimum) |
| Crash rate | 0% |
| `share_of_goflow_guesses_it_flags` | 0.0 |

**Why it adds little:**
- **No GoFlow flags:** none of the 8 guesses per reaction trip it, so it would never change a triage decision.
- **Barely sensitive:** it only just clears the 1% sensitivity minimum on noisy geometries and does almost nothing on mirrored ones.
- **Already covered:** `h_transfer` is an existing rule that does nearly the same thing, and it isn't used in triage either. It passes 99.2% of real TSs, catches 7% of noisy TSs and 1% of mirrored ones. Mine is stricter on the real data but weaker at catching bad geometries, so I'd use `h_transfer` as the base if you want an H-transfer check.
- **Left at default:** I left `use_in_triage` at False, so triage is unchanged.
