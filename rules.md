# Rules for a believable transition state

Plain-English rules. `python scripts/author_all.py` asks Claude to turn each one into a Python check,
then keeps a check only if it accepts at least 99% of 9,540 reference DFT transition states.

- chirality: Stereocentres that are not part of the reaction keep the handedness they have in the reactant.
- product_stereo: Atoms that are stereocentres in the product already have the product's handedness at the transition state.
- clash: Atoms that are never bonded to each other, in the reactant or the product, are not squashed unrealistically close together.
- spectator_bonds: Bonds that survive the reaction keep a length close to their lengths in the reactant and product.
- reacting_bonds: Bonds that form or break have lengths in between their reactant and product values.
- planarity: Double bonds and aromatic rings that do not take part in the reaction stay flat.
