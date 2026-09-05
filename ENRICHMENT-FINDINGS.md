# Enrichment benchmark: what it measured, and why it cannot answer the question yet

Run 2026-09-04 into 2026-09-05. 252 dockings across two arms. AutoDock Vina
1.2.7, receptor 3A4A prepared, box 30x30x30 at (21.52, -7.70, 23.55),
exhaustiveness 32, seed 42 on both sides for symmetry.

**Headline: no verdict on the ranking claim. The benchmark has a construction
flaw that makes one arm uninterpretable, and the other arm lacks the
resolution to decide.** This is not support for the ranking and it is not a
refutation of it.

---

## The question

The accuracy panel already failed: Spearman +0.09 to +0.16 between predicted
score and published yeast IC50, indistinguishable from chance. That tested
rank-ordering by potency, which docking is known to be poor at.

Enrichment asks the easier and more relevant question: can the protocol put
KNOWN ACTIVES above molecules that merely look like them? That is binary
discrimination, it is what docking is decent at, and it is the ability the
thesis's shortlist claim rests on.

## Results as measured

| Set | actives | ROC AUC | 95% CI | verdict |
|---|---|---|---|---|
| Six flavonoids | 6 | 0.588 | 0.456 to 0.728 | could not tell |
| Triterpenes, neutral COOH | 3 | 0.082 | 0.018 to 0.193 | anti-correlated |
| Triterpenes, deprotonated COO- | 3 | 0.102 | 0.068 to 0.153 | anti-correlated |
| All nine, correct forms | 9 | 0.421 | 0.299 to 0.565 | could not tell |

## ⚠️ The triterpene result is an artefact of the decoy set, not a protocol failure

Decoys were matched on heavy atoms, molecular weight, logP, rotatable bonds,
hydrogen-bond donors and acceptors. They were **not** matched on aromaticity,
and across all 248 scored molecules aromaticity is the strongest predictor of
the Vina score in the whole dataset:

    Spearman, aromatic ring count vs score      -0.544
    Spearman, fraction sp3 vs score             +0.118

    0 aromatic rings   n= 25   median  -8.668
    1-2 aromatic rings n=104   median  -9.149
    3+ aromatic rings  n=119   median -10.000

More aromatic rings score better, by roughly 0.85 kcal/mol per step.

The triterpene actives are saturated pentacyclic cages: **median 0 aromatic
rings, fraction sp3 0.90**. Their property-matched decoys, drawn from ZINC
drug-like space, are quinolines, fluorophenyl amides and pyrazolopyrimidines:
**median 4 aromatic rings, fraction sp3 0.19**.

The actives were asked to out-score molecules from a chemical class the
scoring function systematically prefers, on the one axis nobody matched. They
lost 19 times out of 20. That measures Vina's scaffold preference, not whether
ursolic acid inhibits alpha-glucosidase.

**This arm should not be reported as a result about the protocol.**

## The protonation hypothesis was tested and is dead

The first run took SMILES from PubChem, which returns the NEUTRAL molecule, so
the three triterpene acids were docked as COOH rather than the C-28
carboxylate the 2026-08-13 correction adopted. That was a real error and it
was corrected: 60 fresh decoys, disjoint from the first arm, actives
deprotonated and verified at charge -1 (C30H47O3-, logP 7.1 -> 5.8, HBD 2 -> 1).

It changed almost nothing: **AUC 0.082 -> 0.102**. Both conclusively
anti-correlated. The protonation state was not the explanation; the scaffold
mismatch was.

Worth keeping anyway: the corrected arm is the only measurement of the
deprotonated form against matched decoys that exists.

## The flavonoid arm is the trustworthy one, and it is inconclusive

Flavonoid actives carry 3 aromatic rings against a decoy median of 2, a gap of
-1 that runs slightly in the ACTIVES' favour. So this arm is matched on the
axis that mattered, and mildly flattered rather than sabotaged.

It still returns **AUC 0.588, CI 0.456 to 0.728**, spanning chance. Six actives
is too few. The power simulation run before any docking predicted exactly this:
below roughly 0.6 kcal/mol of separation, this design cannot decide, and a null
effect reads conclusive in only 1 percent of runs.

Note the residual: the actives hold a one-ring aromatic advantage worth about
0.85 kcal/mol and still only reach 0.588. That is not encouraging, but with six
actives it is not evidence either.

## What would fix it

1. **Match decoys on aromaticity and fraction sp3**, not only on size and logP.
   This is the missing criterion and it is the whole finding.
2. For the triterpenes specifically, that requires **terpenoid-like decoy
   space**, which ZINC drug-like tranches do not contain. A natural-product
   library is needed, or the triterpene arm cannot be built honestly.
3. More actives. Nine is few; six per arm is fewer.

## Construction notes that must travel with any use of this data

- **Charge matching was dropped in the deprotonated arm.** An exact match was
  unusable: one candidate at charge -1 in the whole 1.7 million molecule
  library. Defensible for Vina, which has no electrostatic term, and its
  physical consequences (HBD count, logP) were matched. Still a relaxation.
- **Decoys are presumed inactive, not measured inactive.** Any real active
  hiding among them lowers the measured enrichment, so this cannot inflate a
  result.
- 4 ligands of 252 could not be embedded by RDKit and are excluded, each with
  a run record saying so.
- Decoy sets are disjoint between arms, so the two are independent.

## Records

    validation-run/enrichment/                 neutral arm, 186 scored
    validation-run/enrichment-deprotonated/    corrected arm, 62 scored
    validation-inputs/enrichment/decoys.json
    validation-inputs/enrichment/decoys_deprotonated.json
