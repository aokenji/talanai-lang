# Validation findings, 2026-08-03

Supersedes the interim notes in `CONTROL-RESULTS.md` and
`RING-CONFORMATION-FINDING.md` where they disagree. Two hypotheses stated in
those files were tested and **refuted**; both refutations are recorded below
rather than removed.

---

## 1. The headline

**The published redocking control passes only because it reuses the crystal
ligand's own experimentally refined geometry.** When the ligand is built
independently, exactly as every screened compound was, the same protocol misses
the crystal pose by 5.9 Å and no amount of search effort recovers it.

This is self-docking bias, a known pitfall: a ligand refined into its own
density carries internal geometry already adapted to the bound state, so
redocking it tests far less than it appears to.

| Input | Score | Position RMSD |
|---|---|---|
| Crystal-derived `glc_ligand.pdbqt` | **−5.905** | **0.518 Å** |
| SMILES-built, single conformer | −5.775 | 5.358 Å |
| SMILES-built, 50 conformers, correct chair | −5.465 | 5.919 Å |

RMSD is symmetry-corrected with no superposition (RDKit `CalcRMS`). The
crystal-derived value reproduces the published `validation.json` figure of
0.519 Å, and Talanai's independent standard-library implementation agrees to
0.518 Å. Three routes, one number, so the published control is correct as far
as it goes.

## 2. What was eliminated, and how

Each of these was a live hypothesis. Each was tested and ruled out.

| Hypothesis | Test | Result |
|---|---|---|
| Different molecule or stereochemistry | canonical SMILES and InChI across all files | **Identical.** Ruled out |
| Wrong ring pucker | Cremer-Pople, ring ordered O5 C1 C2 C3 C4 C5 | crystal θ 13.5° and 50-conformer θ 2.1°, both **⁴C₁ chairs**. Ruled out |
| Insufficient search | exhaustiveness 32, 128, 512, three seeds | rank-1 pinned 5.90 to 5.92 Å, scores varying 0.004 kcal/mol. **Ruled out** |
| Different atom typing or torsion tree | PDBQT columns 78-79 and BRANCH records | all files `C=6 HD=5 OA=6`, 6 BRANCH, TORSDOF 6. **Identical.** Ruled out |

What remains is the internal bond lengths and angles, which Vina holds rigid
along with the ring. Crystal geometry comes from experimental refinement
against density; the prepared ligands come from MMFF94 idealisation. The
former is adapted to the pocket, the latter is not.

The decisive observation is the score. Vina reaches −5.905 from the crystal
geometry and never beats −5.465 from the built one, at sixteen times the search
effort. It is not failing to search. It is finding a different global optimum,
because it has been given a subtly different rigid body.

## 3. Two hypotheses I got wrong, recorded

**Maltose from 3AJ7.** I specified this as the control gate. 3AJ7 contains no
sugar at all. A survey of all five PDB entries of this enzyme found no intact
maltose anywhere; 3AXI is even *titled* as a maltose complex while containing
only glucose. Corrected to isomaltose from 3AXH, which is also the better
control since 3A4A is an isomaltase and isomaltose is its cognate substrate.

**Ring pucker as root cause.** I reported this as the root cause before testing
it. Fixing the pucker did not fix the redock. The single-conformer defect is
real and worth repairing, but it is not why the control fails.

Both were caught by checking rather than by reasoning.

## 4. A real defect, independent of the above

`_prep_assets_local.py:41` embeds **one** conformer with
`AllChem.EmbedMolecule` and minimises it locally. MMFF cannot cross a ring-flip
barrier, Meeko writes rings rigid, and Vina cannot alter ring geometry.
Whatever pucker that single attempt produced is frozen for the whole study.

Measured on glucose: the energy spread across 50 embeddings was **20.78
kcal/mol**, and the single shot the pipeline took was roughly 21 kcal/mol above
the best available conformer, in a twist-boat rather than a chair.

### The rigorous audit

A first pass classified rings by torsion sign pattern, which flags any ring
that is not a clean chair, including rings that **cannot** be chairs because
they contain an sp2 atom. `audit_rings_rigorous.py` replaces it: it determines
saturation from bond orders, computes Cremer-Pople Q and theta, and judges each
ring against what it is capable of adopting. Only a **saturated** ring that is
not a chair is a defect.

| Compound | Saturated rings | Verdict |
|---|---|---|
| **Ursolic acid** | 3 of 5 at theta 83 to 97° | **defective**, boat family |
| **Isovitexin** | 1 of 1 at theta 91.5° | **defective**, the sugar ring |
| **Acarbose (reference)** | 1 at theta 86.8° | **defective**, a sugar ring |
| Spinosin (prepared today) | 1 at theta 82.5° | defective |
| Isomaltose (prepared today) | 1 at theta 89.9° | defective |
| Rutin | 2 of 2, theta 4.9 and 177.9° | clean chairs |
| Betulinic acid | 4 of 4, theta 6 to 180° | clean chairs |
| Vitexin | 1 of 1, theta 177.8° | clean chair |
| Oleanolic acid | 3 of 3, theta 6 to 180° | clean chairs |
| Quercetin, kaempferol, luteolin | none | aromatic only, cannot be affected |

**Both of the first audit's uncertain calls are now resolved, in opposite
directions.** Oleanolic acid is **clean**: its flagged ring is the legitimate
Δ12-13 half-chair, correctly identified here as sp2-containing. Ursolic acid is
**genuinely broken**, in three separate saturated rings, which the Δ12 double
bond does not explain.

The contrast is itself evidence. Oleanolic and ursolic acid are near-identical
pentacyclic scaffolds. One came out entirely correct and the other came out
with three rings in boat conformations. That is exactly what a single-embedding
lottery produces, and it could not happen if the pipeline searched conformers.

**Scope: two of the ten screened compounds, plus the reference.** Ursolic acid
(−8.105, ranked third) and isovitexin (−8.024, ranked fourth) were docked as
distorted molecules. So was acarbose, which is the yardstick every comparison
in the study is measured against.

**The fix is one line:**

```python
cids = AllChem.EmbedMultipleConfs(mol, numConfs=50, params=params)
AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=400)   # keep the lowest
```

## 4b. Ring-aware preparation, and the corrected numbers (2026-08-04)

The one-line fix (embed 50, keep lowest energy) repaired ursolic acid and
isovitexin but made **acarbose worse**, 1 defective ring becoming 2. Diagnosis:
for a large flexible molecule, gas-phase MMFF energy is dominated by folding
and internal hydrogen bonds, so a tightly folded conformer with a strained ring
outscores an extended one with perfect chairs. Energy selects against ring
quality precisely where ring quality matters.

**The corrected rule: among conformers whose saturated six-rings are all
chairs, keep the lowest energy. Ring correctness filters, energy only breaks
ties.** Vina can relax a slightly strained pose; it can never fix a ring.

| Ligand | All-chair conformers | Energy penalty of the chosen one |
|---|---|---|
| **Acarbose** | **14 of 300** (4.7%) | **+8.45 kcal/mol above global minimum** |
| Ursolic acid | 11 of 150 | none, it was the global minimum |
| Isovitexin | 57 of 150 | none, it was the global minimum |

Acarbose is the only molecule where the two criteria disagree, and at a 4.7%
hit rate a 50-conformer search would expect barely two candidates. That is
exactly why the first attempt failed.

### Arm A, thesis configuration, only the ligand preparation differs

| Compound | Published | Old prep | Ring-aware | Δ vs old |
|---|---|---|---|---|
| Ursolic acid | −8.105 | −9.097 | −7.711 | +1.386 |
| Isovitexin | −8.024 | −7.883 | −8.436 | −0.553 |
| **Acarbose** | **−6.660** | −5.942 | **−6.845** | **−0.903** |

**The acarbose result is the notable one.** The regenerated file that failed to
reproduce the published value (−5.942 against −6.660) now returns **−6.845**,
within 0.185 of the published number. That is consistent with the original
thesis acarbose input having had sound ring geometry, and with the earlier
reproduction failure being an artefact of regenerating that file with
single-conformer preparation. It does not prove it, since the thesis input file
no longer exists, but it is the most likely explanation and it is reassuring
for the published value.

### Arm B, Meeko receptor, exhaustiveness 32, three seeds, ring-aware

| Compound | Mean | Range | Spread |
|---|---|---|---|
| Isovitexin | −9.356 | −9.366 to −9.337 | 0.029 |
| Ursolic acid | −8.717 | −8.735 to −8.700 | 0.035 |
| Acarbose | −8.706 | −8.710 to −8.700 | 0.010 |

**Read these by range, not by mean.**

- **Isovitexin beats acarbose.** A 0.650 kcal/mol gap with non-overlapping
  ranges. Real at this sampling.
- **Ursolic acid is indistinguishable from acarbose.** The means differ by
  0.011 kcal/mol and the ranges overlap: ursolic runs to −8.700 and acarbose
  runs to −8.700. There is no separation to report.

The script initially printed "ursolic acid beats acarbose" from a bare mean
comparison. That was wrong by this project's own standard, the same error
flagged earlier in the kaempferol case, and the verdict logic has been changed
to require range separation before it will call a winner.

Note also that correcting acarbose's rings made it **weaker** on the Meeko
receptor, −8.992 to −8.706, which moved ursolic acid from "loses" to "tied".
Ring geometry moves the reference compound enough to reorder the field.

## 4c. The box test, and the settled conclusion (2026-08-04)

The last untested variable was box size. Every independent-ligand failure had
used the 30 Å screening box, while the only configuration that ever passed used
an 18 Å focused box **and** crystal coordinates. Those two had never been
separated. Exhaustiveness 32, seed 42, Meeko receptor throughout.

| Ligand | Box | Score | rank-1 | best of 9 |
|---|---|---|---|---|
| glucose, crystal-derived | 18 Å | −5.885 | **0.513 PASS** | 0.513 |
| glucose, crystal-derived | 24 Å | −5.680 | **4.431 fail** | 3.842 |
| glucose, crystal-derived | 30 Å | −5.905 | **0.518 PASS** | 0.518 |
| glucose, SMILES, chair | 18 Å | −5.465 | 5.918 | 3.891 |
| glucose, SMILES, chair | 24 Å | −5.466 | 5.918 | 4.378 |
| glucose, SMILES, chair | 30 Å | −5.465 | 5.919 | 4.353 |
| isomaltose, SMILES, chairs | 18 Å | −7.199 | 6.863 | **1.397** |
| isomaltose, SMILES, chairs | 24 Å | −7.119 | 9.782 | 6.877 |
| isomaltose, SMILES, chairs | 30 Å | −7.111 | 9.773 | 6.693 |

### Three things this settles

**1. The box matters, but fixing it is not sufficient.** In the 18 Å box the
correct isomaltose pose enters the ensemble at **1.397 Å**, comfortably inside
the bar. At 24 and 30 Å it is not sampled at all, best-of-nine only reaching
6.7 to 6.9 Å. So a focused box genuinely helps the search. But even when the
correct pose is found, it is **not ranked first**: rank-1 is 6.863 Å. The
search located the answer and the scoring function preferred something else.

That is the signature of a **scoring limitation**, and it is the distinction
the orientation experiment was designed to detect. No amount of extra search
fixes it.

**2. Glucose never recovers from an independent build, at any box size.**
5.918, 5.918, 5.919 across three box sizes, with scores agreeing to 0.001
kcal/mol. Converged, reproducible, and wrong. Its correct pose is not in the
ensemble at any configuration tested.

**3. The crystal ligand is not reliable either, and this corrects an earlier
framing of mine.** It passes at 18 Å and 30 Å but **fails at 24 Å** (4.431 Å),
same ligand, same receptor, same settings, one seed each. At exhaustiveness 32
with a single seed, even the crystal ligand's success is partly luck. The
self-docking-bias finding stands, since independently built ligands never
recover while crystal ones sometimes do, but "the crystal ligand reliably
recovers" was too strong and the 24 Å run refutes it.

It also shows why single-seed results should never be trusted. I nearly drew a
conclusion from a configuration that a second seed would have exposed.

### The settled conclusion

**This protocol can rank compounds. It cannot place them.**

Affinity scores reflect pocket occupancy and burial, and on that basis the
screen does what a screen is for: it produces a shortlist. Every control pose,
including the failing ones, sits in the catalytic site and hydrogen-bonds the
full triad.

But pose-level claims are not supportable. The protocol cannot reliably
reproduce a known pose from an independently generated ligand, the correct pose
is often not sampled in the screening box, and when it is sampled it is not
reliably ranked first.

**Consequences, and these are now conclusions rather than concerns:**

1. `/interactions` is correctly gated. Per-compound residue contacts derived
   from a rank-1 pose are not supportable by this protocol.
2. The redocking control must be reported as what it is: a self-docking test
   that passes with crystal coordinates and does not validate the screening
   pipeline.
3. If pose-level work is wanted later, it needs a focused box, ensemble
   reporting rather than rank-1, at least three seeds, and probably rescoring
   with a second function. That is a project, not a parameter change.
4. The affinity ranking survives, with the caveats already documented: size
   bias, the acarbose comparison, and in-silico scope.

## 4d. The full corrected screen (2026-08-04)

All ten screened compounds plus the reference, re-prepared ring-aware and
re-docked on the Meeko receptor at exhaustiveness 32 across three seeds. Run
records in `validation-run/full-screen/`.

**Preparation.** Every ligand's saturated rings verified as chairs. The
all-chair hit rate varied enormously and shows how much of a lottery
single-conformer embedding was: Rutin 10 of 200, Acarbose 14 of 300, Ursolic
11 of 150, against Betulinic 77 of 150 and Vitexin 61 of 150. Luteolin,
quercetin and kaempferol have no saturated ring and cannot be affected.

| Compound | Corrected mean | Spread | Published | Δ |
|---|---|---|---|---|
| Rutin | **−10.317** | 0.020 | −8.857 | −1.460 |
| Betulinic acid | −9.720 | 0.026 | −8.290 | −1.430 |
| Vitexin | −9.490 | 0.043 | −7.469 | −2.021 |
| Isovitexin | −9.356 | 0.029 | −8.024 | −1.332 |
| Spinosin | −9.257 | 0.016 | −7.859 | −1.398 |
| Quercetin | −8.818 | 0.020 | −7.503 | −1.315 |
| Ursolic acid | −8.717 | 0.035 | −8.105 | −0.612 |
| **Acarbose** | **−8.706** | 0.010 | −6.660 | **−2.046** |
| Kaempferol | −8.510 | 0.042 | −7.479 | −1.031 |
| Luteolin | −8.484 | 0.028 | −7.733 | −0.751 |
| Oleanolic acid | −7.919 | 0.022 | −6.922 | −0.997 |

Spreads of 0.010 to 0.043 kcal/mol confirm that exhaustiveness 32 with three
seeds converges properly, unlike the exhaustiveness 8 single-seed original.

### The result

**Six of ten beat acarbose. One is indistinguishable. Three lose to it.**

| Verdict | Compounds |
|---|---|
| Beats acarbose | Rutin, Betulinic acid, Vitexin, Isovitexin, Spinosin, Quercetin |
| Indistinguishable | Ursolic acid, at 0.011 kcal/mol with overlapping ranges |
| Loses to acarbose | Kaempferol, Luteolin, Oleanolic acid |

The published claim is ten of ten. Verdicts here require the three-seed ranges
to separate before a winner is named; a bare mean comparison would have said
seven of ten and would have been wrong about ursolic acid.

### Why the field reordered

**Every compound scored more strongly on the corrected protocol**, by 0.6 to
2.0 kcal/mol, which is expected: the raw receptor's hydrogen-bond term was
measured at exactly zero, so restoring it helps everything.

But **acarbose gained the most of any compound, −2.046**. It is the largest and
most hydrogen-bond-capable molecule in the set, so it had the most to gain from
a receptor that can form hydrogen bonds at all. That single fact is what moved
four compounds from "beats the drug" to "ties or loses".

Vitexin is the other mover, from ninth place to third.

### What did not change

**Rutin remains the strongest binder by a clear margin**, −10.317 against
acarbose's −8.706, a 1.6 kcal/mol gap with non-overlapping ranges. The
headline result of the study survives the correction comfortably.

On ligand efficiency, acarbose is still the worst in the set at 0.198 per heavy
atom, behind every candidate including the three that lose on raw score.
Kaempferol (0.405), Luteolin (0.404) and Quercetin (0.401) lead. So the
size-corrected picture is more favourable than the raw one, and both should be
reported.

### Three cautions, stated plainly

1. **This is not a validated protocol.** The redocking gate does not pass. This
   is a ranking, and section 4c is why no pose or interaction claim may be
   drawn from it. Every run record carries `"validated": false`.
2. **Three variables changed at once** relative to the published screen: the
   receptor preparation, the ligand preparation and the sampling. So the Δ
   column is not attributable to any single cause. Arm A in section 4b is the
   only clean single-variable comparison available.
3. **Ring-aware preparation changed every ligand**, including those whose
   original rings were already chairs, because conformer selection changed for
   all of them. This is a like-for-like set internally, not a patch applied to
   some compounds and not others.

## 4e. External validation: the answer to "is this credible" (2026-08-05)

Redocking validates pose prediction. Consensus validates robustness to the
scoring function. Neither says whether the ranking is **right**. Only measured
data does. Three tests were run.

### Test 1: consensus between scoring functions

Vinardo (Quiroga and Villarreal, PLoS ONE 2016) ships with Vina 1.2.x, so this
needed no extra software. Run twice: rescoring the same Vina poses, and letting
Vinardo drive its own independent search.

**Both give Spearman 0.336.** The agreement is weak, and it is not an artefact
of rescoring poses optimised for another function, because the independent
re-dock reproduces it exactly.

| Compound | Vina rank | Vinardo rank | Moved |
|---|---|---|---|
| **Rutin** | 1 | **1** | — |
| Betulinic acid | 2 | 9 | **+7** |
| Vitexin | 3 | 2 | −1 |
| Isovitexin | 4 | 8 | +4 |
| Spinosin | 5 | 6 | +1 |
| Quercetin | 6 | 4 | −2 |
| Ursolic acid | 7 | 10 | +3 |
| Acarbose | 8 | 7 | −1 |
| Kaempferol | 9 | 5 | −4 |
| Luteolin | 10 | 3 | **−7** |
| **Oleanolic acid** | 11 | **11** | — |

**The extremes are robust; the middle is not.** Rutin is first and oleanolic
acid last under both functions. Between them the ordering is largely a property
of which function you use. Five of eleven compounds change sides on "beats
acarbose".

### Test 2: published IC50 data

Yeast *S. cerevisiae* α-glucosidase only, the assay matching 3A4A. Ten of
eleven compounds have usable data. **Spinosin has none**: its literature covers
sedative and neuroprotective activity, not glycosidase inhibition.

Full sourcing, quotes and citations in `validation-run/ic50/LITERATURE-IC50.md`.

**Acarbose is a weak inhibitor of the yeast enzyme.** Published values span
91 to 841 µM across five sources, against 5 to 17 µM for the triterpenes and
17 to 117 µM for quercetin, measured in the same assays. Primary confirmation,
Oki, Matsui and Osajima, *J Agric Food Chem* 1999, PMID 10563931: "acarbose...
strongly inhibited mammalian AGHs, whereas no or less inhibition was observed
in yeast AGH."

**So "beats acarbose" is a low bar in this system**, and it is one most of
these compounds clear in real wet-lab yeast assays independently of anything
docking says. The comparison anchoring the entire study is against a drug that
barely works on this enzyme.

Note also the within-compound spread: rutin ranges 13 to 196 µM and ursolic
acid 5 to 17 µM across labs, and the raw literature spread before filtering ran
to 60 and 90 fold. Four values were discarded as uncitable, including a rutin
figure that would have been a 350-fold outlier.

### Test 3: does the ranking track measured potency?

Spearman between predicted score and log IC50, computed three ways using the
lowest, geometric-mean and highest published value per compound, so the answer
cannot depend on which paper is cited.

**Sign convention**, which is easy to invert and was inverted in the first
version of this analysis: a potent compound has a **very negative score** and a
**low IC50**. Both low together, so a correct prediction gives a **POSITIVE**
Spearman.

| Scoring function | Range | Verdict |
|---|---|---|
| **Vina** | +0.091 to +0.164 | **No predictive value.** Indistinguishable from chance |
| **Vinardo** | −0.345 to −0.442 | **Anti-predictive.** Systematically prefers the weaker compounds |

The result swings only 0.073 across the three literature choices, so it is not
an artefact of citation selection.

Concretely: ursolic acid is the most potent compound measured, 10.1 µM, and
Vina ranks it seventh. Oleanolic acid is second most potent and receives the
**worst score in the set**. Rutin has the best predicted score and middling
measured potency.

### The settled answer to "is this credible"

**Internally, yes.** Converged sampling, verified ligand geometry, checksummed
inputs, a stored run record for every job. That is better practice than most
published docking work.

**Externally, no demonstrated predictive validity.** The ranking does not track
measured potency, two scoring functions disagree at ρ = 0.336, and the
reference compound is a poor choice for this enzyme.

**What can honestly be claimed:**

> These compounds occupy the catalytic site of a yeast α-glucosidase model, and
> the published literature independently reports several of them as inhibitors
> of that enzyme. The docking did not establish that; it is consistent with it.

**What cannot be claimed:** that the docking ranking predicts potency, that a
compound scoring better will inhibit better, or that beating acarbose in this
system indicates clinical promise.

## 4f. Is rho = 0.12 my failure or the method's ceiling? (2026-08-05)

Section 4e reports the number but not its meaning. A panellist will ask whether
that correlation is bad because the protocol is bad. Three published reference
points answer it, and they do not all point the same way.

**Reference 1. Vina can rank, under favourable conditions.** CASF-2016
(Su et al. 2019, *J Chem Inf Model* 59:895) evaluates scoring functions on 57
targets with 5 ligands each. AutoDock Vina scores **ranking power Spearman
rho = 0.528** and **scoring power Pearson R = 0.604**. So the function is not
inherently incapable of ordering ligands, and I cannot claim it is.

**Reference 2. Under conditions closer to ours, it does not.** A 2026 benchmark
in *Briefings in Bioinformatics* (bbag028) tested potency ranking of **human
acetylcholinesterase** inhibitors: 412 curated compounds, five receptor
structures, Glide rather than Vina. **Single-structure non-covalent docking
gave rho = 0.07 to 0.18.** Our +0.09 to +0.16 sits inside that interval.

⚠️ **Caveat that must be stated whenever this is cited.** That library is
carbamates and organophosphates, which inhibit AChE **covalently**. Docking
them non-covalently leaves the actual mechanism unmodelled, so 0.07 to 0.18 is
a pessimistic bound, not a like-for-like comparison. The authors' point is that
covalent docking lifts it to 0.54 and ML consensus to 0.70. Do not present this
as proof that docking cannot rank. It is corroborating context.

**What the gap between 0.53 and 0.12 actually is.** CASF's favourable
conditions are all absent here, and each absence is a named, defensible reason:

| CASF-2016 ranking power test | This study |
|---|---|
| 5 **congeneric** ligands per target | 10 compounds across 4 chemotypes: flavonol glycoside, flavone C-glycoside, aglycone, pentacyclic triterpene |
| Affinity range selected to span decades | Measured potencies cluster, and overlap their own error bars |
| Affinities from curated single sources | IC50 pooled across different labs, enzyme lots and assay variants |
| Crystal structure of each complex | One receptor, no complex structures for these ligands |

Across chemotypes rather than within a series, Vina's score is dominated by
size and lipophilicity. Rutin, the largest molecule in the set, scores best.
That is the expected behaviour of an empirical scoring function, not a bug.

**The reference data is itself a limitation.** Quercetin's published
alpha-glucosidase and AChE IC50 values span roughly **50-fold** across papers.
A correlation computed against numbers that disagree with each other by that
much has a low ceiling no matter how good the docking is.

**Refinement of the 4e wording, and this matters.** "No demonstrated predictive
validity" is correct and stands: the test was run and it did not demonstrate
validity. But **this test cannot separate a method limitation from reference
data noise**, and 4e should not be read as proving docking cannot predict
potency in general. The honest statement is the narrower one.

**The concrete fix, if it is ever worth doing.** Correlate against a
**single-protocol reference set**: several compounds measured in one lab, one
enzyme source, one assay. Then the test measures the method instead of
measuring the literature's disagreement with itself.

**Why this is good for the thesis, not bad.** Almost no student docking study
measures its own rank correlation at all. Reporting rho = +0.12 with the
benchmark context above is a stronger position than asserting predictive power
and being asked for evidence. The finding to defend is: *this protocol was
tested for potency prediction, it did not demonstrate it, and that outcome is
consistent with published benchmarks for single-structure docking of a
chemically heterogeneous set against pooled literature activity data.*

## 4g. Ligand stereochemistry: three compounds were not the molecules they were named (2026-08-11)

An adversarial examiner panel, assembled to stress-test the thesis before
defence, inspected the actual ligand SMILES used in the corrected screen
rather than taking the compound names on trust. One examiner reported that
Ursolic Acid's SMILES carried no stereo descriptors at all. That claim was
verified independently before anything was acted on: the real molecule
(PubChem CID 64945) has ten stereocenters, and RDKit confirmed the SMILES
used in `run_full_screen.py` specified zero of them. A systematic sweep of
all eleven compounds in the screen followed, checking every ligand rather than
spot-checking the one flagged.

| Compound | Stereocenters specified | Status |
|---|---|---|
| Ursolic Acid | 0 of 10 | Defective |
| Vitexin | 0 of 5 | Defective |
| Isovitexin | 0 of 5 | Defective |
| Rutin, Betulinic Acid, Spinosin, Oleanolic Acid, Acarbose | fully specified | OK |
| Luteolin, Quercetin, Kaempferol | no stereocenters exist | OK, nothing to specify |

Three of eleven, and the missing stereochemistry is not a minor omission.
An unspecified stereocenter is not a smaller version of the correct
structure: RDKit's conformer generator picks an essentially arbitrary 3D
shape at every undefined centre. What was docked under the name "Ursolic
Acid" was, with near certainty, not the natural product. For Vitexin and
Isovitexin the missing centres were on the glucosyl ring, the very feature
that makes them C-glycosides rather than the flavone aglycone.

This is standard ligand-preparation practice, not a judgement call. Assigning
protonation state, tautomer and stereochemistry before conformer generation
is step one of any serious docking protocol; a docking program has no way to
know it was handed the wrong molecule.

**The fix.** Fully stereo-specified isomeric SMILES were re-fetched from
PubChem, using the same CID already recorded in the compound catalogue rather
than a fresh guess. One fetch initially returned the wrong molecule for
Vitexin (a different CID, caught immediately because the formula did not
match C21H20O10) - corrected before proceeding. All three replacement
structures were verified both formula-matched and fully stereo-specified with
RDKit before docking. Re-run under the identical 2026-08-04 protocol: same
receptor, same 30 Å box, exhaustiveness 32, seeds 42/43/44. Exactly one
variable changed relative to that protocol - the ligand structure - which
makes this a cleaner correction than 4d, where three variables moved at once.

**The result.**

| Compound | Prior (defective structure) | Corrected | Verdict |
|---|---|---|---|
| Ursolic Acid | −8.717, indistinguishable from Acarbose | −8.651 | loses to Acarbose |
| Vitexin | −9.490, beats | −9.637 | beats, more strongly |
| Isovitexin | −9.356, beats | −8.763 | beats, but by only ~0.03 kcal/mol edge-to-edge |

The 0.011 kcal/mol margin behind the Ursolic Acid tie call in §4d was itself
an artefact of the broken structure. With the real molecule, the separation
from Acarbose is unambiguous and it loses.

Because Ursolic Acid was a **tie**, not a **beat**, this correction moves the
tie and lose counts, not the beat count. The six compounds beating Acarbose
are the same six as before. **New headline: 6 of 10 beat Acarbose
(unchanged), 0 are indistinguishable (down from 1), 4 lose (up from 3).**

**Ranking after correction**, unaffected compounds unchanged:

| Rank | Compound | ΔG (kcal/mol) | Verdict |
|---|---|---|---|
| 1 | Rutin | −10.317 | beats |
| 2 | Betulinic Acid | −9.720 | beats |
| 3 | Vitexin | −9.637 | beats |
| 4 | Spinosin | −9.257 | beats |
| 5 | Quercetin | −8.818 | beats |
| 6 | Isovitexin | −8.763 | beats, fragile margin |
| ref | Acarbose | −8.706 | — |
| 7 | Ursolic Acid | −8.651 | loses |
| 8 | Kaempferol | −8.510 | loses |
| 9 | Luteolin | −8.484 | loses |
| 10 | Oleanolic Acid | −7.919 | loses |

**What this is worth saying plainly.** This is the second time a protocol
audit has changed the headline, and the second time the correction made the
positive claim more conservative rather than less: 10/10 → 6/1/3 in §4d, then
6/1/3 → 6/0/4 here. Neither audit increased the count of compounds beating
the reference. That trajectory is itself evidence the process is finding real
problems rather than reasoning toward a preferred answer, and it is worth
saying so explicitly rather than leaving it implicit.

**What remains open.** The examiner panel separately flagged that the
pipeline has no protonation-state assignment step, which matters specifically
for Betulinic, Oleanolic and Ursolic acid: each carries a carboxylic acid
(pKa approximately 4.5–5) that should likely be deprotonated at
physiological or assay pH. Not investigated yet. Deliberately not folded into
this correction, to preserve the single-variable discipline this fix was
built on.

Record: `talanai-lang/run_stereo_fix.py`, run files and `summary.json` at
`talanai-lang/validation-run/stereo-fix/`. Applied to the live site and
pushed to `main` at `65ee938`.

## 5. What this does and does not touch

**Does not touch.** That the compounds occupy the catalytic site. Every control
pose, including the failing ones, hydrogen-bonds the full triad: Asp215 2.77 Å,
Glu277 2.87 Å, Asp352 2.75 Å for crystal glucose, and comparable distances for
the others. The box is right, the site is right, nothing docks into a decoy
pocket.

**Does touch.**

1. **The validation status of the protocol.** The redocking control does not
   validate the pipeline that produced the results. It validates redocking a
   crystal-coordinate file.
2. **Pose-level claims.** The PLIP per-compound interaction data describes
   orientations this protocol cannot pin down. For glucose, two orientations
   5.4 Å apart score within 0.13 kcal/mol of each other.
3. **Acarbose specifically**, which was docked as a strained conformer and
   anchors every comparison.

## 6. What has not been done

No compound was re-prepared, re-docked, or re-ranked. No published number was
changed. The validation screen never ran, because its gate did not pass and the
threshold was not relaxed to get through it.

Under the delegation charter this is a section-3 item: it changes what the
study can conclude, so the decision belongs to the author.

## 7. Ordered next steps

1. **Settle acarbose and ursolic acid**, separating legitimate sp2-containing
   half-chairs from genuinely distorted saturated rings. Cheap, no docking, and
   it decides how far the preparation defect reaches.
2. **Fix the preparation** to embed many conformers and keep the lowest.
3. **Re-prepare and re-dock** whatever is genuinely affected, reporting deltas
   beside the published values rather than replacing them.
4. **Decide what the redocking control should be.** Options: report it honestly
   as self-docking-biased and drop the pose-level claims; or adopt a
   cross-docking control, where the ligand geometry cannot be adapted to the
   receptor it is docked into, which is what control 3 was reaching for.
5. **Restrict interaction claims to the pose ensemble** rather than the rank-1
   pose, unless a validated protocol says otherwise.
