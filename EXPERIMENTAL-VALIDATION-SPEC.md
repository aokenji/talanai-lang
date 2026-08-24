# Specification for experimental validation

This document is a deliverable of the thesis, not a wish list appended to it.
The title commits to "a specification for experimental validation" as its
final clause, and this is that specification: precise enough that someone with
a glucose meter, a reference standard and a plant sample could run it without
guessing at a single parameter.

It exists because the docking work in this study cannot, on its own, establish
that any compound inhibits anything. It predicts occupancy of a catalytic
site. Whether that translates into measured inhibition is exactly the
question this specification answers how to ask.

---

## 1. What this experiment does and does not decide

**Decides:** whether crude *Ziziphus talanai* leaf extract inhibits yeast
α-glucosidase activity, and by how much relative to acarbose, under a
protocol whose every parameter is stated below.

**Does not decide:** which specific compound in the extract is responsible.
That requires fractionation and is explicitly out of scope; see §7.

**Does not decide:** whether the docking ranking predicts potency. That
question was already tested and answered in the negative (Spearman +0.09 to
+0.16 against published IC50, see `VALIDATION-FINDINGS.md` §4e-4f). This
experiment tests the plant, not the docking.

---

## 2. Why a personal glucose meter, and why that is not a compromise

The published method is Zhang et al. (2021), *Molecules* 26(15):4638,
"Personal Glucose Meter for α-Glucosidase Inhibitor Screening Based on the
Hydrolysis of Maltose" (PMC8348101), applied in that paper to 34 pure
compounds and 18 medicinal plant extracts.

The chemistry: α-glucosidase hydrolyses maltose to glucose. A personal
glucose meter (PGM) reads glucose directly via its internal
glucose-oxidase/ferricyanide electrode. An inhibitor in the reaction well
lowers the glucose yield, and the meter reads lower. No plate reader,
spectrophotometer, or 405 nm anything is required.

**Why the target enzyme is not a mismatch.** PDB 3A4A, the receptor used
throughout this thesis, is itself a *Saccharomyces cerevisiae* isomaltase.
Testing against commercial baker's-yeast α-glucosidase is not a different
system pressed into service; it is the same organism the computational work
already models. The surrogate-for-human caveat that applies to the docking
(§2 of `VALIDATION-FINDINGS.md`) applies identically here, and no more.

---

## 3. Materials, with substitution rules stated up front

| Item | Specification | Budget source |
|---|---|---|
| Enzyme | α-glucosidase, *S. cerevisiae*. Preferred: Sigma-Aldrich G0660 or equivalent commercial preparation with a stated unit activity. **Fallback**: crude preparation from commercial baker's yeast (see §3a) if no commercial enzyme can be sourced or shared. | Institutional reagent order, or a shared aliquot from a partner lab (see §8) |
| Substrate | Maltose, ACS/USP grade | Chemical supplier |
| Positive control | Acarbose. Preferred: analytical reference standard. **Fallback**: crushed Glucobay tablet, with excipient content disclosed and a blank run to characterise excipient interference (see §3b) | Pharmacy (fallback only) |
| Buffer | Sodium phosphate, pH 6.8, or sodium acetate, pH 5.5 (match to the chosen enzyme's stated optimum) | Chemical supplier |
| Glucose meter | Any mainstream consumer PGM with a stated linear range covering 0-30 mmol/L | Retail, or borrowed from a health/nursing department |
| Glucose test strips | Matched to the meter, sufficient for standard curve + full run in triplicate + repeats | Retail |
| Plant material | *Ziziphus talanai* leaves, with a herbarium voucher specimen and accession number (see §9) | Field collection |

### 3a. If no commercial enzyme can be sourced: crude yeast preparation

Baker's yeast (*S. cerevisiae*) is lysed and the soluble fraction used as a
crude α-glucosidase source. This is an accepted low-resource substitution,
not an improvised one, but it carries real limitations that must be stated
in the results, not discovered by a reader:

- Specific activity cannot be reported in defined units (U/mL), only as a
  fixed dilution of a fixed lysate batch.
- The lysate contains more than one glucosidase isoform and unrelated
  enzymatic activity.
- **The validation gate in §5 exists specifically to certify that a crude
  prep is usable despite these limitations.** A crude prep that passes §5 is
  usable; one that does not, is not, regardless of how much material it took
  to make.

### 3b. If no analytical acarbose standard is available

A crushed Glucobay (acarbose 50 mg or 100 mg) tablet may be used, with two
conditions stated in the results: the nominal label strength is not a
verified purity, and a tablet-excipient blank (inactive-ingredient-only,
prepared by dissolving a matched mass fraction with no acarbose) must be run
alongside to characterise any excipient contribution to the glucose signal.

---

## 4. Extract preparation

1. Collect leaves from a voucher-identified specimen (§9).
2. Air-dry away from direct sunlight; do not oven-dry above 40°C, to avoid
   thermal degradation of glycosides.
3. Powder to pass a 1 mm mesh.
4. Macerate in 70% ethanol, 1:10 w/v, 48 hours at room temperature with
   intermittent agitation. 70% ethanol is chosen because it is the higher-
   yield solvent condition in the only published extraction of this species
   (Reyes et al. 2018) and because most of the candidate compound classes in
   this study (flavonoid glycosides, triterpene acids) are ethanol-soluble.
5. Filter (Whatman No. 1 or equivalent), then evaporate the ethanol at
   reduced pressure and low temperature (rotary evaporator, bath ≤ 40°C, or
   equivalent).
6. Reconstitute the dried extract in the assay buffer, or in a minimal volume
   of DMSO followed by buffer dilution (final DMSO ≤ 1% v/v in the assay
   well, and include a DMSO-only blank at that concentration if used).
7. Record the extract yield (dry mass recovered / dry mass of starting
   material, as a percentage) and the final stock concentration in mg dry
   extract per mL.

---

## 5. The validation gate: run before touching the plant extract

This is the step most protocols of this kind skip, and skipping it is why
uncontrolled assays are not believed. **Do not proceed to §6 until this gate
passes.**

1. **Glucose standard curve, in assay buffer, not water.** Prepare 5-6 known
   glucose concentrations spanning the expected assay range (e.g. 0, 1, 2, 4,
   8, 16 mmol/L) in the SAME buffer used for the enzyme reaction. Read each on
   the PGM in triplicate. **Acceptance criterion:** linear regression R² ≥
   0.98 across the working range. Consumer meters are calibrated for whole
   blood; this step is what certifies the meter also behaves linearly in
   phosphate or acetate buffer, and it is not safe to assume.

2. **Enzyme activity check.** Incubate enzyme with maltose alone (no
   inhibitor) and confirm glucose production is (a) approximately linear with
   incubation time over the intended assay window, and (b) approximately
   linear with enzyme dilution over at least a 2-fold range. This confirms
   the reaction is in an initial-rate regime rather than substrate-depleted
   or saturated.

3. **Acarbose recovery.** Run a 5-6 point acarbose dilution series against
   the enzyme and maltose, fit an IC50. **Acceptance criterion:** the fitted
   IC50 falls within the published range for yeast α-glucosidase, 91-841 µM
   across five literature sources (see `VALIDATION-FINDINGS.md` §4e). If the
   system reproduces a known result, every subsequent measurement on this
   system inherits that credibility. If it does not, do not proceed: diagnose
   the enzyme source, buffer, or meter before spending plant material.

4. **Extract-only blank.** Run the plant extract alone, with buffer and
   maltose but no enzyme. *Z. talanai* leaf extract contains its own free
   sugars, and the PGM cannot distinguish glucose released by the enzyme from
   glucose already present in the extract. **This value must be subtracted
   from every extract-plus-enzyme reading in §6.** Skipping this step is the
   single most likely way this experiment produces a false positive.

---

## 6. The extract assay

1. Prepare a 5-6 point dilution series of the reconstituted extract (§4.6)
   spanning at least two orders of magnitude in concentration, informed by
   the extract yield from §4.7.
2. For each concentration: extract + buffer + enzyme, pre-incubate per the
   enzyme's stated optimum (typically 10-15 min at 37°C), add maltose to
   start the reaction, incubate for the fixed window established in §5.2,
   read glucose on the PGM.
3. Subtract the matched-concentration extract-only blank (§5.4) from each
   reading.
4. Run every concentration in triplicate, on at least two independent days
   (biological or preparation replicates, not just technical repeats within
   one run).
5. Express inhibition as: `% inhibition = 100 × (1 − [glucose with
   extract] / [glucose, enzyme-only control])`, using the blank-subtracted
   values.
6. Fit a dose-response curve and report an IC50 (or, if the extract does not
   reach 50% inhibition within a practical concentration range, report the
   maximum inhibition observed and the concentration at which it occurred,
   and say plainly that no IC50 could be determined).

---

## 7. What is explicitly out of scope, and why that is stated rather than hidden

- **Fractionation and compound-level attribution.** This experiment tests
  the crude extract. It cannot say which compound, or which of the ten
  screened phytochemicals, if any, is responsible for an observed effect.
  Bioassay-guided fractionation is the correct next step if inhibition is
  observed, and is a separate, larger undertaking.
- **Kinetic mechanism.** This assay measures endpoint inhibition, not
  competitive/non-competitive/mixed kinetics. Mechanism determination
  requires a Lineweaver-Burk or equivalent analysis at multiple substrate
  concentrations and is not part of this specification.
- **Human enzyme.** This uses yeast α-glucosidase, matching the docking
  target. It does not test human intestinal α-glucosidase (GH31) and cannot
  resolve the GH13-vs-GH31 surrogate question already disclosed in
  `VALIDATION-FINDINGS.md` §2.

---

## 8. Sourcing the enzyme: the realistic path

Mabalacat City College has no BSL-1 facility. This is not a blocker for this
specific experiment: no infectious agent is involved anywhere in this
protocol, only a purified commercial enzyme, so the facility requirement is
ordinary bench chemistry space, not containment.

The realistic obstacle is instead simply having the enzyme in hand. The
recommended path, in order:

1. Request a small aliquot (10 mg is more than sufficient) from a
   pharmacognosy or biochemistry laboratory at a nearby institution with an
   active research pharmacy programme. This is routine inter-institutional
   courtesy for undergraduate research and costs the lending lab
   effectively nothing.
2. If unavailable, use the crude baker's-yeast preparation (§3a), which
   removes the sourcing obstacle entirely at the cost of the stated
   limitations, all of which are gated by §5.3 before any result is trusted.

---

## 9. Voucher specimen

Collect and deposit a herbarium voucher specimen with an accession number at
a recognised herbarium before or concurrent with extraction. State the
specimen's accession number, collection locality, and collection date in any
report of this experiment's results. This is inexpensive, is expected
practice for any study making a species-level botanical claim, and its
absence is a named, specific weakness a reviewer will look for.

---

## 10. What a result would mean, stated for both outcomes in advance

**If the extract shows dose-dependent inhibition with an IC50 in a plausible
range:** this is the first measurement connecting *Ziziphus talanai* itself,
rather than a congener, to the target this thesis studies. It converts the
docking from "a rationale for a class-level inference" into "a mechanistic
hypothesis for an observed effect," and materially strengthens the honest
claim ceiling stated in `VALIDATION-FINDINGS.md` §4e.

**If the extract shows no inhibition, or inhibition inconsistent with the
predicted mechanism:** this is also a reportable result. It would mean either
that the compounds this study screened are not present in this species at a
relevant concentration, that they are present but the extraction method did
not recover them, or that the class-level inference from congeners does not
hold for this species. Any of these is worth stating plainly, and none of
them is a failure of the computational work, which never claimed to predict
this outcome with confidence in the first place.

---

## Summary of what this specification cost to write and what it would cost to run

Writing this: an afternoon, and no laboratory access.

Running it, per the sourcing path in §8 and the budget in the parent
thesis-status document: under PHP 15,000 excluding the enzyme if commercial
enzyme is shared rather than purchased, roughly one to two weeks of bench
time including the mandatory validation gate in §5, and no equipment beyond
a consumer glucose meter and standard glassware.
