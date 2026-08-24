# Literature IC50 Values for Alpha-Glucosidase Inhibition (Yeast / Saccharomyces cerevisiae)

Compiled 2026-08-05 to validate a molecular docking ranking against published experimental data. Target assay: **yeast (Saccharomyces cerevisiae) alpha-glucosidase**, the enzyme matching PDB 3A4A, almost always assayed against **pNPG** (p-nitrophenyl-alpha-D-glucopyranoside) substrate.

Every value below was retrieved from an identified, citable source. Where a number could not be traced to a specific paper, it was discarded rather than reported (see "Discarded / uncitable values" at the end). No value was estimated, rounded from memory, or inferred.

Full machine-readable version: `D:\THESIS_VSC\talanai-lang\validation-run\ic50\ic50.json`

---

## Summary table

| Compound | Yeast IC50 range found | Independent sources | Enzyme source confirmed as yeast? | Overall confidence |
|---|---|---|---|---|
| Acarbose | 91 - 841.3 uM (quantitative); qualitative "no or less inhibition" also reported | 5 quantitative + 1 qualitative primary source | Yes, in 4 of 5 (Djeujo 2022, Chen 2022, Wu 2017, Li 2009 via MeSH) | Medium-high |
| Quercetin | 2.81 - 117 uM (6.10 ug/mL - 117 uM) | 6 | Yes, in 3 of 6 (He 2019, Nguyen 2023, Li 2009) | Medium-high |
| Kaempferol | 230 uM (single value) | 1 | Yes (Nguyen 2023) | High for the one value, but single-source |
| Luteolin | 27.22 ug/mL (95.1 uM) - 172 uM, i.e. ~32-172 uM range | 4 | Yes, in 2 of 4 (Djeujo 2022, He 2019) | Medium-high |
| Rutin | 13.19 uM - >819 uM (>500 ug/mL) | 3 | Yes, in 2 of 3 (Li 2009, He 2019) | Medium (huge spread) |
| Vitexin | 50.11 - 117.6 uM | 3 | Not independently confirmed in any of the 3 | Medium |
| Isovitexin | 46.3 - 116.14 uM | 2 | Not independently confirmed in either | Medium (sparse, as expected) |
| Spinosin | Not found | 0 | N/A | N/A |
| Betulinic acid | 10.6 - 16.83 uM | 2 | Yes, in 1 of 2 (Chen 2022) | Medium-high |
| Ursolic acid | 5.08 - 466.4 uM | 5 | Yes, in 2 of 5 (Ding 2018, Wu 2017) | Medium (huge spread) |
| Oleanolic acid | 6.35 - 96.6 uM | 4 | Yes, in 1 of 4 (Ding 2018) | Medium |

**Usable yeast-assay data exists for 10 of 11 compounds.** Only spinosin has zero alpha-glucosidase data of any kind.

---

## The acarbose finding

This is the central methodological finding and it is exactly what the brief predicted: **acarbose is a weak inhibitor of yeast alpha-glucosidase**, with IC50 values in the hundreds of micromolar in every yeast assay that reported a number, while several of the plant polyphenols and triterpenes in this list inhibit the same yeast enzyme in the single-digit-to-low-double-digit micromolar range in the very same papers.

Quantitative yeast acarbose IC50 values found:
- 815.4 uM (Djeujo et al. 2022, Pharmaceuticals)
- 841.3 uM (Chen et al. 2022, Molecules)
- 91 uM (Li et al. 2009, J Agric Food Chem) - the one outlier, still 5-9x weaker than most polyphenols in this list
- 780.2 uM (Hong et al. 2013, Chinese Medicine - medium confidence, could not verify full text)
- 569.43 uM (Wu et al. 2017, Scientific Reports)

Qualitative primary confirmation: Oki, Matsui & Osajima (1999, J Agric Food Chem) state plainly: **"Voglibose, acarbose and glucono-1,5-lactone strongly inhibited mammalian AGHs, whereas no or less inhibition was observed in yeast AGH."** This is a direct, primary-literature statement of the exact asymmetry the task description anticipated.

Implication for the thesis: in a docking study built around the yeast/3A4A structure, "beats acarbose" is a low bar in this specific assay system, because acarbose itself performs poorly against yeast alpha-glucosidase. A docking-predicted compound scoring better than acarbose against 3A4A is not obviously more remarkable than the fact that betulinic acid, ursolic acid, oleanolic acid, luteolin, and quercetin all already beat acarbose in real wet-lab yeast assays by one to two orders of magnitude. This should be stated explicitly if "beats acarbose" is used anywhere as a claim of significance.

---

## Per-compound detail

### Acarbose

1. **815.4 uM** - yeast, *Saccharomyces cerevisiae* type I (Sigma-Aldrich), pNPG 2 mM, pH 6.8, 37 degrees C, enzyme 0.075 uM.
   Djeujo FM, Ragazzi E, Urettini M, Sauro B, Cichero E, Tonelli M, Froldi G. (2022). "Magnolol and Luteolin Inhibition of alpha-Glucosidase Activity: Kinetics and Type of Interaction Detected by In Vitro and In Silico Studies." *Pharmaceuticals (Basel)* 15(2):205. DOI: 10.3390/ph15020205. PMC8880268.
   Quote: "The IC50 values of magnolol and luteolin were 32.6 uM and 32.3 uM, respectively, 25 times lower than those of acarbose." Confidence: **high**.

2. **841.3 +/- 29.6 uM** - yeast alpha-glucosidase (Sigma-Aldrich), pNPG 2.5 mM, pH 6.8, 37 degrees C.
   Chen S, Lin B, Gu J, et al. (2022). "Binding Interaction of Betulinic Acid to alpha-Glucosidase and Its Alleviation on Postprandial Hyperglycemia." *Molecules* 27(8):2517. DOI: 10.3390/molecules27082517.
   Quote: "Acarbose ... IC50 value of 841.3 +/- 29.6 uM," ~50-fold weaker than betulinic acid. Confidence: **high**.

3. **91 uM** (reported as 0.091 mmol/L) - yeast, *Saccharomyces cerevisiae* (organism identified via PubMed MeSH indexing, not an explicit abstract sentence).
   Li YQ, Zhou FC, Gao F, Bian JS, Shan F. (2009). "Comparative evaluation of quercetin, isoquercetin and rutin as inhibitors of alpha-glucosidase." *J Agric Food Chem* 57(24):11463-8. DOI: 10.1021/jf903083h. PMID: 19938837.
   Quote: acarbose (control): 0.091 mmol/L, alongside quercetin 0.017, isoquercetin 0.185, rutin 0.196 mmol/L. Confidence: **medium**.

4. **780.2 +/- 1.04 uM** - enzyme source stated as yeast microplate assay in the abstract-level summary, but full text was behind a login/CAPTCHA wall and could not be independently verified.
   Hong HC, Li SL, Zhang XQ, Ye WC, Zhang QW. (2013). "Flavonoids with alpha-glucosidase inhibitory activities and their contents in the leaves of Morus atropurpurea." *Chinese Medicine* 8:19. DOI: 10.1186/1749-8546-8-19. PMC4016240.
   Confidence: **medium**.

5. **569.43 +/- 9.98 uM** - yeast, baker's yeast (*Saccharomyces cerevisiae*), pNPG 1 mmol/L, pH 6.8, 37 degrees C.
   Wu PP, et al. (2017). "Synthesis and biological evaluation of novel ursolic acid analogues as potential alpha-glucosidase inhibitors." *Scientific Reports* 7:45578. DOI: 10.1038/srep45578. PMC5372089.
   Confidence: **high**.

Qualitative primary source: Oki T, Matsui T, Osajima Y. (1999). "Inhibitory effect of alpha-glucosidase inhibitors varies according to its origin." *J Agric Food Chem* 47(2):550-3. DOI: 10.1021/jf980788t. PMID: 10563931.
Quote: "Voglibose, acarbose and glucono-1,5-lactone strongly inhibited mammalian AGHs, whereas no or less inhibition was observed in yeast AGH." Same paper: (+)-catechin IC50 against yeast AGH = 1.3x10^-1 mM (130 uM); voglibose IC50 against yeast AGH = 2.6x10^-2 mM (26 uM). No numeric acarbose IC50 given in the retrievable abstract. Confidence: **high**.

Excluded (organism unconfirmed): 139.4 ug/mL (215.9 uM if converted at MW 645.60) from Subhan et al. 2025, *J Pharm Pharmacogn Res* 13(1):311-323, DOI 10.56499/jppres24.2081_13.1.311 - the PDF could not be rendered to confirm which enzyme source was used, so this number is reported here only as a flagged, low-confidence data point and excluded from the range above.

---

### Quercetin

1. **6.10 ug/mL** (= 20.18 uM at MW 302.236 g/mol) - yeast, *Saccharomyces cerevisiae*, pNPG, pH 6.8, 37 degrees C.
   He C, Liu X, Jiang Z, Geng S, Ma H, Liu B. (2019). "Interaction Mechanism of Flavonoids and alpha-Glucosidase: Experimental and Molecular Modelling Studies." *Foods* 8(9):355. DOI: 10.3390/foods8090355. PMID: 31438605.
   Confidence: **high**.

2. **117 +/- 1.9 uM** - yeast, *Saccharomyces cerevisiae*, pNPG.
   Nguyen NH, Tran NMA, Duong TH, Vo GV. (2023). "alpha-Glucosidase inhibitory activities of flavonoid derivatives isolated from Bouea macrophylla: in vitro and in silico studies." *RSC Advances* 13(12):8190-8201. DOI: 10.1039/D3RA00650F.
   Quote: "Quercetin (compound 5) showed an IC50 of 117 +/- 1.9 uM." Confidence: **high**.

3. **17 uM** (reported as 0.017 mmol/L) - yeast, *Saccharomyces cerevisiae* (MeSH-indicated).
   Li YQ, et al. (2009). J Agric Food Chem 57(24):11463-8. DOI: 10.1021/jf903083h. Confidence: **medium**.

4. **28.7 +/- 1.2 uM** - enzyme organism not independently confirmed (standard pNPG protocol).
   Park MJ, Kang YH. (2020). "Isolation of Isocoumarins and Flavonoids as alpha-Glucosidase Inhibitors from Agrimonia pilosa L." *Molecules* 25(11):2572. DOI: 10.3390/molecules25112572.
   Confidence: **medium**.

5. **0.85 ug/mL** (= 2.81 uM at MW 302.236) - enzyme organism could not be confirmed (PDF unreadable by this search).
   Subhan M, Sanachai K, Sungthong B, Datham S, Ratha J, Puthongking P. (2025). "Comparison in vitro and in silico studies of phenolic acids and flavonoids on alpha-glucosidase inhibition." *J Pharm Pharmacogn Res* 13(1):311-323. DOI: 10.56499/jppres24.2081_13.1.311.
   Quote: "quercetin (QE; IC50 = 0.85 ug/mL) ... than acarbose (IC50 = 139.4 ug/mL)." Confidence: **low**.

6. **7.97 +/- 0.89 uM** - bibliographic metadata incomplete (full author list, volume, DOI could not be confirmed beyond the ScienceDirect identifier); title "Quercetin analogs as alpha-glucosidase inhibitors with antidiabetic activity" (2024), possibly *Food Bioscience* (ScienceDirect PII S2212429224001433, unconfirmed journal). Full text returned HTTP 403. Confidence: **low** - flagged for the reviewer to verify independently before citing in the thesis.

---

### Kaempferol

1. **230 +/- 2.7 uM** - yeast, *Saccharomyces cerevisiae*, pNPG (4-nitrophenyl beta-D-glucopyranoside).
   Nguyen NH, Tran NMA, Duong TH, Vo GV. (2023). "alpha-Glucosidase inhibitory activities of flavonoid derivatives isolated from Bouea macrophylla: in vitro and in silico studies." *RSC Advances* 13(12):8190-8201. DOI: 10.1039/D3RA00650F.
   Quote: "The known compound 3 was identified as kaempferol, which showed a potent inhibitory activity on alpha-glucosidase with an IC50 value of 230 +/- 2.7 uM." Confidence: **high**.

Only one compound-specific yeast value was found for kaempferol despite it being a very common dietary flavonol - notably sparser than quercetin or luteolin. Group-level context (not compound-specific, so not counted in the range above): Tadera K, Minami Y, Takamatsu K, Matsuoka T. (2006). *J Nutr Sci Vitaminol (Tokyo)* 52(2):149-53. DOI: 10.3177/jnsv.52.149. PMID: 16802696. Quote: "Yeast alpha-glucosidase showed potent inhibition by the anthocyanidin, isoflavone and flavonol groups with the IC50 values less than 15 microM." Kaempferol is a flavonol, so this suggests the group average is far more potent than the one compound-specific value found (230 uM) - a discrepancy worth flagging rather than resolving, since the individual per-compound table from Tadera 2006 could not be retrieved (paywalled, only the abstract was accessible).

---

### Luteolin

1. **32.3 +/- 1.17 uM** - yeast, *Saccharomyces cerevisiae* type I, pNPG 2 mM, pH 6.8, 37 degrees C.
   Djeujo FM, et al. (2022). *Pharmaceuticals (Basel)* 15(2):205. DOI: 10.3390/ph15020205.
   Quote: "The IC50 values of magnolol and luteolin were 32.6 uM and 32.3 uM, respectively, 25 times lower than those of acarbose." Confidence: **high**.

2. **27.22 ug/mL** (= 95.1 uM at MW 286.239) - yeast, *Saccharomyces cerevisiae*, pNPG, pH 6.8, 37 degrees C.
   He C, et al. (2019). *Foods* 8(9):355. DOI: 10.3390/foods8090355. Confidence: **high**.

3. **65.8 +/- 1.9 uM** - enzyme organism not independently confirmed.
   Park MJ, Kang YH. (2020). *Molecules* 25(11):2572. DOI: 10.3390/molecules25112572. Confidence: **medium**.

4. **172 uM** (reported as (1.72 +/- 0.05) x 10^-4 mol/L) - enzyme organism not explicitly stated in the retrievable abstract.
   Yan J, Zhang G, Pan J, Wang Y. (2014). "alpha-Glucosidase inhibition by luteolin: kinetics, interaction and molecular docking." *Int J Biol Macromol* 64:213-23. DOI: 10.1016/j.ijbiomac.2013.12.007. PMID: 24333230.
   Quote: "Luteolin reversibly inhibited alpha-glucosidase in a noncompetitive manner with an IC50 value of (1.72 +/- 0.05) x 10^-4 mol/L." Confidence: **medium**.

---

### Rutin

1. **196 uM** (reported as 0.196 mmol/L) - yeast, *Saccharomyces cerevisiae* (MeSH-indicated).
   Li YQ, et al. (2009). *J Agric Food Chem* 57(24):11463-8. DOI: 10.1021/jf903083h. Confidence: **medium**.

2. **>500 ug/mL** (> 819 uM at MW 610.517) - yeast, *Saccharomyces cerevisiae*, pNPG, pH 6.8, 37 degrees C. Rutin was among the weakest of 15 flavonoids tested (alongside naringin, hesperidin, baicalin, all >500 ug/mL).
   He C, et al. (2019). *Foods* 8(9):355. DOI: 10.3390/foods8090355. Confidence: **high**.

3. **13.19 +/- 1.10 uM** - enzyme source described as a yeast microplate assay in the abstract-level summary but not independently confirmed (full text blocked).
   Hong HC, Li SL, Zhang XQ, Ye WC, Zhang QW. (2013). *Chinese Medicine* 8:19. DOI: 10.1186/1749-8546-8-19. Confidence: **medium**.

**This is the widest spread found for any compound in this dataset: roughly 13 uM to over 800 uM, a ~60-fold range**, even though at least two of the three sources both describe the enzyme as yeast alpha-glucosidase with pNPG. This is a striking illustration of the assay-variability warning in the task brief and should be reported as a finding in its own right, not averaged away.

---

### Vitexin

1. **52.80 +/- 1.65 uM** - enzyme likely yeast (docking discussion references the MAL12 active site of *Saccharomyces cerevisiae*), but the abstract text itself could not be directly fetched to confirm the assay enzyme.
   Ni M, Hu X, Gong D, Zhang G. (2020). "Inhibitory mechanism of vitexin on alpha-glucosidase and its synergy with acarbose." *Food Hydrocolloids* 105:105824. DOI: 10.1016/j.foodhyd.2020.105824. Confidence: **medium**.

2. **117.6 +/- 3.6 uM** - enzyme organism not independently confirmed (standard pNPG protocol).
   Park MJ, Kang YH. (2020). *Molecules* 25(11):2572. DOI: 10.3390/molecules25112572. Confidence: **medium**.

3. **50.11 uM** - enzyme organism not independently confirmed.
   Lv R, Liu J, Li S, Gong D, Wang L, Yuan X, Chen X, Li Y. (2025). "Discovery of vitexin, isovitexin and catechin as hypoglycemic factors in mung bean via metabolomics combined with in vitro experiments." *Food Production, Processing and Nutrition* 7:41. DOI: 10.1186/s43014-025-00318-z.
   Quote: "Vitexin inhibited alpha-glucosidase activity via an uncompetitive inhibition mechanism (IC50 = 50.11 uM)." Confidence: **medium**.

---

### Isovitexin

1. **46.3 +/- 1.7 uM** - enzyme organism not independently confirmed (standard pNPG protocol).
   Park MJ, Kang YH. (2020). *Molecules* 25(11):2572. DOI: 10.3390/molecules25112572.
   Quote: "Compounds 1-8 exhibited significant dose-dependent AGI potential, with IC50 values of 24.2-117.6 uM." Confidence: **medium**.

2. **116.14 uM** - enzyme organism not independently confirmed.
   Lv R, et al. (2025). *Food Production, Processing and Nutrition* 7:41. DOI: 10.1186/s43014-025-00318-z.
   Quote: "isovitexin functioned through a competitive inhibition mechanism (IC50 = 116.14 uM)." Confidence: **medium**.

As anticipated in the task brief, isovitexin data is sparse - only two sources located, versus five or six for quercetin/luteolin. Data exists, but the base is thin.

---

### Spinosin

**No alpha-glucosidase IC50 data found**, against yeast or any other enzyme source, after multiple targeted searches (direct compound searches, Ziziphus jujuba var. spinosa-focused searches, and searches for spinosin alongside jujuboside and magnoflorine, the other major Ziziphi Semen constituents). Published spinosin literature concentrates on sedative/hypnotic activity (potentiation of pentobarbital-induced sleep via a serotonergic mechanism), anxiolytic, neuroprotective, cardioprotective, anti-melanogenic, and general pharmacokinetic/toxicity characterization. Alpha-glucosidase inhibition does not appear to have been tested and published for this specific compound. **Write this down as "not found" - do not substitute a value from a structurally similar flavone (e.g. its aglycone or isovitexin) for spinosin itself.**

---

### Betulinic acid

1. **10.6 uM** (reported as (1.06 +/- 0.02) x 10^-5 mol/L) - enzyme organism not explicit in the retrievable abstract, but the same research group's companion 2018 paper on oleanolic/ursolic acid used an explicitly-yeast assay of identical design, making yeast probable but not confirmed for this specific paper.
   Ding H, Wu X, Pan J, Hu X, Gong D, Zhang G. (2018). "New Insights into the Inhibition Mechanism of Betulinic Acid on alpha-Glucosidase." *J Agric Food Chem* 66(27):7065-7075. DOI: 10.1021/acs.jafc.8b02992. PMID: 29902001.
   Quote: "stronger inhibition of alpha-glucosidase than acarbose" (mixed-type inhibition; no numeric acarbose IC50 given in the abstract). Confidence: **medium**.

2. **16.83 +/- 1.16 uM** - yeast, yeast alpha-glucosidase (Sigma-Aldrich), pNPG 2.5 mM, pH 6.8, 37 degrees C. Acarbose IC50 = 841.3 +/- 29.6 uM in the same assay (~50-fold weaker than betulinic acid).
   Chen S, Lin B, Gu J, et al. (2022). *Molecules* 27(8):2517. DOI: 10.3390/molecules27082517. Confidence: **high**.

---

### Ursolic acid

1. **16.9 uM** (reported as (1.69 +/- 0.03) x 10^-5 mol/L) - yeast, *Saccharomyces cerevisiae* (MeSH-indicated).
   Ding H, Hu X, Xu X, Zhang G, Gong D. (2018). "Inhibitory mechanism of two allosteric inhibitors, oleanolic acid and ursolic acid on alpha-glucosidase." *Int J Biol Macromol* 107(Pt B):1844-1855. DOI: 10.1016/j.ijbiomac.2017.10.040. PMID: 29030193. Confidence: **medium**.

2. **12.1 +/- 1.0 uM** - enzyme organism not stated in the abstract.
   Zhang BW, Xing Y, Wen C, Yu XX, Sun WL, Xiu ZL, Dong YS. (2017). "Pentacyclic triterpenes as alpha-glucosidase and alpha-amylase inhibitors: Structure-activity relationships and the synergism with acarbose." *Bioorg Med Chem Lett* 27(22):5065-5070. DOI: 10.1016/j.bmcl.2017.09.027. PMID: 28964635. Confidence: **medium**.

3. **5.08 +/- 0.70 uM** - yeast, baker's yeast (*Saccharomyces cerevisiae*), pNPG 1 mmol/L, pH 6.8, 37 degrees C. Acarbose IC50 = 569.43 +/- 9.98 uM in the same assay.
   Wu PP, et al. (2017). *Scientific Reports* 7:45578. DOI: 10.1038/srep45578. Confidence: **high**.

4. **3.1 ug/mL** (= 6.79 uM at MW 456.70) - enzyme organism not stated in the abstract.
   Poongunran J, Perera HKI, Jayasinghe L, Fernando IT, Sivakanesan R, Araya H, Fujimoto Y. (2017). "Bioassay-guided fractionation and identification of alpha-amylase inhibitors from Syzygium cumini leaves." *Pharmaceutical Biology* 55(1):206-211. DOI: 10.1080/13880209.2016.1257031. PMID: 27927056. Confidence: **medium**.

5. **0.213 +/- 0.042 mg/mL** (= 466.4 uM at MW 456.70) - enzyme organism not stated in the abstract. This is a substantial outlier relative to the other four ursolic acid values.
   Wang J, Zhao J, Yan Y, Liu D, Wang C, Wang H. (2020). "Inhibition of glycosidase by ursolic acid: in vitro, in vivo and in silico study." *J Sci Food Agric* 100(3):986-994. DOI: 10.1002/jsfa.10098. PMID: 31650545. Confidence: **medium**.

**Ursolic acid shows an approximately 90-fold spread (5.08 to 466.4 uM)** - comparable to rutin's spread, and the second-widest range in this dataset. Only two of five sources confirm the enzyme as yeast; the true yeast-only range may be narrower, but this cannot be established from the sources located.

---

### Oleanolic acid

1. **6.35 uM** (reported as (6.35 +/- 0.02) x 10^-6 mol/L) - yeast, *Saccharomyces cerevisiae* (MeSH-indicated).
   Ding H, Hu X, Xu X, Zhang G, Gong D. (2018). *Int J Biol Macromol* 107(Pt B):1844-1855. DOI: 10.1016/j.ijbiomac.2017.10.040. Confidence: **medium**.

2. **35.6 +/- 2.6 uM** - enzyme organism not stated.
   Zhang BW, et al. (2017). *Bioorg Med Chem Lett* 27(22):5065-5070. DOI: 10.1016/j.bmcl.2017.09.027. Confidence: **medium**.

3. **10.11 +/- 0.30 uM** - enzyme organism not stated.
   Castellano JM, Guinda A, Macias L, Santos-Lozano JM, Lapetra J, Rada M. (2016). "Free radical scavenging and alpha-glucosidase inhibition, two potential mechanisms involved in the anti-diabetic activity of oleanolic acid." *Grasas y Aceites* 67(3):e142. DOI: 10.3989/gya.1237153.
   Quote: "OA may capture ... peroxyl radicals, and exert a strong and non-competitive inhibition of alpha-glucosidase (IC50 10.11 +/- 0.30 uM)." Confidence: **medium**.

4. **44.1 ug/mL** (= 96.6 uM at MW 456.70) - enzyme organism not stated.
   Poongunran J, et al. (2017). *Pharmaceutical Biology* 55(1):206-211. DOI: 10.1080/13880209.2016.1257031. Confidence: **medium**.

Approximately 15-fold spread (6.35 to 96.6 uM).

---

## Compounds with no yeast data found

**Spinosin.** Zero alpha-glucosidase IC50 values of any kind (yeast, rat, human, or other) were located for spinosin despite multiple targeted search strategies. See the detailed note under "Spinosin" above.

All other 10 compounds have at least one yeast-attributed or yeast-plausible value.

---

## Caveats

1. **Assay variability is real and large, not a search artifact.** For several compounds (rutin, ursolic acid, and to a lesser extent quercetin and luteolin), independently published IC50 values for the same compound against nominally the same enzyme (yeast alpha-glucosidase, pNPG substrate) span one to two orders of magnitude. Rutin ranges from ~13 uM to over 800 uM (~60-fold); ursolic acid from ~5 uM to ~466 uM (~90-fold). This is consistent with the well-documented sensitivity of alpha-glucosidase inhibition assays to enzyme batch/supplier, enzyme concentration, substrate concentration, pre-incubation time, temperature, buffer pH, and even how percent inhibition is fit to derive IC50. Any single value quoted in the thesis as "the" IC50 for a compound should be flagged as one point drawn from a wide published range, not a fixed physical constant.

2. **Yeast and mammalian (rat intestinal, human, porcine) alpha-glucosidase values are not interchangeable and were not pooled anywhere in this document.** The literature is explicit that acarbose and related iminosugar/pseudo-tetrasaccharide inhibitors are potent against mammalian intestinal alpha-glucosidase but weak against the yeast enzyme, while many plant polyphenols show the reverse pattern of relative potency, or at least a much smaller organism-dependent swing. Mixing the two enzyme sources in a single ranking or correlation would be a category error. Every value in this document was checked (where possible) against the enzyme organism stated in the source, and values whose organism could not be confirmed are explicitly flagged as such throughout (they are the majority of the medium-confidence entries).

3. **Enzyme-organism confirmation was frequently blocked by paywalls, CAPTCHAs, or PDF-rendering failures**, not just an oversight. Several important primary papers (e.g. Ni et al. 2020 Food Hydrocolloids on vitexin, the Morus atropurpurea Chinese Medicine 2013 paper on rutin, the Park & Kang 2020 Molecules paper on isovitexin/vitexin/kaempferide) could only be accessed through search-engine summaries or partial abstracts rather than full methods sections, so the specific enzyme organism (Saccharomyces cerevisiae strain, supplier, catalog number) could not be independently verified even though the numeric IC50 value itself was confirmed as belonging to that specific paper. These are marked "medium confidence" throughout, meaning: the number is real and attributable to a specific citable paper, but the enzyme-organism claim rests on the paper's likely standard practice (nearly all use commercial S. cerevisiae alpha-glucosidase with pNPG) rather than a directly-quoted confirmation.

4. **A small number of numeric claims surfaced by search were discarded outright** because no specific citable source paper could be identified for them after follow-up searching (see "Discarded / uncitable values" below). This includes one rutin value (0.037 uM) that would have been a dramatic outlier - about 350-fold more potent than any other rutin value found - which is exactly the kind of number that should not enter a thesis without a traceable citation.

5. **Kaempferol and, to a lesser extent, betulinic acid, oleanolic acid, and ursolic acid rest on very few independent sources** (one to two for kaempferol; two to four for the triterpenes). A "range" built from two or three points is not the same evidentiary weight as the five-to-six-point ranges for quercetin, luteolin, or ursolic acid. Treat single-source values (kaempferol in particular) with appropriate caution in any correlation analysis - a single yeast value of 230 uM for kaempferol is not corroborated by a second independent lab.

6. **Group-level statements from review-style primary papers (e.g. Tadera et al. 2006) were not converted into compound-specific numbers.** Tadera 2006 states that the flavonol group as a whole inhibits yeast alpha-glucosidase with IC50 <15 uM, and kaempferol is a flavonol - but the only compound-specific kaempferol value actually retrieved (230 uM, Nguyen et al. 2023) is over 15-fold higher than that group ceiling. This discrepancy is reported, not resolved; the Tadera 2006 full table with individual compound values could not be accessed (paywalled at J-STAGE), so it is not possible to say whether kaempferol was even one of the flavonols driving that <15 uM group statement.

7. **Units were preserved from source and converted transparently where useful, never silently.** Every ug/mL to uM conversion in this document states the molecular weight used (see `ic50.json`, "molecular_weights_used_for_conversion_g_per_mol"). Standard PubChem molecular weights were used: quercetin 302.236, rutin 610.517, kaempferol/luteolin 286.239 (isomers, same formula), vitexin/isovitexin 432.38 (isomers, same formula), betulinic/ursolic/oleanolic acid 456.70 (isomers, same formula), acarbose 645.60 g/mol.

8. **This document does not compute or suggest any correlation between docking scores and the IC50 values gathered here.** That step was explicitly out of scope and was not performed.

---

## Discarded / uncitable values

These numeric claims surfaced during search but could not be traced to a specific, identifiable primary paper after follow-up searching, so per the no-invention rule they were excluded from the tables and ranges above rather than reported as fact:

- **Rutin, 0.037 uM** - claimed via a search-engine summary ("Rutin had a stronger inhibition of alpha-glucosidase (IC50 = 0.037 uM) activities than quercetin"), no source paper identified. Would be a ~350-fold outlier versus every other rutin value found - discarded rather than reported.
- **Acarbose, 177.47 ug/mL and 200 ug/mL** - claimed via an unattributed search-engine synthesis, no source paper identified. Discarded.
- **Isovitexin, 68.71 ug/mL** - claimed via an unattributed search-engine synthesis, no source paper identified. Discarded.
- **Quercetin, 4.80 ug/mL** - search summary said this value "appears in multiple research studies" without naming one; could not attribute to a specific paper. Discarded.

If any of these turn out to matter later, they would need to be re-derived from a named, retrievable paper before being cited.
