# spink1-gbm-hypoxia

This repository has the code for the bioinformatic part of our article:

**SPINK1 Links Glioblastoma to Hypoxia-Associated Transcriptional Programs**

All analyses use public data. The cell culture experiments of the article (siRNA silencing and RT-qPCR in U87MG cells, Figure 6) are not part of this repository.

The result tables are not stored here, they are in the Supplementary Tables file of the article. When you run the notebooks, the same tables and figures are written again to the `results/` folder.

Some analyses were added or changed during peer review: Tables S15b and S15c, Figure 5, Table S16, the analyses without tumors of unknown IDH status, the Freedman-Lane permutation test in the adjusted models (Table S10) and the two-sided p-value of the random gene set test (Table S14).


## Main findings

- SPINK1 RNA is detected in 61% of primary IDH-wildtype glioblastomas in TCGA. In the other 39% it is not detected at all.
- Of the 50 hallmark pathway scores, the hypoxia score has the strongest link with SPINK1, before adjustment and after adjustment for the myeloid score. The link stays when we adjust for the myeloid, proliferation and mesenchymal scores together. These scores come from the same bulk RNA data, so this does not prove that the link is fully independent of tissue composition.
- Two hypoxia signatures made by other groups (Buffa, Winter) give the same result.
- The result repeats in two independent CGGA cohorts, also gene by gene.
- Removing the tumors with unknown IDH status does not change these results.
- We did not detect an association between SPINK1 and overall survival in any of the three glioblastoma cohorts.


## Files

| File | What it does | Output |
|---|---|---|
| `spink1_utils.py` | Gene panels and small functions used by all notebooks | - |
| `01_cohort.ipynb` | TCGA cohort and expression matrix | Tables S1, S2 |
| `02_cell_context.ipynb` | SPINK1 distribution, cell type scores, list of all gene panels | Tables S3, S16, Figure 1 |
| `03_transcriptome_gsea.ipynb` | Correlation of all genes with SPINK1, GSEA | Tables S4, S5 |
| `04_pathway_adjustment.ipynb` | Pathway scores, adjustment, transcriptional subtypes | Tables S6-S8, Figure 2 |
| `05_hypoxia_validation.ipynb` | Published hypoxia signatures, robustness tests | Tables S9, S10, Figure 3 |
| `06_cgga_replication.ipynb` | Replication in two CGGA cohorts, controls | Tables S11-S14, Figure 4 |
| `07_survival.ipynb` | Overall survival in the three cohorts | Table S15 (a-c), Figure 5 |


## Data

Nothing is stored in this repository. All files are free to download.

**TCGA (discovery cohort).** Glioblastoma PanCancer Atlas study on cBioPortal, file `gbm_tcga_pan_can_atlas_2018.tar.gz`
(https://www.cbioportal.org/study/summary?id=gbm_tcga_pan_can_atlas_2018).
We use RNA-seq (RSEM), clinical data, mutation calls, copy number (GISTIC) and the hypoxia scores (Buffa, Winter, Ragnum) that come with this study. We downloaded this file on 24 August 2026.

**CGGA (validation cohorts).** From http://www.cgga.org.cn/download.jsp:

- `CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip`
- `CGGA.mRNAseq_693_clinical.20200506.txt.zip`
- `CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip`
- `CGGA.mRNAseq_325_clinical.20200506.txt.zip`

We downloaded these files on 27 July 2026. cBioPortal and CGGA can update their files later, so our results come from the files and dates written here.

Why read counts? SPINK1 is missing from the normalized expression files of CGGA, because genes with low expression were removed from those files. The read count files keep all genes, so we normalized the counts ourselves (notebook 06).

**Gene sets.** The 50 MSigDB hallmark gene sets, downloaded inside the notebooks with gseapy (library `MSigDB_Hallmark_2020`). All other gene panels are in `spink1_utils.py` and in Table S16.


## How to run

We wrote the notebooks for Google Colab. Put the files in Google Drive like this:

```
MyDrive/SPINK1/
    gbm_tcga_pan_can_atlas_2018.tar.gz
    CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip
    CGGA.mRNAseq_693_clinical.20200506.txt.zip
    CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip
    CGGA.mRNAseq_325_clinical.20200506.txt.zip
    repo/            <- this repository
        results/     <- create empty, tables and figures are written here
        processed/   <- create empty, notebook 01 writes the expression matrix here
```

If you use another folder name, change `BASE` (or `REPO`) at the top of each notebook.

Run the notebooks in order, 01 to 07. Notebook 01 has to be first, because the others read its output. Notebook 06 also reads Table S4 from notebook 03.

Every random step has a fixed seed, so you get the same numbers in every run.

We used Python 3.13 in Google Colab with these package versions (the same as in Section 4.11 of the article): pandas 2.2.3, numpy 2.1.3, scipy 1.16.3, statsmodels 0.15.0, gseapy 1.3.2, lifelines 0.30.3 and matplotlib 3.10.0. The other packages come with Colab, and the notebooks install gseapy and lifelines with `pip` without a version number, so a new run can get newer versions (Colab also updates its own packages from time to time). To use the same versions, write `gseapy==1.3.2` and `lifelines==0.30.3` in the `pip install` line at the top of the notebooks.


## A few terms we use a lot

**Spearman correlation (ρ).** It shows if two values go up and down together, from -1 to +1. It works on ranks, not on raw values. We use it everywhere because SPINK1 has many zeros and a skewed distribution. With ranks, a few extreme values have little effect and the data do not need to be normally distributed.

**p and FDR.** If you do thousands of tests, some small p-values appear just by chance. The Benjamini-Hochberg method corrects for this. FDR < 0.05 means that among the results called significant, less than 5% are expected to be false.

**Score.** Many analyses need one number per tumor for a cell type or a pathway. Each gene of the panel is turned into a z-score across tumors (mean 0, SD 1), and the score is the mean of these z-scores. Higher score, stronger signal of that cell type or pathway.

**Partial correlation.** The correlation between two values after removing what a third value explains. For example, if SPINK1 and a pathway only go up together because some tumors have more macrophages, their partial correlation with the myeloid score as covariate drops to near zero. We use a rank based version (partial Spearman): all values are ranked, the covariates are removed by linear regression, and the rest is correlated.


## Notebooks and results

### 01_cohort: TCGA discovery cohort (Tables S1, S2)

TCGA has RNA-seq for 160 glioblastoma samples. We kept only primary tumors, one sample per patient, which gives 153 tumors. Then we removed the 9 tumors with an IDH1 or IDH2 mutation, since in the 2021 WHO classification glioblastoma is IDH-wildtype by definition. Final cohort: **144 tumors**.

3 of these 144 tumors have no mutation data in the study and no IDH-mutant subtype label. We kept them and wrote "Not available" for their IDH status in Table S1. Results without these 3 tumors are given in notebooks 05 and 07.

Expression values are log2(RSEM + 1). Negative values are set to 0, and if a gene name appears twice we keep the row with the higher mean. A gene is kept if its value is above 1 in at least 20% of all 160 samples. This leaves 17,754 genes, SPINK1 is one of them.

**SPINK1 groups.** SPINK1 is not detected (value 0) in 56 of the 144 tumors (38.9%). With so many zeros we made three groups:

- Undetectable: SPINK1 = 0 (n = 56)
- Low: above 0 and up to the median of the detectable tumors, 2.15 (n = 44)
- High: above this median (n = 44)

We use the groups in figures and some group tests. The main analyses use SPINK1 as a continuous value.

**Copy number (Table S2).** When a gene is amplified in the tumor DNA, its RNA often goes up as well, so a positive correlation between copy number and expression would suggest that copy number changes drive SPINK1. We found no relation (ρ = 0.055, p = 0.52, n = 140 tumors with copy number data).


### 02_cell_context: which cells? (Tables S3, S16, Figure 1)

Bulk RNA-seq measures a mix of cells: tumor cells, immune cells, vessels and normal brain. A gene can look linked to a pathway only because it is made by one cell type, so we looked at this first.

We made scores for 10 cell types and cell states (myeloid/macrophage, neutrophil, T/NK cell, endothelial, stromal, astrocyte, neuron, oligodendrocyte, OPC/proneural, proliferation) and corrected the 10 SPINK1 tests with FDR. To see what a strong correlation looks like, we tested three known marker genes against the same scores: CD68 (macrophage), SYT1 (neuron) and MKI67 (proliferation). Each gene is removed from its own panel first, so it is not correlated with itself.

**Result.** The markers behave as expected: CD68 with the myeloid score ρ = 0.91, SYT1 with the neuron score 0.80, MKI67 with the proliferation score 0.60. SPINK1 is positively linked with the myeloid (ρ = 0.43, FDR = 6.6 × 10⁻⁷), neutrophil (0.37) and stromal (0.30) scores, and negatively with proliferation (-0.25, FDR = 0.007). The T/NK link (0.18) has p = 0.035 but does not pass FDR correction (0.07). Neuron, astrocyte and endothelial scores show almost nothing (|ρ| < 0.04).

SPINK1 follows the myeloid score much less closely than CD68 does (0.43 vs 0.91). This suggests that SPINK1 is not just a marker of macrophage content. Still, a panel score is not a cell count, and bulk data can not show which cells actually make SPINK1. Because myeloid cells may affect the results, we adjust for the myeloid score in all later steps.

**Table S16** lists every gene panel used in the study (cell types, subtypes, the clean mesenchymal panel and the proliferation covariate) and marks which genes are present in the TCGA data. The myeloid covariate in later notebooks is the Myeloid / macrophage panel.


### 03_transcriptome_gsea: genes and pathways (Tables S4, S5)

**Correlation with every gene (Table S4).** We correlated SPINK1 with the other 17,753 expressed genes and corrected for multiple testing. 4,158 genes have FDR < 0.05 (2,581 positive, 1,577 negative). The top positive ones are inflammatory and hypoxia genes, for example IL8 (CXCL8, ρ = 0.66), TREM1 (0.62), CCL20 (0.59) and NDRG1 (0.58). The top negative ones are DNA replication and repair genes such as ATAD2 (-0.48) and MCM4 (-0.47). TCGA still uses some old gene names, IL8 is CXCL8.

**GSEA (Table S5).** Here all genes are sorted from the most positive to the most negative correlation with SPINK1. GSEA asks whether the genes of a pathway sit at the top of this list (positive enrichment) or at the bottom (negative enrichment) more than by chance. NES, the normalized enrichment score, is the size of this effect corrected for the size of the gene set, so different pathways can be compared. p and FDR come from 1,000 permutations. When no permutation reaches the real value, GSEA gives 0. We write this as "<0.001", because that is the resolution of 1,000 permutations. Settings: pre-ranked gseapy, gene set size 10 to 500, seed 42.

**Result.** 40 of the 50 hallmark sets have FDR < 0.05. The strongest positive ones are TNF-α signaling via NF-κB (NES 3.06), inflammatory response (2.96), hypoxia (2.80) and IL-6/JAK/STAT3 signaling (2.76). The strongest negative ones are E2F targets (-3.90), G2M checkpoint (-3.71) and mitotic spindle (-2.69). Tumors with more SPINK1 show more inflammatory and hypoxic activity and less proliferation.

GSEA can not tell cell composition apart from real pathway activity. Notebook 04 looks at this.

You will see a gseapy warning that 1.68% of the genes have tied values. Tied genes get an arbitrary but fixed order, so the result is the same in every run.


### 04_pathway_adjustment: is it only cell composition? (Tables S6-S8, Figure 2)

**Pathway scores (Table S6, Figure 2).** For each tumor and each hallmark pathway we made a score and correlated it with SPINK1, three times: unadjusted, adjusted for the myeloid score, and adjusted for the proliferation score (E2F targets plus G2M checkpoint genes).

- Unadjusted, 42 of 50 pathways are significant (FDR < 0.05).
- After myeloid adjustment only 14 of 50 stay. Most of the lost ones are immune pathways, for example allograft rejection (ρ = 0.37 to -0.08) and IL-6/JAK/STAT3 signaling (0.43 to 0.11). These immune signals seem to follow the amount of immune cells in the tumor.
- Hypoxia stays: 0.54 unadjusted, 0.40 after myeloid adjustment (FDR < 10⁻⁴) and 0.48 after proliferation adjustment.
- After proliferation adjustment 31 of 48 pathways stay. E2F targets and G2M checkpoint are left out here, because they are the covariate itself.
- Hypoxia is the strongest pathway before adjustment and after myeloid adjustment. After proliferation adjustment, glycolysis (0.49) and hypoxia (0.48) are together at the top. Glycolysis is a main part of the hypoxic response anyway.

**Transcriptional subtypes (Table S7).** Glioblastoma has transcriptional subtypes, and mesenchymal tumors are known to have more hypoxia and more immune cells. So the subtype could be a hidden reason for the SPINK1-hypoxia link. Each subtype score is the mean z-score of its marker genes (Verhaak et al. 2010, Wang et al. 2017). The three scores are then standardized across tumors, and every tumor gets the subtype with the highest standardized score. We did not use the neural subtype, because Wang et al. showed that it mostly reflects normal brain tissue in the sample. These short marker panels approximate the published classifiers, they are not the same classifier.

To compare SPINK1 between subtypes we used the Kruskal-Wallis test, which compares more than two groups using ranks and does not need normal data.

**Result.** 45 proneural, 51 classical and 48 mesenchymal tumors. SPINK1 differs between subtypes (Kruskal-Wallis p = 8.0 × 10⁻⁶): highest in mesenchymal and lowest in classical tumors (61% of classical tumors have no SPINK1). But inside every subtype SPINK1 is still positively correlated with all three hypoxia scores (9 of 9 tests positive). These tests use the 143 tumors with published hypoxia scores, so the subtypes have 45, 50 and 48 tumors. With the hallmark hypoxia score ρ is between 0.41 and 0.48 and p < 0.01 in each subtype. Some Buffa and Winter values do not reach p < 0.05 inside a subtype, which is not surprising with so few tumors per subtype.

**All covariates together, and one circularity problem (Tables S7, S8).** Next we adjusted for myeloid, mesenchymal and proliferation scores at the same time. There is a catch here. Four genes of the mesenchymal signature (COL5A1, LOX, PLAUR, SERPINE1) are also in the hallmark hypoxia set. Adjusting for a score that contains hypoxia genes would remove part of the hypoxia signal itself. So we also made a "clean" mesenchymal score without these four genes and report both versions.

**Result** (n = 143, tumors with published hypoxia scores): the hallmark hypoxia link is ρ = 0.38 with the original mesenchymal score and 0.43 with the clean one. Buffa gives 0.30 and 0.33, Winter 0.32 and 0.36, all with p < 0.001. So the link is not explained by these measured covariates. This does not prove that it is fully independent of tissue composition, because the covariates are themselves estimated from bulk expression. The clean score is used in all later notebooks.


### 05_hypoxia_validation: is the hypoxia score right? (Tables S9, S10, Figure 3)

We built the hallmark hypoxia score ourselves, so it could measure something else. To check this we used three hypoxia signatures made by other groups. Their scores come with the TCGA data and exist for 143 of the 144 tumors. Buffa and Winter are our main validation scores. Ragnum was made in prostate cancer and we include it for comparison.

**Result (Table S9).**

- Our hallmark score agrees with Buffa (ρ = 0.63, p = 1.9 × 10⁻¹⁷) and Winter (0.63, p = 7.5 × 10⁻¹⁷), so it picks up largely the same signal.
- SPINK1 correlates with Buffa (ρ = 0.35, p = 1.7 × 10⁻⁵) and Winter (0.36, p = 8.5 × 10⁻⁶). After myeloid adjustment it is 0.27 (p = 0.001) and 0.34 (p = 4.3 × 10⁻⁵).
- SPINK1 does not correlate with Ragnum (ρ = -0.03, p = 0.74). Ragnum also agrees only weakly with our hallmark score in these tumors (ρ = 0.21). A hypoxia signature made in one tumor type does not always work the same way in another one.
- The median Buffa score goes up step by step from Undetectable to Low to High (3, 9 and 21; Kruskal-Wallis p = 2.3 × 10⁻⁵).
- The published scores use other gene lists and other scoring methods, so their ρ values can not be compared one to one with the hallmark score.

**Robustness (Table S10).** We tested the Buffa and Winter results in several ways, both unadjusted and fully adjusted (myeloid, clean mesenchymal, proliferation):

- *Bootstrap confidence interval.* We draw 143 tumors with replacement and compute ρ again, 2,000 times. The middle 95% of these values is the 95% confidence interval. If it stays above 0, the result does not depend on which tumors happen to be in the sample. The per-tumor scores are kept fixed in this step. Result: all four intervals are above 0 (lowest bound 0.15, adjusted Buffa, 0.15 to 0.49).
- *Permutation test.* We shuffle the data 2,000 times so that any real link is broken, and count how often the shuffled ρ is as strong as the real one. Unadjusted, we shuffle SPINK1 between tumors. In the adjusted models this would also break the real link between SPINK1 and the covariates, so there we use the Freedman-Lane method: the covariates stay fixed and only the part of the hypoxia score that they do not explain is shuffled. Result: at most 1 of the 2,000 shuffles reached the real value (p ≤ 0.0005 in all four tests).
- *Leave-one-out.* We remove one tumor at a time and compute ρ again, 143 times, to see if a single tumor drives the result. ρ changes by at most 0.03.
- *Only tumors with detectable SPINK1 (n = 88).* This checks whether the result only comes from zero vs non-zero tumors. The link gets stronger, not weaker (adjusted Buffa 0.52, Winter 0.47), so the relation is graded.
- *Trimmed.* We remove SPINK1 values outside the 2.5 to 97.5 percentile range. Because 39% of the values are zero, in practice only the highest values are removed. Almost no change (ρ = 0.33 to 0.35).
- *Only IDH-wildtype tumors (n = 140).* We remove the 3 tumors with unknown IDH status. The results stay the same or get a bit stronger: Buffa 0.37 (adjusted 0.35), Winter 0.37 (adjusted 0.37). Below the table the notebook also prints the hallmark hypoxia result for these 140 tumors: ρ = 0.55 (p = 1.2 × 10⁻¹²), adjusted 0.45 (p = 2.5 × 10⁻⁸).


### 06_cgga_replication: does it repeat in other patients? (Tables S11-S14, Figure 4)

**Cohorts.** Same rules as in TCGA: primary tumor, WHO grade IV, glioblastoma histology, IDH not mutant. CGGA-693 has 116 tumors and CGGA-325 has 74, with no patient in both. In CGGA-693, 7 of the 116 tumors have no IDH information (109 are IDH-wildtype). We kept them, like the 3 TCGA tumors, and Table S11 also shows the result without them.

We analyze every cohort on its own and do not merge them. TCGA and CGGA give different data types (RSEM values vs read counts), and keeping the cohorts apart gives three independent tests without batch correction.

**Normalization.** CPM (counts per million) is the read count of a gene divided by all read counts of the sample, times one million. It corrects for different sequencing depth. Then we take log2(CPM + 1). Duplicated gene names keep the row with the higher mean count. A gene is kept if CPM > 1 in at least 20% of the tumors of the cohort (17,895 genes in CGGA-693, 18,004 in CGGA-325).

In CGGA-693 SPINK1 does not pass this filter because its level is low, but it is still detected (count above 0) in 87 of 116 tumors (75%). It is the gene we study, so we kept it. Low counts make the measurement noisier, and we can not fully rule out technical effects of this. But the same result in CGGA-325, where SPINK1 passes the filter, and in TCGA, which uses a different pipeline, makes a technical artifact less likely.

Scores are made the same way as in TCGA, inside each cohort. Hallmark hypoxia: 178 and 177 of 200 genes found. Myeloid: 9 of 11. Clean mesenchymal: 10 of 10. Proliferation: 311 of 327. SPINK1 groups follow the TCGA rule (Undetectable/Low/High: 29/44/43 and 16/29/29).

**Replication (Table S11, Figure 4).**

- Unadjusted: ρ = 0.59 (p = 4.8 × 10⁻¹²) in CGGA-693 and 0.63 (p = 1.3 × 10⁻⁹) in CGGA-325. In TCGA it was 0.54.
- Fully adjusted (myeloid, clean mesenchymal, proliferation): 0.21 (p = 0.025) and 0.45 (p = 8.8 × 10⁻⁵). In TCGA 0.43. The adjusted value is smaller in CGGA-693, but it stays positive and significant.
- Only tumors with detectable SPINK1: 0.57 (n = 87) and 0.62 (n = 58), similar to the full cohorts.
- CGGA-693 without the 7 tumors of unknown IDH status (n = 109): 0.62 (p = 6.4 × 10⁻¹³), adjusted 0.20 (p = 0.044).
- The 95% confidence intervals in Figure 4D use the Fisher z method with the variance correction for Spearman correlation of Fieller, Hartley and Pearson (1957): SE = sqrt(1.06 / (n - 3 - k)), where k is the number of covariates.

**Gene by gene (Table S12).** If the link is real, the same hypoxia genes should go with SPINK1 in every cohort. For each hallmark hypoxia gene we computed its correlation with SPINK1 in each cohort. 176 genes have data in all three, and we compared these 176 values between cohorts.

The gene profiles agree well (ρ = 0.65 between TCGA and CGGA-693, 0.86 between TCGA and CGGA-325, 0.73 between the two CGGA cohorts). 130 of the 176 genes are positive in all three cohorts. Of the 89 genes with FDR < 0.05 in TCGA, 81 (91%) are positive in all three. The most consistent ones are known HIF targets: NDRG1, LOX, HMOX1, ANGPTL4 and PLIN2.

**Pipeline controls (Table S13).** Before trusting a new result, we checked that the same pipeline finds things we already know. A positive control gene should correlate clearly with its own score, a negative control should not. Each gene is removed from the score before the test.

- CA9, the classic HIF target, correlates with the hypoxia score (0.68 and 0.75). CA9 is not in the hallmark hypoxia set, so this is an outside check.
- PTPRC (CD45) with the myeloid score: 0.67 and 0.78. CD68 is not in the CGGA gene list, so we used PTPRC.
- MKI67 and TOP2A with the proliferation score: 0.84 to 0.93.
- CHI3L1 (YKL-40) with the mesenchymal score: 0.73 and 0.62.
- Neuronal genes (SNAP25, SYT1, RBFOX3) with the hypoxia score: close to zero and not significant in CGGA-693 (0.06 to 0.16), negative in CGGA-325 (-0.26 to -0.40).
- SPINK1 itself, for comparison: 0.59 and 0.63.

**Random gene sets (Table S14).** In bulk tissue almost any large gene set shares some general signal, for example the share of tumor cells in the sample. So maybe SPINK1 would correlate with any big gene set, not only with hypoxia. To test this we drew 1,000 random gene sets of the same size as the hypoxia set (without hypoxia genes), made a score from each one and correlated it with SPINK1. This shows what we get by chance. The table has:

- z-score = (real ρ - mean of random ρ) / SD of random ρ, so how many standard deviations the real value is away from chance.
- `p_stronger`: share of random sets with |ρ| at least as large as the real |ρ|. It answers "is the link stronger than for random genes?".
- `p_two_sided`: two times the smaller tail of the random distribution. It answers "is the real value different from random gene sets in either direction?".

Random sets also correlate a little with SPINK1 (mean ρ = 0.24 and 0.11), as expected. The hypoxia set is much higher: 10.0 and 8.0 SD above the random mean (both p < 0.001). CA9 gives a similar picture (z = 7.9 and 10.6). The negative control SNAP25 is far below the random mean (z = -5.8 and -4.3, `p_two_sided` < 0.001 in both). In CGGA-693 its `p_stronger` is 1.000, because every random set correlates at least as strongly with SNAP25 as the hypoxia set does, which is what we expect from a negative control. Random sets are matched only by size, not by expression level or gene-gene correlation, so this test shows that the hypoxia link is much stronger than for size-matched random genes. It is not a proof that hypoxia is the only program linked to SPINK1.


### 07_survival: overall survival (Table S15, Figure 5)

**Survival terms.**

- *Kaplan-Meier curve.* The share of patients still alive over time. Every step down is a death. The small vertical marks are censored patients, alive at the last follow-up, so their full survival time is not known.
- *Log-rank test.* Compares whole survival curves between groups.
- *Cox model and hazard ratio (HR).* HR above 1 means higher risk of death, below 1 lower risk, 1 means no effect. If the 95% confidence interval includes 1, the effect is not significant. A Cox model can include other factors, like age and sex, at the same time.
- *HR per 1 SD.* SPINK1 has different units in TCGA (RSEM) and CGGA (CPM), so we divided SPINK1 by its standard deviation in each analysis. The HR is then the change in risk for one SD more SPINK1, and cohorts can be compared.
- *Proportional hazards test.* The Cox model assumes that the effect of SPINK1 is the same during the whole follow-up. We tested this for SPINK1 with Schoenfeld residuals. p > 0.05 means no sign that the assumption is broken.

**Table S15a: the TCGA analysis of the first submission.** In TCGA, age and sex are missing for 45 of the 143 patients with follow-up. To compare models with and without covariates on the same patients, this part uses the 98 patients with complete data. HR per log2 unit of SPINK1: 1.00 (95% CI 0.87 to 1.15, p = 0.98). Adding age and sex, then the Buffa score, then the myeloid, mesenchymal and proliferation scores does not change this (HR 0.95 to 1.02, all p > 0.6).

**Tables S15b, S15c and Figure 5: three cohorts.**

- Cox, unadjusted: TCGA HR 1.10 (0.92 to 1.31), CGGA-693 0.89 (0.71 to 1.12), CGGA-325 1.01 (0.80 to 1.26). With age and sex: 1.04, 0.86 and 1.02. All p > 0.2. The HRs point in different directions and all intervals include 1.
- Without tumors of unknown IDH status: TCGA HR 1.12 (0.94 to 1.34, n = 140), CGGA-693 0.87 (0.70 to 1.08, n = 105). Same picture.
- Proportional hazards tests for SPINK1: all p > 0.05 (lowest 0.096).
- Kaplan-Meier with the three SPINK1 groups: log-rank p = 0.40 (TCGA), 0.50 (CGGA-693) and 0.41 (CGGA-325). Median survival for Undetectable/Low/High: 14.9/11.7/10.6 months in TCGA, 12.0/13.3/16.2 in CGGA-693 and 13.9/10.3/15.9 in CGGA-325. The order of the groups is not the same in the three cohorts.
- Best cut-off scan. We also tried every cut-off that leaves at least 20% of the patients in each group (59, 62 and 44 cut-offs). This is close to the method behind the survival plots of the Human Protein Atlas. Picking the smallest p out of many tests is too optimistic, because some small p-values show up by chance, so this scan is only exploratory. Even so, the smallest p was 0.097, 0.148 and 0.103, and no cut-off reached 0.05.
- All primary gliomas in CGGA (grade 2 to 4, both cohorts together, IDH status known, n = 588). Here more SPINK1 goes with worse survival (HR 1.26, 1.16 to 1.38, p = 2.4 × 10⁻⁷). After adding grade and IDH status to the model the association is no longer seen (HR 0.96, 0.87 to 1.06, p = 0.41). This is a simple reference model, it is not stratified by cohort.

In short, within IDH-wildtype glioblastoma we did not detect an association between SPINK1 and overall survival, in any of the three cohorts and with any of the methods above. This does not prove that there is no effect at all, but the estimates are close to 1 in all three cohorts.


## Why the number of patients changes

- **144**: the TCGA cohort.
- **143** (notebooks 04 and 05): one tumor (TCGA-14-0871) has no published hypoxia scores, so all analyses with Buffa, Winter or Ragnum use 143 tumors.
- **140** (notebooks 05 and 07): the 143 tumors without the 3 tumors of unknown IDH status.
- **143** (notebook 07): another tumor (TCGA-28-2510) has no follow-up data.
- **98** (Table S15a, and the age and sex models in Table S15b): patients with follow-up, age and sex.
- **140** (Table S2): tumors with copy number data.
- **116 and 74**: the CGGA cohorts. **109**: CGGA-693 without the 7 tumors of unknown IDH status.
- **112 and 105**: CGGA-693 in the survival analysis, because 4 patients have no follow-up data
## Contact

Code: Melisa Ağrı (agrim22@itu.edu.tr)
Corresponding author: Zeliha Emrence (zeliha.emrence@istanbul.edu.tr)
