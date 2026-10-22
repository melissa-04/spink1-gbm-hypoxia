# spink1-gbm-hypoxia

Code for the bioinformatic part of the article:

**SPINK1 Links Glioblastoma to Hypoxia-Associated Transcriptional Programs**


All analyses use public data. The cell culture part of the article (siRNA silencing and RT-qPCR in U87MG cells, Figure 5) is not in this repository.

Result tables are not stored here. They are in the Supplementary Tables file of the article. If you run the notebooks, the same tables and figures are written again to the `results/` folder.


## Main findings

- SPINK1 RNA is detected in 61% of primary IDH-wildtype glioblastomas in TCGA. In the other 39% it is not detected at all.
- Among the 50 hallmark pathways, hypoxia has the strongest link with SPINK1. The link stays after we control for myeloid cells, proliferation and mesenchymal state.
- Two hypoxia signatures published by other groups (Buffa, Winter) give the same result.
- The result repeats in two independent CGGA cohorts, also gene by gene.
- SPINK1 expression is not associated with overall survival in any of the three glioblastoma cohorts.


## Files

| File | What it does | Output |
|---|---|---|
| `spink1_utils.py` | Gene panels and small functions used by all notebooks | - |
| `01_cohort.ipynb` | Builds the TCGA cohort and the expression matrix | Tables S1, S2 |
| `02_cell_context.ipynb` | SPINK1 distribution, cell type scores | Table S3, Figure 1 |
| `03_transcriptome_gsea.ipynb` | Correlation of all genes with SPINK1, GSEA | Tables S4, S5 |
| `04_pathway_adjustment.ipynb` | Pathway scores, adjustment, transcriptional subtypes | Tables S6-S8, Figure 2 |
| `05_hypoxia_validation.ipynb` | Published hypoxia signatures, robustness tests | Tables S9, S10, Figure 3 |
| `06_cgga_replication.ipynb` | Replication in two CGGA cohorts, controls | Tables S12-S15, Figure 4 |
| `07_survival.ipynb` | Overall survival in the three cohorts | Table S11, Figure S1 |


## Data

Nothing is stored in this repository. All files can be downloaded for free.

**TCGA (discovery cohort).** Glioblastoma PanCancer Atlas study on cBioPortal, file `gbm_tcga_pan_can_atlas_2018.tar.gz`
(https://www.cbioportal.org/study/summary?id=gbm_tcga_pan_can_atlas_2018).
We use RNA-seq (RSEM), clinical data, mutation calls, copy number (GISTIC) and the hypoxia scores (Buffa, Winter, Ragnum) that come with this study.

**CGGA (validation cohorts).** From http://www.cgga.org.cn/download.jsp:

- `CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt.zip`
- `CGGA.mRNAseq_693_clinical.20200506.txt.zip`
- `CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt.zip`
- `CGGA.mRNAseq_325_clinical.20200506.txt.zip`

Why read counts? SPINK1 is not in the normalized expression files of CGGA, because genes with low expression were removed from those files. The read count files keep all genes. So we normalized the counts ourselves (notebook 06).

**Gene sets.** The 50 MSigDB hallmark gene sets. The notebooks download them with gseapy (library `MSigDB_Hallmark_2020`).


## How to run

The notebooks were written for Google Colab. Put the files in Google Drive like this:

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

Run the notebooks in order, 01 to 07. Notebook 01 must be first, the others read its output. Notebook 06 also reads Table S4 from notebook 03.

All random steps have a fixed seed, so the numbers are the same in every run.

Main packages: pandas, numpy, scipy, statsmodels, gseapy (1.3.1), lifelines, matplotlib.


## Some words we use many times

**Spearman correlation (ρ).** Shows if two values go up and down together. It goes from -1 to +1. It uses ranks, not the raw values. We use it everywhere because SPINK1 has many zeros and a skewed distribution. With ranks, a few extreme values have little influence, and we do not need normally distributed data.

**p and FDR.** When we do thousands of tests, some small p-values appear only by chance. The Benjamini-Hochberg method corrects for this. FDR < 0.05 means that among the results we call significant, less than 5% are expected to be false.

**Score.** Many analyses need one number per tumor for a cell type or a pathway. For this, each gene of the panel is changed to a z-score across tumors (mean 0, SD 1). Then we take the mean of these z-scores. A higher score means a stronger signal of this cell type or pathway.

**Partial correlation.** The correlation between two values after removing the part that a third value explains. Example: if SPINK1 and a pathway both go up only because a tumor has more macrophages, their partial correlation, controlling for the myeloid score, falls close to zero. We use a rank based version (partial Spearman): all values are ranked, the effect of the covariates is removed with linear regression, and the rest is correlated.


## Notebooks and results

### 01_cohort: TCGA discovery cohort (Tables S1, S2)

**Which tumors.** TCGA has RNA-seq for 160 glioblastoma samples. We kept primary tumors only, one sample per patient: 153 tumors. We removed the 9 tumors with an IDH1 or IDH2 mutation. In the 2021 WHO classification glioblastoma is IDH-wildtype by definition. IDH-mutant tumors are now a separate tumor type. Final cohort: **144 tumors**.

3 of the 144 tumors have no mutation data in this study and no IDH-mutant subtype label. We kept them. In Table S1 their IDH status is written as "Not available".

**Expression.** Values are log2(RSEM + 1). If a gene name appears twice, we keep the row with the higher mean. A gene is kept if its value is above 1 in at least 20% of the 160 samples. This removes genes that are almost never expressed. 17,754 genes pass, SPINK1 is one of them.

**SPINK1 groups.** SPINK1 is not detected (value 0) in 56 of 144 tumors (38.9%). Because so many values are zero, we made three groups:

- Undetectable: SPINK1 = 0 (n = 56)
- Low: above 0, up to the median of the detectable tumors (2.15) (n = 44)
- High: above this median (n = 44)

The groups are used in figures and in some group tests. The main analyses use SPINK1 as a continuous value.

**Copy number (Table S2).** If a gene is amplified in tumor DNA, its RNA often goes up too. So a positive correlation between copy number and expression would suggest that copy number changes drive the expression. Result: no relation (ρ = 0.055, p = 0.52, n = 140 tumors with copy number data). SPINK1 levels are not mainly explained by gains or losses of the gene.


### 02_cell_context: which cells? (Table S3, Figure 1)

Bulk RNA-seq measures a mix of cells: tumor cells, immune cells, vessels, normal brain. A gene can look linked to a pathway only because it comes from one cell type. So we checked this first.

We made scores for 10 cell types and cell states: myeloid/macrophage, neutrophil, T/NK cell, endothelial, stromal, astrocyte, neuron, oligodendrocyte, OPC/proneural and proliferation. The marker genes are in `spink1_utils.py`.

To know what a "strong" correlation looks like, we also tested known marker genes against the same scores: CD68 (macrophage), SYT1 (neuron) and MKI67 (proliferation). A gene is removed from its own panel before the test, so it is not correlated with itself.

**Result.**

- The markers work as expected: CD68 with myeloid score ρ = 0.91, SYT1 with neuron score 0.80, MKI67 with proliferation score 0.60.
- SPINK1 is positively linked with the myeloid (ρ = 0.43, p = 6.6 × 10⁻⁸), neutrophil (0.37) and stromal (0.30) scores. It is negatively linked with proliferation (-0.25, p = 0.003). It has almost no link with neuron, astrocyte and endothelial scores (|ρ| < 0.04).
- SPINK1 follows the myeloid score much less than CD68 does (0.43 against 0.91). So SPINK1 is not simply a macrophage gene. But myeloid cells can still affect the results, so we control for the myeloid score in all later steps.


### 03_transcriptome_gsea: genes and pathways (Tables S4, S5)

**Correlation with every gene (Table S4).** We correlated SPINK1 with all other 17,753 expressed genes and corrected for multiple testing (FDR).

Result: 4,158 genes have FDR < 0.05 (2,581 positive, 1,577 negative). The top positive genes are inflammatory and hypoxia genes, for example IL8 (CXCL8, ρ = 0.66), TREM1 (0.62), CCL20 (0.59) and NDRG1 (0.58). The top negative genes are DNA replication and repair genes, for example ATAD2 (-0.48) and MCM4 (-0.47). (TCGA uses some older gene names, IL8 is CXCL8.)

**GSEA (Table S5).** We sort all genes from the most positive to the most negative correlation with SPINK1. GSEA asks if the genes of one pathway are collected at the top of this list (positive enrichment) or at the bottom (negative enrichment), more than by chance. NES (normalized enrichment score) is the size of this effect, corrected for the size of the gene set, so pathways can be compared. FDR comes from 1,000 permutations. Settings: pre-ranked mode of gseapy, gene set size 10 to 500, seed 42.

Result: 40 of the 50 hallmark sets have FDR < 0.05. Strongest positive: TNF-α signaling via NF-κB (NES 3.06), inflammatory response (2.96), hypoxia (2.80), IL-6/JAK/STAT3 signaling (2.76). Strongest negative: E2F targets (-3.90), G2M checkpoint (-3.71), mitotic spindle (-2.69). So tumors with more SPINK1 show more inflammatory and hypoxic activity and less proliferation.

GSEA alone can not separate cell composition from real pathway activity. Notebook 04 does this.

gseapy prints a warning that 1.68% of the genes have tied values. The order of tied genes is arbitrary, but it is the same in every run, so the result does not change.


### 04_pathway_adjustment: is it only cell composition? (Tables S6-S8, Figure 2)

**Pathway scores (Table S6, Figure 2).** For every tumor and every hallmark pathway we made a score and correlated it with SPINK1. We did this three times: unadjusted, adjusted for the myeloid score, and adjusted for the proliferation score (E2F targets plus G2M checkpoint genes).

Result:

- Unadjusted: 42 of 50 pathways are significant (FDR < 0.05).
- Myeloid adjusted: only 14 of 50 stay significant. Most of the lost ones are immune pathways, for example allograft rejection (ρ = 0.37 to -0.08) and IL-6/JAK/STAT3 signaling (0.43 to 0.11). These immune signals came from the number of immune cells in the tumor.
- Hypoxia stays: 0.54 unadjusted, 0.40 after myeloid adjustment (FDR < 10⁻⁴), 0.48 after proliferation adjustment.
- Proliferation adjusted: 31 of 48 pathways stay significant. E2F targets and G2M checkpoint are not tested here, because they are the covariate itself.
- Hypoxia is the strongest pathway before adjustment and after myeloid adjustment. After proliferation adjustment, glycolysis (0.49) and hypoxia (0.48) are together at the top. Glycolysis is a main part of the hypoxic response.

**Transcriptional subtypes (Table S7).** Glioblastoma has transcriptional subtypes. Mesenchymal tumors are known to have more hypoxia and more immune cells. So subtype could be a hidden reason for the SPINK1-hypoxia link. We gave each tumor the subtype with the highest signature score (proneural, classical or mesenchymal; marker genes from Verhaak et al. 2010 and Wang et al. 2017, listed in `spink1_utils.py`). The neural subtype is not used, because Wang et al. showed that it mostly reflects normal brain tissue in the sample.

To compare SPINK1 between subtypes we used the Kruskal-Wallis test. It compares more than two groups using ranks and does not need normal data.

Result: 45 proneural, 51 classical and 48 mesenchymal tumors. SPINK1 differs between subtypes (Kruskal-Wallis p = 8.0 × 10⁻⁶). It is highest in mesenchymal tumors and lowest in classical tumors (61% of classical tumors have no SPINK1). But inside every subtype SPINK1 is still positively correlated with all three hypoxia scores (9 of 9 tests positive). With the hallmark hypoxia score, ρ is 0.42 to 0.48 and p < 0.01 in each subtype. Some Buffa and Winter values do not reach p < 0.05 inside a subtype; each subtype has only 45 to 50 tumors.

**All covariates together, and a circularity problem (Tables S7, S8).** Then we adjusted for myeloid, mesenchymal and proliferation scores at the same time. Here there is a problem. Four genes of the mesenchymal signature (COL5A1, LOX, PLAUR, SERPINE1) are also in the hallmark hypoxia set. If we adjust for a score that contains hypoxia genes, we remove part of the hypoxia signal itself. So we also made a "clean" mesenchymal score without these four genes. Both versions are reported.

Result (n = 143, tumors with published hypoxia scores): the hallmark hypoxia link is ρ = 0.38 with the original mesenchymal score and 0.43 with the clean one. Buffa: 0.30 and 0.33. Winter: 0.32 and 0.36. All p ≤ 0.0003. So these three factors together do not explain the link. The clean score is used in all later notebooks.


### 05_hypoxia_validation: is the hypoxia score right? (Tables S9, S10, Figure 3)

Our hallmark hypoxia score is made by us. Maybe it measures something else. To check this, we used three hypoxia signatures that other groups made. Their scores come with the TCGA data and are available for 143 of the 144 tumors. Buffa and Winter are our main validation scores. Ragnum was made in prostate cancer, we include it for comparison.

**Result (Table S9).**

- Our hallmark score agrees with Buffa (ρ = 0.63, p = 1.9 × 10⁻¹⁷) and Winter (0.63, p = 7.5 × 10⁻¹⁷). So it captures largely the same signal as these two published scores.
- SPINK1 correlates with Buffa (ρ = 0.35, p = 1.7 × 10⁻⁵) and Winter (0.36, p = 8.5 × 10⁻⁶). After myeloid adjustment: 0.27 and 0.34.
- SPINK1 does not correlate with Ragnum (ρ = -0.03, p = 0.74). Ragnum was made in prostate cancer and it agrees only weakly with our hallmark score in these tumors (ρ = 0.21). A hypoxia signature made in one tumor type does not always work the same way in another.
- The median Buffa score goes up step by step from Undetectable to Low to High (3, 9 and 21; Kruskal-Wallis p = 2.3 × 10⁻⁵).
- The published scores use other gene lists and other calculation methods. So their ρ values are not directly comparable with the ρ of the hallmark score.

**Robustness (Table S10).** We tested the Buffa and Winter results in five ways, unadjusted and fully adjusted (myeloid, clean mesenchymal, proliferation):

- *Bootstrap confidence interval.* We draw 143 tumors with replacement and compute ρ again, 2,000 times. The middle 95% of these values is the 95% confidence interval. If it does not include 0, the result does not depend on which tumors happen to be in the sample. Result: all four intervals are above 0 (lowest bound 0.15, Buffa adjusted, 0.15 to 0.49).
- *Permutation test.* We shuffle SPINK1 values between tumors 2,000 times. This breaks any real link. The p-value is the share of shuffles that give a ρ as strong as the real one. Result: none of the 2,000 shuffles reached the real value (p < 0.0005).
- *Leave-one-out.* We remove one tumor at a time and compute ρ again, 143 times. This shows if a single tumor drives the result. Result: ρ changes by at most 0.03.
- *Only tumors with detectable SPINK1 (n = 88).* This checks if the result comes only from the difference between zero and non-zero tumors. Result: the link gets stronger, not weaker (adjusted Buffa 0.52, Winter 0.47). So the relation is graded, not only zero against non-zero.
- *Trimmed.* We remove the most extreme SPINK1 values (outside the 2.5 to 97.5 percentile range). Result: almost no change (ρ = 0.33 to 0.35).


### 06_cgga_replication: does it repeat in other patients? (Tables S12-S15, Figure 4)

**Cohorts.** We used the same rules as in TCGA: primary tumor, WHO grade IV, glioblastoma histology, IDH not mutant. CGGA-693: 116 tumors. CGGA-325: 74 tumors. No patient is in both cohorts. In CGGA-693, 7 of the 116 tumors have no IDH information (109 are IDH-wildtype). We kept them, in the same way as the 3 TCGA tumors without mutation data.

We analyze every cohort on its own and do not merge them. TCGA and CGGA give different data types (RSEM values and read counts). Keeping the cohorts apart gives three independent tests and needs no batch correction.

**Normalization.** CPM (counts per million) is the read count of a gene divided by all read counts of the sample, times one million. It corrects for different sequencing depth between samples. Then we used log2(CPM + 1). If a gene name appears twice, we keep the row with the higher mean count. A gene is kept if CPM > 1 in at least 20% of the tumors of the cohort (17,895 genes in CGGA-693, 18,004 in CGGA-325).

In CGGA-693, SPINK1 does not pass this filter because its level is low. But it is detected (count above 0) in 87 of 116 tumors (75%). It is the gene we study, so we kept it. Our tests use ranks, so a low level does not change the method. Random noise from low counts makes a correlation weaker, not stronger. In CGGA-325 SPINK1 passes the filter.

Scores are made in the same way as in TCGA, inside each cohort. Hallmark hypoxia: 178 and 177 of 200 genes found. Myeloid: 9 of 11. Clean mesenchymal: 10 of 10. Proliferation: 311 of 327. SPINK1 groups use the same rule as in TCGA (Undetectable/Low/High: 29/44/43 and 16/29/29).

**Replication (Table S12, Figure 4).**

- Unadjusted: ρ = 0.59 (p = 4.8 × 10⁻¹²) in CGGA-693 and 0.63 (p = 1.3 × 10⁻⁹) in CGGA-325. In TCGA it was 0.54.
- Fully adjusted (myeloid, clean mesenchymal, proliferation): 0.21 (p = 0.025) and 0.45 (p = 8.8 × 10⁻⁵). In TCGA 0.43. The adjusted value is smaller in CGGA-693, but it stays positive and significant.
- Only tumors with detectable SPINK1: 0.57 (n = 87) and 0.62 (n = 58), similar to the full cohorts.
- The 95% confidence intervals in Figure 4D use the Fisher z method with the variance for Spearman correlation from Bonett and Wright (2000): SE = sqrt(1.06 / (n - 3 - k)), where k is the number of covariates.

**Gene by gene (Table S13).** If the link is real, the same hypoxia genes should go with SPINK1 in every cohort. For each hallmark hypoxia gene we computed its correlation with SPINK1 in each cohort. 176 genes have data in all three. Then we compared these 176 values between cohorts.

Result: the gene profiles agree well (ρ = 0.65 between TCGA and CGGA-693, 0.86 between TCGA and CGGA-325, 0.73 between the two CGGA cohorts). 130 of the 176 genes are positive in all three cohorts. Of the 89 genes with FDR < 0.05 in TCGA, 81 (91%) are positive in all three. The most consistent genes are known HIF targets: NDRG1, LOX, HMOX1, ANGPTL4 and PLIN2.

**Pipeline controls (Table S14).** Before trusting a new result, we check that the same pipeline finds what we already know. A positive control gene should correlate clearly with its own score. A negative control gene should not. Each gene is removed from the score before the test.

- CA9, the classic HIF target, correlates with the hypoxia score (0.68 and 0.75). CA9 is not in the hallmark hypoxia set, so this is an outside check.
- PTPRC (CD45) with the myeloid score: 0.67 and 0.78. (CD68 is not in the CGGA gene list, so PTPRC was used.)
- MKI67 and TOP2A with the proliferation score: 0.84 to 0.93.
- CHI3L1 (YKL-40) with the mesenchymal score: 0.73 and 0.62.
- Neuronal genes (SNAP25, SYT1, RBFOX3) with the hypoxia score: close to zero and not significant in CGGA-693 (0.06 to 0.16), negative in CGGA-325 (-0.26 to -0.40).
- SPINK1 itself, for comparison: 0.59 and 0.63.

**Random gene sets (Table S15).** In bulk tissue almost any large gene set shares some general signal (for example the share of tumor cells in the sample). So maybe SPINK1 correlates with any big gene set, not only with hypoxia. To test this we drew 1,000 random gene sets of the same size as the hypoxia set (without hypoxia genes), made a score from each, and correlated it with SPINK1. This shows what we can get by chance.

- z-score = (real ρ - mean of random ρ) / SD of random ρ. It says how many standard deviations the real value is above chance.
- Empirical p = share of random sets with |ρ| at least as large as the real |ρ|.

Result: random sets also correlate a little with SPINK1 (mean ρ = 0.24 and 0.11), as expected. But the hypoxia set is much higher: 10.0 and 8.0 SD above the random mean (empirical p < 0.001 in both cohorts). CA9 gives a similar picture (z = 7.9 and 10.6). The negative control SNAP25 is far below the random mean (z = -5.8 and -4.3). In CGGA-693 its empirical p is 1.000: every random set correlates more strongly with SNAP25 than the hypoxia set does. This is the expected result for a negative control.


### 07_survival: overall survival (Table S11, Figure S1)

Parts b and c of Table S11 and Figure S1 were added during peer review.

**Some survival terms.**

- *Kaplan-Meier curve.* The share of patients still alive over time. Every step down is a death. The small vertical marks are censored patients: alive at the last follow-up, so their full survival time is not known.
- *Log-rank test.* Compares whole survival curves between groups.
- *Cox model, hazard ratio (HR).* HR above 1 means higher risk of death, below 1 lower risk, 1 means no effect. If the 95% confidence interval includes 1, the effect is not significant. A Cox model can include other factors, like age and sex, at the same time.
- *HR per 1 SD.* SPINK1 has different units in TCGA (RSEM) and CGGA (CPM). So we divided SPINK1 by its standard deviation in each cohort. Then the HR is the change in risk for one SD more SPINK1, and the cohorts can be compared.
- *Proportional hazards test.* The Cox model assumes that the effect of SPINK1 is the same over the whole follow-up time. We tested this with Schoenfeld residuals. p > 0.05 means no sign that the assumption is broken.

**Table S11a: the TCGA analysis of the first submission.** In TCGA, age and sex are missing for 45 of the 143 patients with follow-up data. To compare models with and without covariates on the same patients, this part uses the 98 patients with complete data. HR per log2 unit of SPINK1: 1.00 (95% CI 0.87 to 1.15, p = 0.98). Adding age and sex, then the Buffa hypoxia score, then the myeloid, mesenchymal and proliferation scores does not change this (HR 0.95 to 1.02, all p > 0.6).

**Table S11b, S11c and Figure S1: three cohorts.**

- Cox, unadjusted: TCGA HR 1.10 (0.92 to 1.31), CGGA-693 0.89 (0.71 to 1.12), CGGA-325 1.01 (0.80 to 1.26). With age and sex: 1.04, 0.86 and 1.02. All p > 0.2. The HRs point in different directions and all intervals include 1.
- Proportional hazards tests: all p > 0.05 (lowest 0.096).
- Kaplan-Meier with the three SPINK1 groups: log-rank p = 0.40 (TCGA), 0.50 (CGGA-693) and 0.41 (CGGA-325). Median survival in Undetectable/Low/High: 14.9/11.7/10.6 months in TCGA, 12.0/13.3/16.2 in CGGA-693, 13.9/10.3/15.9 in CGGA-325. The order of the groups is not the same in the three cohorts.
- Best cut-off scan. We also tried every cut-off that leaves at least 20% of the patients in each group (59, 62 and 44 cut-offs). This is close to the method of the survival plots in the Human Protein Atlas. Taking the smallest p of many tests is too optimistic, because some small p-values come by chance. Even so, the smallest p was 0.097, 0.148 and 0.103. No cut-off reached 0.05.
- All primary gliomas in CGGA (grade 2 to 4, both cohorts together, IDH status known, n = 588). Here more SPINK1 goes with worse survival (HR 1.26, 1.16 to 1.38, p = 2.4 × 10⁻⁷). But when grade and IDH status are added to the model, the effect is gone (HR 0.96, 0.87 to 1.06, p = 0.41). So in a mixed glioma group the survival difference is explained by grade and IDH status, not by SPINK1 itself.

In short: inside IDH-wildtype glioblastoma, we see no link between SPINK1 and overall survival, in any of the three cohorts and with any of the methods above.


## Why some numbers of patients differ

- **144**: the TCGA cohort.
- **143** (notebooks 04 and 05): one tumor (TCGA-14-0871) has no published hypoxia scores. All analyses that use Buffa, Winter or Ragnum use 143 tumors.
- **143** (notebook 07): one other tumor (TCGA-28-2510) has no follow-up data.
- **98** (Table S11a): patients with follow-up, age and sex.
- **140** (Table S2): tumors with copy number data.
- **116 and 74**: the CGGA cohorts. In survival analysis CGGA-693 has 112 patients, because 4 have no follow-up data.


## Citation

If you use this code, please cite the article. The full reference will be added after publication.


## Contact

Code: Melisa Ağrı (agrim22@itu.edu.tr)
Corresponding author: Zeliha Emrence (zeliha.emrence@istanbul.edu.tr)
