# HYDRA-ccRCC High-Confidence Gene Dossiers

These dossiers are generated from the reproducible pipeline. They are interpretation scaffolds, not final biological claims.
Candidates are listed alphabetically. The discovery evidence score is a selection heuristic, not an external-validation or biological rank.

## ACADM

- Gene name: acyl-CoA dehydrogenase medium chain
- Pathway class: Metabolism
- Tumor-normal signal: TCGA log2FC -1.545, GSE40435 log2FC -1.5, GSE53757 log2FC -1.637
- Survival signal: HR 0.639 (95% CI 0.55-0.741), FDR 9.38e-07, PH p 0.239
- Cell-type sanity: Mitochondrial fatty-acid metabolism; plausible renal epithelial metabolism.
- Literature prior: Metabolic cancer relevance likely, ccRCC specificity needs review.
- Manual review note: Assess metabolic adaptation against normal-kidney retention and external outcomes.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=ACADM%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## C1QTNF6

- Gene name: C1q and TNF related 6
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC 2.327, GSE40435 log2FC 0.277, GSE53757 log2FC 1.222
- Survival signal: HR 1.545 (95% CI 1.321-1.807), FDR 4.42e-06, PH p 0.126
- Cell-type sanity: Secreted/metabolic-inflammatory signal; cell type unresolved.
- Literature prior: Cancer literature likely; ccRCC-specific role needs review.
- Manual review note: Check immune, adipokine, and tumor-cell sources.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=C1QTNF6%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## CADPS2

- Gene name: calcium dependent secretion activator 2
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.628, GSE40435 log2FC -0.452, GSE53757 log2FC -1.286
- Survival signal: HR 0.664 (95% CI 0.572-0.77), FDR 4.93e-06, PH p 0.478
- Cell-type sanity: Vesicle/secretion-related signal; not kidney-cell-type specific here.
- Literature prior: ccRCC literature requires manual confirmation.
- Manual review note: Lower-confidence biology unless supported externally.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CADPS2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## CLCN5

- Gene name: Cl-/H+ antiporter 5
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.434, GSE40435 log2FC -0.752, GSE53757 log2FC -1.8
- Survival signal: HR 0.659 (95% CI 0.566-0.766), FDR 4.71e-06, PH p 0.825
- Cell-type sanity: Renal proximal tubule/endocytic biology prior.
- Literature prior: Strong kidney biology prior; ccRCC prognostic agreement needs review.
- Manual review note: May mark retained normal tubule differentiation rather than aggressive tumor biology.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CLCN5%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## CRYL1

- Gene name: crystallin lambda 1
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.066, GSE40435 log2FC -1.696, GSE53757 log2FC -1.718
- Survival signal: HR 0.634 (95% CI 0.555-0.723), FDR 1.98e-08, PH p 0.181
- Cell-type sanity: Metabolic enzyme signal compatible with renal epithelial biology.
- Literature prior: Kidney cancer literature requires manual confirmation.
- Manual review note: Check whether downregulation is a renal identity/metabolic-loss signal.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CRYL1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## CYFIP2

- Gene name: cytoplasmic FMR1 interacting protein 2
- Pathway class: Metabolism
- Tumor-normal signal: TCGA log2FC -2.281, GSE40435 log2FC -1.132, GSE53757 log2FC -1.686
- Survival signal: HR 0.644 (95% CI 0.568-0.73), FDR 1.65e-08, PH p 0.11
- Cell-type sanity: Broad cytoskeletal/regulatory signal; not kidney-cell-type specific from this pipeline alone.
- Literature prior: ccRCC-specific role requires manual confirmation.
- Manual review note: Prioritize whether the survival signal reflects tumor biology or tissue composition.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CYFIP2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## DDC

- Gene name: dopa decarboxylase
- Pathway class: Metabolism
- Tumor-normal signal: TCGA log2FC -1.7, GSE40435 log2FC -3.084, GSE53757 log2FC -3.226
- Survival signal: HR 0.659 (95% CI 0.573-0.758), FDR 9.45e-07, PH p 0.157
- Cell-type sanity: Metabolic/decarboxylase signal compatible with epithelial metabolic rewiring.
- Literature prior: Kidney cancer literature requires manual confirmation.
- Manual review note: Check whether expression is kidney-lineage, neuroendocrine-like, or stromal.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=DDC%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## GJB1

- Gene name: gap junction protein beta 1
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.467, GSE40435 log2FC -1.369, GSE53757 log2FC -1.853
- Survival signal: HR 0.652 (95% CI 0.567-0.751), FDR 7.55e-07, PH p 0.057
- Cell-type sanity: Gap-junction protein with renal epithelial expression but stronger expression in other normal tissues.
- Literature prior: Connexin biology is established; ccRCC-specific prognostic evidence needs review.
- Manual review note: Treat as a compartment-sensitive epithelial hypothesis and report the non-proportional-hazards diagnostic.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=GJB1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## GRAMD1A

- Gene name: GRAM domain containing 1A
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC 1.317, GSE40435 log2FC 0.914, GSE53757 log2FC 0.572
- Survival signal: HR 1.5 (95% CI 1.273-1.769), FDR 3.77e-05, PH p 0.551
- Cell-type sanity: Cholesterol transport/contact-site biology; cell type unresolved.
- Literature prior: ccRCC literature requires manual confirmation.
- Manual review note: Review with lipid/metabolic adaptation lens.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=GRAMD1A%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## HHLA2

- Gene name: HHLA2 member of B7 family
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC 3.119, GSE40435 log2FC 0.726, GSE53757 log2FC 1.189
- Survival signal: HR 0.627 (95% CI 0.547-0.719), FDR 3.06e-08, PH p 0.056
- Cell-type sanity: Immune-checkpoint-like signal; may reflect tumor-immune interaction.
- Literature prior: Known cancer immunology gene; ccRCC-specific direction needs review.
- Manual review note: Interpret as immune/microenvironmental until cell source is resolved.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=HHLA2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## HIBCH

- Gene name: 3-hydroxyisobutyryl-CoA hydrolase
- Pathway class: Metabolism
- Tumor-normal signal: TCGA log2FC -1.285, GSE40435 log2FC -1.918, GSE53757 log2FC -1.758
- Survival signal: HR 0.664 (95% CI 0.56-0.786), FDR 5.25e-05, PH p 0.305
- Cell-type sanity: Mitochondrial branched-chain amino-acid metabolism.
- Literature prior: Metabolic cancer relevance likely, ccRCC specificity needs review.
- Manual review note: Candidate metabolic-adaptation gene.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=HIBCH%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## IFFO1

- Gene name: intermediate filament family orphan 1
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC 1.779, GSE40435 log2FC 0.416, GSE53757 log2FC 0.886
- Survival signal: HR 1.505 (95% CI 1.273-1.778), FDR 4.32e-05, PH p 0.279
- Cell-type sanity: Poorly resolved structural/nuclear signal.
- Literature prior: ccRCC literature likely sparse.
- Manual review note: Keep as computational association unless external evidence appears.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=IFFO1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## KL

- Gene name: klotho
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.639, GSE40435 log2FC -2.45, GSE53757 log2FC -2.426
- Survival signal: HR 0.661 (95% CI 0.576-0.758), FDR 7.83e-07, PH p 0.195
- Cell-type sanity: Kidney-enriched aging/mineral-axis gene; plausible renal epithelial context.
- Literature prior: Strong kidney biology prior; ccRCC prognostic agreement needs review.
- Manual review note: A lower tumor expression signal may reflect loss of normal kidney program.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=KL%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## KNTC1

- Gene name: kinetochore associated 1
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC 1.156, GSE40435 log2FC 0.61, GSE53757 log2FC 1.113
- Survival signal: HR 1.512 (95% CI 1.304-1.752), FDR 4.08e-06, PH p 0.647
- Cell-type sanity: NA
- Literature prior: NA
- Manual review note: NA
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=KNTC1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## LRBA

- Gene name: LPS responsive beige-like anchor protein
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.251, GSE40435 log2FC -0.334, GSE53757 log2FC -0.716
- Survival signal: HR 0.655 (95% CI 0.557-0.77), FDR 1.3e-05, PH p 0.372
- Cell-type sanity: Immune-regulatory signal.
- Literature prior: Immune literature strong; ccRCC-specific direction needs review.
- Manual review note: Likely immune-composition sensitivity.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=LRBA%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## LTB4R

- Gene name: leukotriene B4 receptor
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC 2.586, GSE40435 log2FC 0.298, GSE53757 log2FC 1.038
- Survival signal: HR 1.515 (95% CI 1.285-1.785), FDR 2.36e-05, PH p 0.029
- Cell-type sanity: Leukotriene-receptor inflammatory signal with renal epithelial and immune expression.
- Literature prior: Inflammatory signaling is biologically plausible; ccRCC-specific cell source needs review.
- Manual review note: Treat as a composition-sensitive risk association and report the non-proportional-hazards diagnostic.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=LTB4R%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## PANK1

- Gene name: pantothenate kinase 1
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.487, GSE40435 log2FC -0.974, GSE53757 log2FC -1.539
- Survival signal: HR 0.623 (95% CI 0.535-0.726), FDR 4.45e-07, PH p 0.301
- Cell-type sanity: Metabolism/coenzyme A pathway signal; plausible metabolic adaptation.
- Literature prior: Kidney cancer literature requires manual confirmation.
- Manual review note: Review whether signal is tumor-intrinsic metabolism or normal-tissue retention.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=PANK1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## PODXL

- Gene name: podocalyxin like
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.369, GSE40435 log2FC -0.485, GSE53757 log2FC -1.256
- Survival signal: HR 0.662 (95% CI 0.567-0.772), FDR 8.77e-06, PH p 0.399
- Cell-type sanity: Podocyte/endothelial/cell-adhesion signal.
- Literature prior: Cancer invasion literature likely; kidney compartment matters.
- Manual review note: Cell-type sanity flag: glomerular/endothelial biology may confound bulk signal.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=PODXL%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## RBM47

- Gene name: RNA binding motif protein 47
- Pathway class: Unclassified
- Tumor-normal signal: TCGA log2FC -1.103, GSE40435 log2FC -0.979, GSE53757 log2FC -0.965
- Survival signal: HR 0.666 (95% CI 0.588-0.755), FDR 1.21e-07, PH p 0.04
- Cell-type sanity: RNA-binding and splicing-regulatory signal with broad normal-tissue expression.
- Literature prior: Cancer-regulatory literature exists; ccRCC-specific evidence needs review.
- Manual review note: Interpret cautiously because HPA mapping is immune-dominant and the TCGA proportional-hazards diagnostic is nominally significant.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=RBM47%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## TCIRG1

- Gene name: T cell immune regulator 1, ATPase H+ transporting V0 subunit a3
- Pathway class: Metabolism
- Tumor-normal signal: TCGA log2FC 1.455, GSE40435 log2FC 0.413, GSE53757 log2FC 1.719
- Survival signal: HR 1.531 (95% CI 1.297-1.807), FDR 1.9e-05, PH p 0.459
- Cell-type sanity: Immune/lysosomal/proton-pump-associated signal; cell type unresolved.
- Literature prior: Immune/cancer literature likely, ccRCC specificity needs review.
- Manual review note: Check immune-cell composition before tumor-cell interpretation.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TCIRG1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## TMEM140

- Gene name: transmembrane protein 140
- Pathway class: Immune
- Tumor-normal signal: TCGA log2FC 1.387, GSE40435 log2FC 1.043, GSE53757 log2FC 1.336
- Survival signal: HR 0.665 (95% CI 0.575-0.77), FDR 4.42e-06, PH p 0.536
- Cell-type sanity: NA
- Literature prior: NA
- Manual review note: NA
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TMEM140%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## TNFAIP2

- Gene name: TNF alpha induced protein 2
- Pathway class: Immune
- Tumor-normal signal: TCGA log2FC 1.29, GSE40435 log2FC 0.352, GSE53757 log2FC 0.437
- Survival signal: HR 1.533 (95% CI 1.271-1.849), FDR 1.34e-04, PH p 0.748
- Cell-type sanity: TNF/inflammatory-response signal.
- Literature prior: Inflammation/cancer literature likely; ccRCC specificity needs review.
- Manual review note: Interpret as immune/inflammatory until cell source is resolved.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TNFAIP2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## TNIP1

- Gene name: TNFAIP3 interacting protein 1
- Pathway class: Immune
- Tumor-normal signal: TCGA log2FC 1.047, GSE40435 log2FC 0.674, GSE53757 log2FC 0.293
- Survival signal: HR 0.66 (95% CI 0.567-0.769), FDR 6.17e-06, PH p 0.196
- Cell-type sanity: NA
- Literature prior: NA
- Manual review note: NA
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TNIP1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

