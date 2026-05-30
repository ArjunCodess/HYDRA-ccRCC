# HYDRA-ccRCC High-Confidence Gene Dossiers

These dossiers are generated from the reproducible pipeline. They are interpretation scaffolds, not final biological claims.

## 1. DDC

- Gene name: dopa decarboxylase
- Pathway class: Metabolism
- Final rank score: 0.667
- Tumor-normal signal: TCGA log2FC -1.712, GSE40435 log2FC -3.084, GSE53757 log2FC -3.226
- Survival signal: HR 0.659 (95% CI 0.572-0.76), FDR 1.59e-06, PH p 0.14
- Cell-type sanity: Metabolic/decarboxylase signal compatible with epithelial metabolic rewiring.
- Literature prior: Kidney cancer literature requires manual confirmation.
- Manual review note: Check whether expression is kidney-lineage, neuroendocrine-like, or stromal.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=DDC%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 2. CYFIP2

- Gene name: cytoplasmic FMR1 interacting protein 2
- Pathway class: Metabolism
- Final rank score: 0.656
- Tumor-normal signal: TCGA log2FC -2.278, GSE40435 log2FC -1.132, GSE53757 log2FC -1.686
- Survival signal: HR 0.643 (95% CI 0.568-0.729), FDR 1.55e-08, PH p 0.114
- Cell-type sanity: Broad cytoskeletal/regulatory signal; not kidney-cell-type specific from this pipeline alone.
- Literature prior: ccRCC-specific role requires manual confirmation.
- Manual review note: Prioritize whether the survival signal reflects tumor biology or tissue composition.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CYFIP2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 3. KL

- Gene name: klotho
- Pathway class: Unclassified
- Final rank score: 0.579
- Tumor-normal signal: TCGA log2FC -1.643, GSE40435 log2FC -2.45, GSE53757 log2FC -2.426
- Survival signal: HR 0.661 (95% CI 0.576-0.759), FDR 9.35e-07, PH p 0.191
- Cell-type sanity: Kidney-enriched aging/mineral-axis gene; plausible renal epithelial context.
- Literature prior: Strong kidney biology prior; ccRCC prognostic agreement needs review.
- Manual review note: A lower tumor expression signal may reflect loss of normal kidney program.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=KL%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 4. ACADM

- Gene name: acyl-CoA dehydrogenase medium chain
- Pathway class: Metabolism
- Final rank score: 0.572
- Tumor-normal signal: TCGA log2FC -1.542, GSE40435 log2FC -1.5, GSE53757 log2FC -1.637
- Survival signal: HR 0.638 (95% CI 0.549-0.741), FDR 8.63e-07, PH p 0.246
- Cell-type sanity: Mitochondrial fatty-acid metabolism; plausible renal epithelial metabolism.
- Literature prior: Metabolic cancer relevance likely, ccRCC specificity needs review.
- Manual review note: Strong candidate for metabolic-adaptation interpretation.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=ACADM%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 5. CRYL1

- Gene name: crystallin lambda 1
- Pathway class: Unclassified
- Final rank score: 0.526
- Tumor-normal signal: TCGA log2FC -1.075, GSE40435 log2FC -1.696, GSE53757 log2FC -1.718
- Survival signal: HR 0.636 (95% CI 0.556-0.726), FDR 4.59e-08, PH p 0.169
- Cell-type sanity: Metabolic enzyme signal compatible with renal epithelial biology.
- Literature prior: Kidney cancer literature requires manual confirmation.
- Manual review note: Check whether downregulation is a renal identity/metabolic-loss signal.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CRYL1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 6. CLCN5

- Gene name: Cl-/H+ antiporter 5
- Pathway class: Unclassified
- Final rank score: 0.519
- Tumor-normal signal: TCGA log2FC -1.433, GSE40435 log2FC -0.752, GSE53757 log2FC -1.8
- Survival signal: HR 0.658 (95% CI 0.567-0.765), FDR 4.37e-06, PH p 0.82
- Cell-type sanity: Renal proximal tubule/endocytic biology prior.
- Literature prior: Strong kidney biology prior; ccRCC prognostic agreement needs review.
- Manual review note: May mark retained normal tubule differentiation rather than aggressive tumor biology.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CLCN5%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 7. HIBCH

- Gene name: 3-hydroxyisobutyryl-CoA hydrolase
- Pathway class: Metabolism
- Final rank score: 0.479
- Tumor-normal signal: TCGA log2FC -1.275, GSE40435 log2FC -1.918, GSE53757 log2FC -1.758
- Survival signal: HR 0.661 (95% CI 0.557-0.783), FDR 4.29e-05, PH p 0.32
- Cell-type sanity: Mitochondrial branched-chain amino-acid metabolism.
- Literature prior: Metabolic cancer relevance likely, ccRCC specificity needs review.
- Manual review note: Candidate metabolic-adaptation gene.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=HIBCH%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 8. HHLA2

- Gene name: HHLA2 member of B7 family
- Pathway class: Unclassified
- Final rank score: 0.479
- Tumor-normal signal: TCGA log2FC 3.107, GSE40435 log2FC 0.726, GSE53757 log2FC 1.189
- Survival signal: HR 0.626 (95% CI 0.544-0.72), FDR 6.33e-08, PH p 0.054
- Cell-type sanity: Immune-checkpoint-like signal; may reflect tumor-immune interaction.
- Literature prior: Known cancer immunology gene; ccRCC-specific direction needs review.
- Manual review note: Interpret as immune/microenvironmental until cell source is resolved.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=HHLA2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 9. TCIRG1

- Gene name: T cell immune regulator 1, ATPase H+ transporting V0 subunit a3
- Pathway class: Metabolism
- Final rank score: 0.466
- Tumor-normal signal: TCGA log2FC 1.448, GSE40435 log2FC 0.413, GSE53757 log2FC 1.719
- Survival signal: HR 1.536 (95% CI 1.301-1.814), FDR 1.51e-05, PH p 0.478
- Cell-type sanity: Immune/lysosomal/proton-pump-associated signal; cell type unresolved.
- Literature prior: Immune/cancer literature likely, ccRCC specificity needs review.
- Manual review note: Check immune-cell composition before tumor-cell interpretation.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TCIRG1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 10. PANK1

- Gene name: pantothenate kinase 1
- Pathway class: Unclassified
- Final rank score: 0.452
- Tumor-normal signal: TCGA log2FC -1.48, GSE40435 log2FC -0.974, GSE53757 log2FC -1.539
- Survival signal: HR 0.621 (95% CI 0.532-0.724), FDR 3.92e-07, PH p 0.318
- Cell-type sanity: Metabolism/coenzyme A pathway signal; plausible metabolic adaptation.
- Literature prior: Kidney cancer literature requires manual confirmation.
- Manual review note: Review whether signal is tumor-intrinsic metabolism or normal-tissue retention.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=PANK1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 11. ACAT1

- Gene name: acetyl-CoA acetyltransferase 1
- Pathway class: Metabolism
- Final rank score: 0.42
- Tumor-normal signal: TCGA log2FC -1.194, GSE40435 log2FC -1.829, GSE53757 log2FC -0.839
- Survival signal: HR 0.664 (95% CI 0.559-0.789), FDR 7.18e-05, PH p 0.4
- Cell-type sanity: Mitochondrial acetyl-CoA/ketone and amino-acid metabolism.
- Literature prior: Metabolic cancer relevance likely, ccRCC specificity needs review.
- Manual review note: Strong candidate for metabolic-adaptation interpretation.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=ACAT1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 12. FUT6

- Gene name: fucosyltransferase 6
- Pathway class: Unclassified
- Final rank score: 0.412
- Tumor-normal signal: TCGA log2FC -1.31, GSE40435 log2FC -1.125, GSE53757 log2FC -2.368
- Survival signal: HR 0.666 (95% CI 0.574-0.773), FDR 6.13e-06, PH p 0.218
- Cell-type sanity: Glycosylation signal; cell type unresolved.
- Literature prior: Cancer glycosylation literature likely, ccRCC specificity needs review.
- Manual review note: Treat as pathway hypothesis until kidney/cancer literature is checked.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=FUT6%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 13. TEK

- Gene name: TEK receptor tyrosine kinase
- Pathway class: Unclassified
- Final rank score: 0.399
- Tumor-normal signal: TCGA log2FC -1.438, GSE40435 log2FC -0.67, GSE53757 log2FC -0.858
- Survival signal: HR 0.666 (95% CI 0.577-0.77), FDR 3.62e-06, PH p 0.464
- Cell-type sanity: Endothelial/angiogenesis-associated signal.
- Literature prior: Angiogenesis prior is relevant to ccRCC, but cell source is endothelial.
- Manual review note: Interpret as vascular composition or angiogenic state, not necessarily tumor-cell intrinsic.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TEK%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 14. IQGAP2

- Gene name: IQ motif containing GTPase activating protein 2
- Pathway class: Unclassified
- Final rank score: 0.373
- Tumor-normal signal: TCGA log2FC -1.106, GSE40435 log2FC -1.102, GSE53757 log2FC -0.712
- Survival signal: HR 0.66 (95% CI 0.579-0.753), FDR 2.21e-07, PH p 0.104
- Cell-type sanity: Scaffold/cytoskeletal signaling; not cell-type-resolved here.
- Literature prior: Cancer literature exists; ccRCC-specific evidence needs review.
- Manual review note: Check direction against tumor suppressor versus composition explanations.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=IQGAP2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 15. TNFAIP2

- Gene name: TNF alpha induced protein 2
- Pathway class: Immune
- Final rank score: 0.373
- Tumor-normal signal: TCGA log2FC 1.284, GSE40435 log2FC 0.352, GSE53757 log2FC 0.437
- Survival signal: HR 1.538 (95% CI 1.275-1.855), FDR 1.22e-04, PH p 0.75
- Cell-type sanity: TNF/inflammatory-response signal.
- Literature prior: Inflammation/cancer literature likely; ccRCC specificity needs review.
- Manual review note: Interpret as immune/inflammatory until cell source is resolved.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=TNFAIP2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 16. CADPS2

- Gene name: calcium dependent secretion activator 2
- Pathway class: Unclassified
- Final rank score: 0.369
- Tumor-normal signal: TCGA log2FC -1.6, GSE40435 log2FC -0.452, GSE53757 log2FC -1.286
- Survival signal: HR 0.654 (95% CI 0.562-0.762), FDR 4.29e-06, PH p 0.502
- Cell-type sanity: Vesicle/secretion-related signal; not kidney-cell-type specific here.
- Literature prior: ccRCC literature requires manual confirmation.
- Manual review note: Lower-confidence biology unless supported externally.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=CADPS2%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 17. EMCN

- Gene name: endomucin
- Pathway class: Unclassified
- Final rank score: 0.359
- Tumor-normal signal: TCGA log2FC -1.124, GSE40435 log2FC -0.633, GSE53757 log2FC -1.316
- Survival signal: HR 0.666 (95% CI 0.574-0.773), FDR 6.38e-06, PH p 0.42
- Cell-type sanity: Endothelial marker-like signal.
- Literature prior: Kidney vascular literature likely; ccRCC survival interpretation needs review.
- Manual review note: Cell-type sanity flag: likely vascular composition.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=EMCN%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 18. C1QTNF6

- Gene name: C1q and TNF related 6
- Pathway class: Unclassified
- Final rank score: 0.339
- Tumor-normal signal: TCGA log2FC 2.315, GSE40435 log2FC 0.277, GSE53757 log2FC 1.222
- Survival signal: HR 1.566 (95% CI 1.335-1.836), FDR 3.62e-06, PH p 0.131
- Cell-type sanity: Secreted/metabolic-inflammatory signal; cell type unresolved.
- Literature prior: Cancer literature likely; ccRCC-specific role needs review.
- Manual review note: Check immune, adipokine, and tumor-cell sources.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=C1QTNF6%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 19. PODXL

- Gene name: podocalyxin like
- Pathway class: Unclassified
- Final rank score: 0.335
- Tumor-normal signal: TCGA log2FC -1.372, GSE40435 log2FC -0.485, GSE53757 log2FC -1.256
- Survival signal: HR 0.664 (95% CI 0.569-0.775), FDR 1.08e-05, PH p 0.4
- Cell-type sanity: Podocyte/endothelial/cell-adhesion signal.
- Literature prior: Cancer invasion literature likely; kidney compartment matters.
- Manual review note: Cell-type sanity flag: glomerular/endothelial biology may confound bulk signal.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=PODXL%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 20. DBT

- Gene name: dihydrolipoamide branched chain transacylase E2
- Pathway class: Metabolism
- Final rank score: 0.325
- Tumor-normal signal: TCGA log2FC -1.558, GSE40435 log2FC -0.53, GSE53757 log2FC -0.705
- Survival signal: HR 0.659 (95% CI 0.552-0.787), FDR 8.64e-05, PH p 0.284
- Cell-type sanity: Branched-chain amino-acid metabolism.
- Literature prior: Metabolic cancer relevance likely, ccRCC specificity needs review.
- Manual review note: Strong candidate for metabolic-adaptation interpretation.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=DBT%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 21. LRBA

- Gene name: LPS responsive beige-like anchor protein
- Pathway class: Unclassified
- Final rank score: 0.322
- Tumor-normal signal: TCGA log2FC -1.199, GSE40435 log2FC -0.334, GSE53757 log2FC -0.716
- Survival signal: HR 0.636 (95% CI 0.536-0.753), FDR 9.71e-06, PH p 0.393
- Cell-type sanity: Immune-regulatory signal.
- Literature prior: Immune literature strong; ccRCC-specific direction needs review.
- Manual review note: Likely immune-composition sensitivity.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=LRBA%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 22. GRAMD1A

- Gene name: GRAM domain containing 1A
- Pathway class: Unclassified
- Final rank score: 0.3
- Tumor-normal signal: TCGA log2FC 1.306, GSE40435 log2FC 0.914, GSE53757 log2FC 0.572
- Survival signal: HR 1.55 (95% CI 1.302-1.845), FDR 2.6e-05, PH p 0.564
- Cell-type sanity: Cholesterol transport/contact-site biology; cell type unresolved.
- Literature prior: ccRCC literature requires manual confirmation.
- Manual review note: Review with lipid/metabolic adaptation lens.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=GRAMD1A%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 23. IFFO1

- Gene name: intermediate filament family orphan 1
- Pathway class: Unclassified
- Final rank score: 0.268
- Tumor-normal signal: TCGA log2FC 1.769, GSE40435 log2FC 0.416, GSE53757 log2FC 0.886
- Survival signal: HR 1.523 (95% CI 1.285-1.804), FDR 3.32e-05, PH p 0.289
- Cell-type sanity: Poorly resolved structural/nuclear signal.
- Literature prior: ccRCC literature likely sparse.
- Manual review note: Keep as computational association unless external evidence appears.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=IFFO1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

## 24. FHOD1

- Gene name: formin homology 2 domain containing 1
- Pathway class: Unclassified
- Final rank score: 0.149
- Tumor-normal signal: TCGA log2FC 1.29, GSE40435 log2FC 0.506, GSE53757 log2FC 0.536
- Survival signal: HR 1.517 (95% CI 1.26-1.826), FDR 1.7e-04, PH p 0.105
- Cell-type sanity: Actin/cytoskeletal remodeling signal.
- Literature prior: Cancer migration literature likely; ccRCC specificity needs review.
- Manual review note: Check whether survival signal reflects invasion biology or stromal composition.
- PubMed query: https://pubmed.ncbi.nlm.nih.gov/?term=FHOD1%20%28ccRCC%20OR%20clear%20cell%20renal%20cell%20carcinoma%20OR%20kidney%20cancer%29

