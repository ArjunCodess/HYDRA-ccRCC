# HYDRA-ccRCC Manual Biological Review of 24 High-Confidence Candidate Genes

Date: 2026-05-31

Scope: This is a literature-anchored biological review of the 24 final candidate prognostic associations. It does not add new analyses. The review uses the existing HYDRA-ccRCC outputs, especially `high_confidence_ranked_shortlist.csv`, and evaluates whether each association is biologically plausible enough to discuss in a manuscript.

Interpretation of direction:

- `Tumor-normal`: positive log2FC means higher in ccRCC tumor than adjacent normal kidney; negative means lower in tumor.
- `Survival`: HR < 1 means higher expression is associated with lower hazard, after stage/grade/age/sex adjustment in the existing pipeline; HR > 1 means higher expression is associated with higher hazard.
- Bulk RNA-seq cannot distinguish tumor-cell-intrinsic expression from renal epithelial retention, endothelial content, stromal admixture, necrosis, or immune composition.

## Executive Ranking After Biological Review

| Review rank | Gene | Pipeline direction | Evidence quality | Main interpretation | Confidence |
|---:|---|---|---|---|---:|
| 1 | KL | Down in tumor; higher expression protective | Strong | Renal epithelial/tumor-suppressive marker; likely prognostic but also normal-kidney retention | 88 |
| 2 | ACADM | Down in tumor; higher expression protective | Strong | Fatty-acid oxidation/mitochondrial metabolic retention | 84 |
| 3 | CRYL1 | Down in tumor; higher expression protective | Strong | Renal metabolic/tumor-suppressive marker | 84 |
| 4 | ACAT1 | Down in tumor; higher expression protective | Strong | Ketone/lipid metabolism tumor suppressor-like signal | 83 |
| 5 | DDC | Down in tumor; higher expression protective | Strong | Amino-acid decarboxylase; validated ccRCC prognostic marker, mechanism less canonical | 82 |
| 6 | PANK1 | Down in tumor; higher expression protective | Strong | CoA metabolism and mitochondrial function; plausible metabolic marker | 78 |
| 7 | TCIRG1 | Up in tumor; higher expression risky | Strong | V-ATPase/lysosomal/glycolytic-immune risk marker; immune confounding likely | 78 |
| 8 | DBT | Down in tumor; higher expression protective | Strong | Branched-chain amino-acid metabolism; functional ccRCC support | 77 |
| 9 | CLCN5 | Down in tumor; higher expression protective | Moderate-strong | Proximal tubule/endosomal transporter; likely retained renal differentiation/metabolic restraint | 75 |
| 10 | CYFIP2 | Down in tumor; higher expression protective | Moderate | p53/cytoskeleton/Rho-GTPase-associated tumor suppressor-like signal | 73 |
| 11 | C1QTNF6 | Up in tumor; higher expression risky | Moderate | Secreted inflammatory/adipokine-like risk marker | 70 |
| 12 | IQGAP2 | Down in tumor; higher expression protective | Moderate | Cytoskeletal scaffold; tumor-suppressor-like but not ccRCC-specific enough | 68 |
| 13 | TNFAIP2 | Up in tumor; higher expression risky | Moderate | TNF/inflammatory marker; plausible but composition-sensitive | 65 |
| 14 | GRAMD1A | Up in tumor; higher expression risky | Moderate | Cholesterol/ER contact-site biology; ccRCC evidence exists but mechanism immature | 62 |
| 15 | HHLA2 | Up in tumor; higher expression protective | Conflicting | Immune checkpoint with ccRCC-specific prognostic paradox | 60 |
| 16 | HIBCH | Down in tumor; higher expression protective | Weak-moderate | Valine/BCAA metabolism; plausible as metabolic-retention marker | 58 |
| 17 | TEK | Down in tumor; higher expression protective | Moderate but confounded | Endothelial/Tie2 angiogenesis marker; likely vascular composition | 55 |
| 18 | EMCN | Down in tumor; higher expression protective | Moderate but confounded | Endothelial glycocalyx marker; likely vascular composition | 52 |
| 19 | FHOD1 | Up in tumor; higher expression risky | Weak-moderate | EMT/actin-remodeling plausibility but little ccRCC-specific support | 50 |
| 20 | PODXL | Down in tumor; higher expression protective | Moderate but confounded | Kidney/endothelial/podocyte marker; renal survival direction plausible but tumor-cell story weak | 45 |
| 21 | FUT6 | Down in tumor; higher expression protective | Weak/conflicting | Glycosylation hypothesis; direction conflicts with many cancer fucosylation priors | 44 |
| 22 | LRBA | Down in tumor; higher expression protective | Weak | Immune trafficking/CTLA4 biology; likely immune-composition signal | 42 |
| 23 | CADPS2 | Down in tumor; higher expression protective | Weak | Vesicle/secretion gene; insufficient ccRCC rationale | 38 |
| 24 | IFFO1 | Up in tumor; higher expression risky | Weak | Sparse renal/cancer evidence; computational association only | 28 |

## Gene-by-Gene Review

### 1. KL

- Function: Klotho is a kidney-enriched membrane/secreted protein involved in phosphate/FGF23 biology, aging, oxidative stress, insulin/IGF signaling, Wnt/TGF-beta signaling, and renal tubular homeostasis.
- ccRCC/kidney/cancer roles: Strong kidney prior. Multiple ccRCC/RCC studies report KL downregulation in tumor and better outcome with higher KL expression. Klotho has tumor-suppressive literature in several cancers through growth-factor and EMT-related pathways.
- Observed direction: Tumor downregulation plus protective higher expression is highly plausible.
- Pathways: Kidney differentiation, mineral metabolism, PI3K/Akt, Wnt, EMT suppression, oxidative stress; indirectly linked to hypoxia/metabolic stress.
- Supporting literature: ccRCC/RCC studies report low KL in tumor and favorable prognosis with higher KL expression: https://pmc.ncbi.nlm.nih.gov/articles/PMC4860372/, https://pmc.ncbi.nlm.nih.gov/articles/PMC7657207/, https://pmc.ncbi.nlm.nih.gov/articles/PMC4912252/.
- Contradictions/limitations: A recent small protein-level pilot reported consistently low KL in ccRCC and limited grade separation, while TCGA mRNA survival was favorable: https://pubmed.ncbi.nlm.nih.gov/41240794/.
- Confounding: Strong candidate for retained normal renal epithelial program rather than active tumor-cell suppression in bulk data.
- Assessment: Likely tumor marker and likely prognostic marker; not enough for therapeutic-target language.
- Validation priority: HPA/IHC, CPTAC/FUSCC proteomics, E-MTAB-1980, independent FFPE TMA with adjacent normal.
- Evidence quality: Strong.
- Confidence: 88/100.

### 2. ACADM

- Function: Mitochondrial medium-chain acyl-CoA dehydrogenase required for fatty-acid beta-oxidation.
- ccRCC/kidney/cancer roles: ccRCC is characterized by lipid accumulation and suppressed mitochondrial oxidative metabolism. Low ACADM has been reported as a poor-prognosis and suppressive-TME-associated feature in ccRCC. ACADM also appears in fatty-acid metabolism risk signatures.
- Observed direction: Tumor downregulation plus protective higher expression is highly plausible.
- Pathways: Fatty-acid oxidation, mitochondrial metabolism, lipid handling, immune-metabolic coupling.
- Supporting literature: Low ACADM predicts poor prognosis in ccRCC: https://pmc.ncbi.nlm.nih.gov/articles/PMC11045743/. Fatty-acid metabolism signatures include ACADM/ACAT1: https://pmc.ncbi.nlm.nih.gov/articles/PMC9304894/. General ccRCC metabolic context: https://pmc.ncbi.nlm.nih.gov/articles/PMC6561085/.
- Contradictions/limitations: Fatty-acid oxidation can support survival in some cancers, so the direction is cancer-context dependent.
- Confounding: May mark retained proximal tubular mitochondrial function rather than a causal anti-aggressive program.
- Assessment: Likely prognostic marker; plausible tumor metabolic state marker.
- Validation priority: CPTAC protein abundance, metabolomics, isotope/FAO functional assays, external transcript cohorts.
- Evidence quality: Strong.
- Confidence: 84/100.

### 3. CRYL1

- Function: Crystallin lambda 1 is an enzyme-like protein linked to carbonyl metabolism and renal/hepatic metabolic programs.
- ccRCC/kidney/cancer roles: ccRCC-focused studies report CRYL1 downregulation, better survival with higher expression, and functional tumor-suppressive behavior in renal cancer cells. It has also been tied to metabolic and cuproptosis-related signatures.
- Observed direction: Tumor downregulation plus protective higher expression is very plausible.
- Pathways: Sugar/fat/amino-acid metabolism, apoptosis/proliferation, possibly cuproptosis-associated metabolic state.
- Supporting literature: CRYL1 reported as independent/prognostic ccRCC marker and tumor suppressor-like gene: https://pmc.ncbi.nlm.nih.gov/articles/PMC9563328/, https://pmc.ncbi.nlm.nih.gov/articles/PMC10946081/, https://pubmed.ncbi.nlm.nih.gov/38497139/.
- Contradictions/limitations: Some CRYL1 work is bioinformatics-heavy and partially overlaps public datasets similar to this pipeline.
- Confounding: Retained renal epithelial metabolic identity.
- Assessment: Likely prognostic marker; plausible tumor marker.
- Validation priority: IHC/TMA, CPTAC, E-MTAB-1980, functional rescue/knockdown in ccRCC cell lines.
- Evidence quality: Strong.
- Confidence: 84/100.

### 4. ACAT1

- Function: Mitochondrial acetyl-CoA acetyltransferase involved in ketone-body metabolism, branched-chain amino-acid catabolism, and acetyl-CoA handling.
- ccRCC/kidney/cancer roles: ccRCC studies identify decreased ACAT1 and related ketone metabolism enzymes in tumor; lower expression associates with worse outcome, and overexpression suppresses proliferation/migration in renal cancer models.
- Observed direction: Tumor downregulation plus protective higher expression is highly plausible.
- Pathways: Ketone body metabolism, fatty-acid/acetyl-CoA metabolism, mitochondrial function, AMPK/lipid metabolism.
- Supporting literature: Ketone metabolism study: https://pmc.ncbi.nlm.nih.gov/articles/PMC6928137/. ACAT1 co-expression/functional support: https://pmc.ncbi.nlm.nih.gov/articles/PMC6795108/. Fatty-acid enzyme prognostic analysis: https://pmc.ncbi.nlm.nih.gov/articles/PMC6856888/.
- Contradictions/limitations: ACAT1 biology differs across cancers; ketone utilization can be tumor-promoting in some contexts.
- Confounding: Renal mitochondrial differentiation and tumor purity.
- Assessment: Likely prognostic marker and biologically coherent metabolic candidate.
- Validation priority: Protein/metabolite validation, matched normal/tumor IHC, external transcript/proteomic cohorts.
- Evidence quality: Strong.
- Confidence: 83/100.

### 5. DDC

- Function: Aromatic L-amino-acid decarboxylase that converts L-DOPA to dopamine and 5-hydroxytryptophan to serotonin; broader amino-acid/amine metabolism.
- ccRCC/kidney/cancer roles: DDC is strongly downregulated in RCC/ccRCC, has diagnostic value, and a multi-omics ccRCC paper supports prognostic and immune-microenvironment relevance.
- Observed direction: Tumor downregulation plus protective higher expression is plausible and externally supported.
- Pathways: Amino-acid metabolism, neurotransmitter/biogenic amine metabolism, PI3K/Akt and immune-microenvironment associations reported in ccRCC.
- Supporting literature: Multi-omics ccRCC DDC study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9760914/. Earlier RCC diagnostic report: https://www.sciencedirect.com/science/article/abs/pii/S0009912015000612.
- Contradictions/limitations: Mechanism in renal epithelial cancer is not as canonical as KL/ACAT1/ACADM; some immune associations may be secondary.
- Confounding: Kidney-lineage/metabolic retention, tumor purity, endocrine-like expression artifacts.
- Assessment: Likely prognostic marker; likely tumor marker. Mechanistic claims should remain cautious.
- Validation priority: CPTAC/FUSCC proteomics, external RNA cohorts, IHC by grade/stage, multivariable survival with immune covariates.
- Evidence quality: Strong.
- Confidence: 82/100.

### 6. PANK1

- Function: Pantothenate kinase 1 catalyzes the rate-limiting step in coenzyme A biosynthesis.
- ccRCC/kidney/cancer roles: ccRCC studies report lower PANK1 expression in tumors, poor survival with low expression, and functional effects on apoptosis, invasion, migration, EMT, and mitochondrial metabolism.
- Observed direction: Tumor downregulation plus protective higher expression is plausible.
- Pathways: CoA metabolism, mitochondrial metabolism, apoptosis, EMT/invasion, immune infiltration.
- Supporting literature: TCGA-based prognostic study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9372198/. Functional KIRC study: https://pubmed.ncbi.nlm.nih.gov/39196459/.
- Contradictions/limitations: Much evidence is public-data-driven; functional data should be checked for model specificity.
- Confounding: Retained normal metabolic state, stage/grade correlation, tumor purity.
- Assessment: Likely prognostic marker; biologically plausible metabolic candidate.
- Validation priority: Protein-level validation, external survival cohorts, promoter methylation, functional rescue.
- Evidence quality: Strong.
- Confidence: 78/100.

### 7. TCIRG1

- Function: V-ATPase V0 subunit a3; involved in proton transport, organelle acidification, lysosomal function, osteoclast biology, and immune-cell biology.
- ccRCC/kidney/cancer roles: ccRCC studies report overexpression, poor outcome, immune infiltration associations, glycolysis/AKT/mTOR links, and reduced migration after knockdown.
- Observed direction: Tumor upregulation plus higher-risk expression is plausible.
- Pathways: Lysosomal acidification, glycolysis, AKT/mTOR, immune infiltration, Treg/CD8 contexts, immunotherapy relevance.
- Supporting literature: Integrative ccRCC TCIRG1 study: https://pmc.ncbi.nlm.nih.gov/articles/PMC9558535/. Glycolysis/immune study: https://pmc.ncbi.nlm.nih.gov/articles/PMC10468907/.
- Contradictions/limitations: TCIRG1 may be coming from immune/myeloid/osteoclast-like cells, not tumor cells.
- Confounding: Immune infiltration, macrophage content, necrosis, grade, therapy-response subgrouping.
- Assessment: Likely prognostic marker; uncertain tumor-cell marker.
- Validation priority: Single-cell/spatial RNA, multiplex IHC with myeloid/T cell markers, IMmotion/ICI cohorts.
- Evidence quality: Strong, but cell-source caveat is important.
- Confidence: 78/100.

### 8. DBT

- Function: E2 subunit of branched-chain alpha-ketoacid dehydrogenase complex; required for branched-chain amino-acid catabolism.
- ccRCC/kidney/cancer roles: Multiple ccRCC-focused reports describe DBT downregulation, favorable prognosis with higher expression, lipid accumulation effects, and Hippo/YAP pathway involvement.
- Observed direction: Tumor downregulation plus protective higher expression is plausible.
- Pathways: Branched-chain amino-acid metabolism, lipid accumulation, m6A regulation, Hippo/YAP signaling, immune infiltration.
- Supporting literature: DBT prognostic ccRCC analysis: https://pmc.ncbi.nlm.nih.gov/articles/PMC9589218/. Functional m6A/ANXA2/YAP/Hippo study: https://pmc.ncbi.nlm.nih.gov/articles/PMC10091108/.
- Contradictions/limitations: Some evidence is derived from public databases; BCAA metabolism can have mixed cancer roles.
- Confounding: Renal metabolic differentiation and tumor purity.
- Assessment: Likely prognostic marker; plausible metabolic tumor-state marker.
- Validation priority: Proteomics/metabolomics, independent transcript cohorts, cell-line functional perturbation.
- Evidence quality: Strong.
- Confidence: 77/100.

### 9. CLCN5

- Function: Cl-/H+ exchanger enriched in renal proximal tubule endosomes; critical for receptor-mediated endocytosis and Dent disease biology.
- ccRCC/kidney/cancer roles: Strong kidney biology prior. A recent ccRCC study reports CLCN5 inhibits tumorigenesis and fatty-acid accumulation through EHHADH regulation.
- Observed direction: Tumor downregulation plus protective higher expression is plausible.
- Pathways: Proximal tubule endocytosis, fatty-acid metabolism via EHHADH, renal epithelial differentiation.
- Supporting literature: Kidney CLCN5 function review: https://pmc.ncbi.nlm.nih.gov/articles/PMC7402612/. ccRCC CLCN5/EHHADH study: https://pmc.ncbi.nlm.nih.gov/articles/PMC12320790/.
- Contradictions/limitations: As a proximal tubule marker, lower tumor expression may reflect dedifferentiation rather than specific suppressor activity.
- Confounding: Normal-tubule contamination, tumor purity, retained renal epithelium.
- Assessment: Likely tumor marker; plausible prognostic marker, but mechanism needs independent replication.
- Validation priority: Spatial transcriptomics/IHC separating tumor cells from entrapped tubules, CPTAC protein, external survival cohorts.
- Evidence quality: Moderate-strong.
- Confidence: 75/100.

### 10. CYFIP2

- Function: Cytoplasmic FMR1-interacting protein involved in actin cytoskeleton regulation, WAVE regulatory complex biology, Rac/Rho signaling, apoptosis, and p53-responsive programs.
- ccRCC/kidney/cancer roles: A ccRCC study reports CYFIP2 downregulation and prognostic/immune-infiltration relevance. Rho-GTPase signature work also reports CYFIP2 underexpression in ccRCC.
- Observed direction: Tumor downregulation plus protective higher expression is plausible, but less kidney-specific.
- Pathways: Cytoskeleton, apoptosis, Rho-GTPase signaling, immune infiltration; weak direct hypoxia/metabolism connection despite pipeline annotation.
- Supporting literature: CYFIP2 ccRCC paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC8523908/. Rho-GTPase ccRCC signature includes CYFIP2/IQGAP2 underexpression: https://pmc.ncbi.nlm.nih.gov/articles/PMC11861022/.
- Contradictions/limitations: Biology is broad; tumor-suppressor interpretation is not yet as mature as KL/ACAT1/CRYL1.
- Confounding: Cell composition, cytoskeletal state, immune/stromal admixture.
- Assessment: Likely prognostic marker; tumor-marker status uncertain.
- Validation priority: Protein/IHC, single-cell localization, independent multivariable cohorts.
- Evidence quality: Moderate.
- Confidence: 73/100.

### 11. C1QTNF6

- Function: Secreted C1q/TNF-related protein with adipokine-like, inflammatory, metabolic, and immune-modulatory associations.
- ccRCC/kidney/cancer roles: ccRCC-specific biomarker work and HPA survival data support high C1QTNF6 as unfavorable in KIRC. Other cancers also link C1QTNF6 to progression and immune infiltration.
- Observed direction: Tumor upregulation plus higher-risk expression is plausible.
- Pathways: Secreted inflammatory signaling, immune infiltration, metabolic inflammation, possible EMT/proliferation depending on tumor context.
- Supporting literature: ccRCC biomarker report page: https://scholars.lib.ntu.edu.tw/entities/publication/0a0acd11-c321-4de0-8d29-60123735d9bc. HPA KIRC survival page: https://www.proteinatlas.org/ENSG00000133466-C1QTNF6/cancer/renal%2Bcancer.
- Contradictions/limitations: Mechanistic ccRCC literature appears thinner than for the strongest metabolic genes.
- Confounding: Secreted-cell source, immune/stromal content, adiposity/metabolic host factors.
- Assessment: Likely prognostic marker; tumor marker uncertain.
- Validation priority: Spatial expression, serum/protein assays, external cohorts with BMI/metabolic covariates.
- Evidence quality: Moderate.
- Confidence: 70/100.

### 12. IQGAP2

- Function: Scaffold protein regulating cytoskeleton, adhesion, small GTPase signaling, and cell-cell signaling. Often considered tumor suppressor-like relative to IQGAP1/3 in several cancers.
- ccRCC/kidney/cancer roles: Rho-GTPase ccRCC work reports IQGAP2 underexpression in tumors and higher expression in normal renal epithelial cells; broader cancer literature supports tumor-suppressor roles.
- Observed direction: Tumor downregulation plus protective higher expression is plausible.
- Pathways: Cytoskeleton, Rho-GTPase signaling, cell adhesion, possible EMT restraint.
- Supporting literature: ccRCC Rho-GTPase signature including IQGAP2: https://pmc.ncbi.nlm.nih.gov/articles/PMC11861022/.
- Contradictions/limitations: Direct ccRCC functional evidence is less deep than the metabolic candidates.
- Confounding: Renal epithelial differentiation, stromal/immune composition.
- Assessment: Likely prognostic marker; tumor-cell mechanism plausible but not proven.
- Validation priority: IHC/spatial, protein-level validation, functional perturbation in ccRCC models.
- Evidence quality: Moderate.
- Confidence: 68/100.

### 13. TNFAIP2

- Function: TNF-alpha inducible protein implicated in inflammation, endothelial/immune activation, angiogenesis, migration, and cancer progression.
- ccRCC/kidney/cancer roles: Review-level evidence reports high TNFAIP2 mRNA associated with shorter survival in KIRC. The inflammatory direction matches high-risk expression in this pipeline.
- Observed direction: Tumor upregulation plus higher-risk expression is plausible.
- Pathways: TNF/inflammation, immune regulation, angiogenesis, migration/invasion.
- Supporting literature: TNFAIP2 cancer review notes KIRC survival association: https://pmc.ncbi.nlm.nih.gov/articles/PMC6201362/.
- Contradictions/limitations: Inflammation markers may be non-specific; TNFAIP2 can have different prognostic directions across cancers.
- Confounding: Myeloid/endothelial content, inflammation, necrosis, treatment-related immune state.
- Assessment: Likely prognostic marker; uncertain tumor-cell marker.
- Validation priority: Single-cell/spatial RNA, macrophage/endothelial covariate adjustment, ICI datasets.
- Evidence quality: Moderate.
- Confidence: 65/100.

### 14. GRAMD1A

- Function: ER/plasma-membrane contact-site protein involved in cholesterol transport/sensing.
- ccRCC/kidney/cancer roles: A KIRC/ccRCC biomarker paper reports GRAMD1A as associated with immune infiltration and unfavorable prognosis; HCC studies support pro-growth roles in another cancer context.
- Observed direction: Tumor upregulation plus higher-risk expression is plausible given ccRCC lipid/cholesterol biology.
- Pathways: Cholesterol transport, ER contact sites, lipid metabolism, immune infiltration; possible STAT5 link from non-kidney cancer.
- Supporting literature: KIRC GRAMD1A biomarker/immune infiltration paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9293538/. HCC growth/stemness support: https://pmc.ncbi.nlm.nih.gov/articles/PMC5009375/.
- Contradictions/limitations: Mechanistic ccRCC evidence is immature; lipid-contact biology is plausible but not established.
- Confounding: Lipid-rich tumor state, immune infiltration, grade.
- Assessment: Likely prognostic marker; biologically interesting lipid-adaptation candidate.
- Validation priority: CPTAC lipid/protein correlations, spatial RNA, independent survival cohorts.
- Evidence quality: Moderate.
- Confidence: 62/100.

### 15. HHLA2

- Function: B7-family immune checkpoint ligand interacting with KIR3DL3 and other receptors; can regulate T-cell function.
- ccRCC/kidney/cancer roles: HHLA2 is high in KIRC/ccRCC. Literature is contradictory: one ccRCC study reports high HHLA2 associated with favorable survival and CD8 trends, while another reports overexpression associated with poor survival, advanced stage, and EMT.
- Observed direction: Tumor upregulation plus protective higher expression is plausible only under an immune-inflamed/favorable-T cell interpretation; it is not universally consistent.
- Pathways: Immune checkpoint biology, T-cell regulation, immune infiltration, possible EMT.
- Supporting literature: Favorable-prognosis report: https://pmc.ncbi.nlm.nih.gov/articles/PMC7248229/. Conflicting poor-prognosis report: https://pmc.ncbi.nlm.nih.gov/articles/PMC6469208/. PD-L1/HHLA2 co-expression poor-prognosis context: https://pmc.ncbi.nlm.nih.gov/articles/PMC7057441/. Regulation in kidney cancer/myeloid cells: https://pubmed.ncbi.nlm.nih.gov/37891555/.
- Contradictions/limitations: The literature is genuinely conflicting and probably depends on compartment, assay, cutoff, PD-L1 coexpression, TIL state, and cohort.
- Confounding: TIL abundance, tumor-cell versus myeloid expression, checkpoint-exhaustion state, PH p near threshold in this pipeline.
- Assessment: Likely tumor immune marker; prognostic direction uncertain.
- Validation priority: Multiplex IHC/spatial profiling with CD8, PD-L1, KIR3DL3, myeloid markers; ICI cohorts.
- Evidence quality: Conflicting.
- Confidence: 60/100.

### 16. HIBCH

- Function: Mitochondrial 3-hydroxyisobutyryl-CoA hydrolase in valine catabolism.
- ccRCC/kidney/cancer roles: Direct ccRCC evidence is limited, but HIBCH appears in metabolic/immunologic prognostic signatures and fits the broader loss of mitochondrial/BCAA metabolism seen in ccRCC.
- Observed direction: Tumor downregulation plus protective higher expression is plausible but not strongly validated.
- Pathways: Branched-chain amino-acid/valine metabolism, mitochondrial metabolic state.
- Supporting literature: m6A/lactylation-related ccRCC prognostic work includes HIBCH context: https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1225023/pdf. General ccRCC metabolic hallmarks: https://pmc.ncbi.nlm.nih.gov/articles/PMC6561085/.
- Contradictions/limitations: Little direct single-gene ccRCC biology.
- Confounding: Renal metabolic differentiation and normal-tubule contamination.
- Assessment: Plausible metabolic passenger/prognostic marker; not a standalone manuscript lead.
- Validation priority: External RNA cohorts, CPTAC protein, metabolomic correlation with BCAA pathway.
- Evidence quality: Weak-moderate.
- Confidence: 58/100.

### 17. TEK

- Function: Tie2 receptor tyrosine kinase, a central endothelial angiopoietin receptor involved in vascular maturation, endothelial quiescence, angiogenesis, and Tie2-expressing monocytes.
- ccRCC/kidney/cancer roles: Angiogenesis is core to VHL/HIF ccRCC biology. TEK/Tie2 is endothelial rather than renal tumor-cell intrinsic. A ccRCC tumor-vasculature study supports Tie2 relevance.
- Observed direction: Tumor downregulation plus protective higher expression is biologically possible if it marks organized/less aggressive vasculature, but it is also composition-heavy.
- Pathways: Angiogenesis, endothelial biology, vascular maturation, immune-myeloid angiogenic cells.
- Supporting literature: Tie2/B7-H3 tumor vasculature ccRCC study: https://pmc.ncbi.nlm.nih.gov/articles/PMC5692197/. General ccRCC VHL/HIF/VEGF context: https://pmc.ncbi.nlm.nih.gov/articles/PMC8113694/.
- Contradictions/limitations: Angiogenesis can support tumor growth, so "higher TEK protective" is not straightforward.
- Confounding: Endothelial fraction, vessel normalization, normal kidney vascular content, stromal admixture.
- Assessment: Likely vascular-composition marker; uncertain prognostic marker.
- Validation priority: Spatial endothelial quantification, microvessel density, CD31/EMCN/TEK co-staining.
- Evidence quality: Moderate but confounded.
- Confidence: 55/100.

### 18. EMCN

- Function: Endomucin is an endothelial sialomucin/glycocalyx protein involved in vascular identity and leukocyte-endothelial adhesion restraint.
- ccRCC/kidney/cancer roles: EMCN appears in VHL-mutant ccRCC subcluster/prognosis work and mucin-signature analyses. It is most plausibly an endothelial marker.
- Observed direction: Tumor downregulation plus protective higher expression may indicate preserved vascular architecture or lower aggressive dedifferentiation, but bulk interpretation is weak.
- Pathways: Endothelial biology, vascular integrity, leukocyte adhesion, angiogenesis/immune trafficking.
- Supporting literature: PBX1/EMCN/ERG VHL-mutant ccRCC prognostic paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC9142578/. MUCIN signature including EMCN: https://pmc.ncbi.nlm.nih.gov/articles/PMC8419780/.
- Contradictions/limitations: Not a tumor-cell gene; prognosis could be vascular composition.
- Confounding: Endothelial content, tissue sampling, normal-vessel admixture.
- Assessment: Likely vascular marker; prognostic value uncertain without spatial validation.
- Validation priority: Spatial transcriptomics, endothelial normalization, IHC with CD31/ERG.
- Evidence quality: Moderate but confounded.
- Confidence: 52/100.

### 19. FHOD1

- Function: Formin-family actin regulator involved in cytoskeletal remodeling, stress fibers, migration, invasion, and invadopodia.
- ccRCC/kidney/cancer roles: Strong general cancer/EMT plausibility; limited direct ccRCC evidence. EMT and invasion are relevant to aggressive RCC biology.
- Observed direction: Tumor upregulation plus higher-risk expression is plausible.
- Pathways: EMT, actin cytoskeleton, ECM degradation, invasion, PI3K/ZEB/Snail-like programs.
- Supporting literature: FHOD1 EMT/cancer migration study: https://pmc.ncbi.nlm.nih.gov/articles/PMC3784416/. FHOD1 tumor/TME review: https://pmc.ncbi.nlm.nih.gov/articles/PMC12069282/.
- Contradictions/limitations: Lack of ccRCC-specific support; stromal/mesenchymal cells may drive bulk expression.
- Confounding: Stromal content, sarcomatoid/EMT fraction, fibroblasts, endothelial/pericyte contamination.
- Assessment: Biologically plausible risk marker, but weak manuscript candidate without ccRCC validation.
- Validation priority: ccRCC spatial/IHC, sarcomatoid annotation, external survival cohorts.
- Evidence quality: Weak-moderate.
- Confidence: 50/100.

### 20. PODXL

- Function: Podocalyxin-like protein, a CD34-related sialomucin expressed in podocytes, endothelium, and some stem/cancer cells; regulates cell adhesion and anti-adhesion.
- ccRCC/kidney/cancer roles: In many cancers, high PODXL is associated with invasion and poor prognosis. Kidney is a special case because PODXL is also a normal glomerular/endothelial marker, and HPA reports favorable KIRC prognosis with higher PODXL transcript expression.
- Observed direction: Tumor downregulation plus protective higher expression is plausible as a renal/endothelial differentiation or organized-vasculature signal, but not convincing as a tumor-cell-intrinsic protective mechanism.
- Pathways: Cell adhesion, EMT/invasion, endothelial/podocyte biology, glycosylation.
- Supporting literature: HPA renal cancer page reports favorable KIRC prognosis: https://www.proteinatlas.org/ENSG00000128567-PODXL/cancer/renal%2Bcancer. General PODXL poor-prognosis meta-analysis in many cancers: https://pmc.ncbi.nlm.nih.gov/articles/PMC5581042/.
- Contradictions/limitations: General cancer literature often treats PODXL as invasion-associated; renal bulk direction likely depends on compartment.
- Confounding: Podocyte/endothelial admixture, glomerular content in normal/adjacent tissue, tumor-cell versus endothelial expression.
- Assessment: Possible prognostic marker, but should not be presented as a tumor-cell marker without spatial/protein validation.
- Validation priority: Spatial/IHC compartment separation; compare tumor-cell PODXL vs endothelial/podocyte PODXL.
- Evidence quality: Moderate but confounded.
- Confidence: 45/100.

### 21. FUT6

- Function: Alpha-1,3/1,4 fucosyltransferase involved in terminal fucosylation and Lewis antigen synthesis.
- ccRCC/kidney/cancer roles: Glycosylation/fucosylation is relevant to urologic cancers, cell adhesion, immune recognition, and metastasis. Direct ccRCC evidence for FUT6 is thin.
- Observed direction: Tumor downregulation plus protective higher expression is not strongly intuitive because fucosylation is often pro-adhesive/pro-metastatic in cancer.
- Pathways: Glycosylation, selectin ligands, cell adhesion, PI3K/Akt in other urologic contexts.
- Supporting literature: Urologic cancer fucosylation review: https://pmc.ncbi.nlm.nih.gov/articles/PMC8708646/.
- Contradictions/limitations: Many fucosylation studies imply tumor-promoting roles for FUT-family activity; direct FUT6 ccRCC evidence is lacking.
- Confounding: Blood/endothelial/immune glycan biology, adjacent normal kidney composition, platform artifacts.
- Assessment: Uncertain; likely passenger/composition marker unless independently validated.
- Validation priority: Glycoproteomics, spatial expression, validation against FUT3/4/7 and selectin ligand markers.
- Evidence quality: Weak/conflicting.
- Confidence: 44/100.

### 22. LRBA

- Function: Intracellular trafficking/adaptor protein important for immune regulation, CTLA4 recycling, lysosomal trafficking, and immune tolerance.
- ccRCC/kidney/cancer roles: Strong immune biology generally, but ccRCC-specific evidence is weak and appears in broader immune-signature contexts.
- Observed direction: Tumor downregulation plus protective higher expression could reflect beneficial immune competence, but the cell source is unresolved.
- Pathways: Immune regulation, CTLA4 trafficking, lysosomal vesicle biology, cytotoxic lymphocyte function.
- Supporting literature: ccRCC immune-function signature context mentions LRBA: https://pmc.ncbi.nlm.nih.gov/articles/PMC12611398/. TME context for ccRCC immune confounding: https://pmc.ncbi.nlm.nih.gov/articles/PMC6774890/.
- Contradictions/limitations: LRBA has complex immunosuppressive and immune-function roles; direction is not straightforward.
- Confounding: Lymphocyte subset composition, CTLA4/Treg state, immune activation/exhaustion.
- Assessment: Uncertain immune-composition marker; not a core tumor biomarker.
- Validation priority: Single-cell/spatial immune localization, CTLA4/Treg/cytotoxic-cell adjustment, ICI datasets.
- Evidence quality: Weak.
- Confidence: 42/100.

### 23. CADPS2

- Function: Calcium-dependent secretion activator involved in dense-core vesicle exocytosis and regulated secretion.
- ccRCC/kidney/cancer roles: Direct ccRCC evidence is sparse. Cancer relevance is not well established for renal tumors.
- Observed direction: Tumor downregulation plus protective higher expression is possible but weakly grounded.
- Pathways: Vesicle secretion, calcium-regulated exocytosis; no strong direct link to hypoxia, angiogenesis, metabolism, ECM/EMT, or immune regulation in ccRCC.
- Supporting literature: No compelling ccRCC-specific support identified in this review.
- Contradictions/limitations: Lack of kidney/cancer mechanism makes it fragile.
- Confounding: Neural/endocrine-like expression, stromal/normal tissue signal, annotation artifacts.
- Assessment: Likely passenger or retained tissue-state association.
- Validation priority: First validate expression/protein localization; do not discuss mechanistically before that.
- Evidence quality: Weak.
- Confidence: 38/100.

### 24. IFFO1

- Function: Intermediate filament family orphan protein; limited functional characterization. Emerging evidence links it to nuclear/mitochondrial organization and metabolism in non-renal cancer contexts.
- ccRCC/kidney/cancer roles: No strong ccRCC-specific evidence identified. Very recent non-renal work suggests tumor-suppressive biology in breast cancer, which does not support a confident renal interpretation.
- Observed direction: Tumor upregulation plus higher-risk expression is computationally coherent but biologically under-supported.
- Pathways: Structural/nuclear/mitochondrial organization; no strong ccRCC pathway placement.
- Supporting literature: Non-renal breast cancer mechanistic paper: https://www.nature.com/articles/s41389-026-00609-1.pdf.
- Contradictions/limitations: Sparse literature and no clear kidney biology.
- Confounding: Annotation artifacts, low-expression instability, stromal/structural content, unmeasured tumor state.
- Assessment: Likely passenger/computational association until independently validated.
- Validation priority: Confirm expression robustness, protein localization, external survival cohorts.
- Evidence quality: Weak.
- Confidence: 28/100.

## Top 5 Strongest Manuscript Candidates

1. KL
2. ACADM
3. CRYL1
4. ACAT1
5. DDC

Rationale: These have the best combination of coherent direction, renal/metabolic plausibility, and ccRCC-focused literature. PANK1, TCIRG1, and DBT are close alternates. If the manuscript wants a risk-direction candidate, TCIRG1 is the strongest risk gene.

## Top 5 Most Biologically Interesting Candidates

1. HHLA2: Immune checkpoint with contradictory ccRCC survival direction; worth discussing as a compartment/cutoff-dependent immune paradox, not as a simple protective gene.
2. TCIRG1: Acidification/glycolysis/immune infiltration gene with risk direction and functional support.
3. GRAMD1A: Lipid/cholesterol contact-site biology, aligning with the lipid-rich ccRCC phenotype.
4. C1QTNF6: Secreted inflammatory/metabolic signal with unfavorable direction.
5. CLCN5: Proximal tubule endosomal/lipid metabolism connection; potentially bridges retained renal identity and metabolic restraint.

## Candidates to Remove or Downgrade From Final Manuscript Claims

Remove from the final highlighted candidate list unless new independent validation is added:

- IFFO1: too little ccRCC/kidney/cancer support.
- CADPS2: weak mechanistic and ccRCC literature support.
- LRBA: immune-composition explanation dominates; not enough ccRCC-specific support.
- FUT6: glycosylation is interesting, but FUT6-specific ccRCC evidence is weak and direction is not convincing.

Downgrade to "composition/pathway flags" rather than core candidates:

- TEK and EMCN: useful vascular biology signals, but likely endothelial composition rather than tumor-cell intrinsic prognosis.
- PODXL: useful renal/endothelial marker, but heavily compartment-confounded and not safe as a tumor-cell biomarker.
- FHOD1: plausible EMT/invasion risk marker, but needs ccRCC-specific validation.

## Recurring Biological Themes

1. Loss of renal epithelial metabolic differentiation: KL, CLCN5, DDC, CRYL1, ACADM, ACAT1, DBT, HIBCH, PANK1.
2. Suppression of mitochondrial oxidative, fatty-acid, ketone, CoA, and branched-chain amino-acid metabolism: ACADM, ACAT1, DBT, HIBCH, PANK1.
3. Immune/checkpoint/inflammatory risk or paradox: HHLA2, TCIRG1, TNFAIP2, LRBA, C1QTNF6.
4. Vascular/endothelial composition and angiogenic state: TEK, EMCN, PODXL.
5. Cytoskeletal/EMT/invasion signals: FHOD1, IQGAP2, CYFIP2, PODXL.
6. Bulk RNA-seq composition risk: many protective downregulated genes may mark retained normal kidney/proximal tubule biology rather than causal tumor suppression.

## Assessment of the Central Hypothesis

Central hypothesis: "Are the most reproducibly dysregulated genes also the most prognostically important, or is survival driven by a narrower hypoxia-associated adaptation program?"

Conclusion: The results support the broader hypothesis that reproducible dysregulation and prognostic importance are not equivalent. However, the evidence does not support a simple "hypoxia-only" program.

The strongest biological interpretation is narrower and more precise:

- Survival signal appears enriched for a renal metabolic differentiation/adaptation axis: preserved mitochondrial, fatty-acid, ketone, CoA, and BCAA metabolism is associated with lower hazard.
- This axis is compatible with VHL/HIF-driven ccRCC biology because hypoxia signaling suppresses oxidative metabolism and reshapes lipid handling, but most top candidates are not canonical HIF target genes.
- A second layer involves immune and vascular microenvironment states, including checkpoint biology and endothelial composition.
- Therefore, the manuscript should frame the finding as "hypoxia-linked metabolic dedifferentiation/adaptation plus immune/vascular composition" rather than "hypoxia-associated adaptation program" alone.

## Independent Validation Resources

Recommended validation resources before strong claims:

- E-MTAB-1980 ccRCC survival transcriptomics.
- CPTAC ccRCC proteogenomics and Fudan/FUSCC proteomics where available.
- Human Protein Atlas renal cancer pathology/IHC pages.
- cBioPortal TCGA-KIRC for copy-number/mutation/correlation checks.
- External GEO ccRCC expression cohorts not used in model training, such as GSE53000 and GSE105261, where appropriate.
- Spatial/single-cell resources: TISCH2, Kidney Tumor Atlas/single-cell ccRCC studies, and multiplex IHC/TMA validation.
- Immunotherapy cohorts for immune candidates: IMmotion and CheckMate-derived expression datasets if accessible.

## Bottom Line for Manuscript Language

Use "candidate prognostic association" for all 24. Use stronger language only for KL, ACADM, CRYL1, ACAT1, DDC, PANK1, TCIRG1, DBT, and CLCN5, and even there avoid causal language. The strongest manuscript message is not that these 24 genes are all biomarkers. It is that the evidence funnel compresses broad reproducible dysregulation into a smaller set dominated by renal metabolic retention and selected immune/vascular states.
