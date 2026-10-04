# Scientific source verification

Sources below were checked through the browsing tool on **2026-10-01**, the user's local date. This is a literature verification date, not the unknown original retrieval date of cached study data. Primary sources support the added related-work claims. No citations were invented.

| Citation | Verified source | Supported use |
| --- | --- | --- |
| McShane et al., REMARK, 2005, BJC 93:387–391 | [Publisher](https://www.nature.com/articles/6602678), DOI 10.1038/sj.bjc.6602678 | Reporting patient selection, assay methods, covariates, and validation. |
| Simon et al., 2003, JNCI 95:14–18 | [NCI full paper](https://brb.nci.nih.gov/techreport/JNCICommentary.pdf), [PubMed](https://pubmed.ncbi.nlm.nih.gov/12509396/), DOI 10.1093/jnci/95.1.14 | Feature selection belongs inside validation of a research procedure. |
| Brooks et al., ClearCode34, 2014, Eur Urol 66:77–84 | [Full paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4058355/), [PubMed](https://pubmed.ncbi.nlm.nih.gov/24613583/), DOI 10.1016/j.eururo.2014.02.035 | Original subtype classifier and clinical-specimen evaluation; HYDRA's signed mean is an adaptation. Missing author Oishee Sen restored in bibliography. |
| Love et al., DESeq2, 2014 | [Publisher](https://link.springer.com/article/10.1186/s13059-014-0550-8) | Count-based differential-expression method in retained supplement. |
| Colaprico et al., TCGAbiolinks, 2016 | [Full paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4856967/) | Retained TCGA retrieval workflow. |
| Zhu et al., apeglm, 2019 | [Full paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6581436/) | Shrinkage of differential-expression effects; not the survival-test engine. |
| Ritchie et al., limma, 2015 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/25605792/) | Retained microarray differential-expression method. |
| Leek et al., sva, 2012 | [Full paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3307112/) | Assessment of unwanted variation; zero surrogate variables does not establish successful batch removal. |
| Simon et al., glmnet Cox, 2011 | [Journal](https://www.jstatsoft.org/article/view/v039i05) | Retained ridge-Cox method. |
| Aran et al., purity, 2015 | [Publisher](https://www.nature.com/articles/ncomms9971) | Consensus purity and transcriptomic confounding; retained bibliography article number 8971 is correct. |
| Fernández-Sanromán et al., TRACERx, 2025 | [Full paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11873726/) | Retained multiregion source. |
| Braun et al., CheckMate RNA, 2020 | [PubMed](https://pubmed.ncbi.nlm.nih.gov/32472114/) | Retained treatment-interaction input. First ten author names corrected against the bibliographic record. |
| Cox, 1972 | [Publisher](https://academic.oup.com/jrsssb/article/34/2/187/7027194), [original paper](https://web.stanford.edu/~lutian/coursepdf/cox1972paper.pdf), DOI 10.1111/j.2517-6161.1972.tb00899.x | Retained proportional-hazards regression method. |

Bibliographic existence checks do not certify every interpretation in older project notes. The main scientific results are supported by study artifacts and the executable claim/ledger checks. The retained supplement describes the all-gene Cox sensitivity accurately: apeglm estimates are differential-expression annotations, not its survival-test engine.

## Limitation-audit source checks, 2026-10-04

The [official GDC barcode documentation](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/) identifies the two-character field after `TCGA-` as tissue source site. Stage 46 uses that field for collection-site baseline-hazard strata. It does not treat site as a measured technical batch. The same documentation defines the four-character plate field used by stage 49 as a technical-processing proxy, with full-rank designs checked locally.

The [DESeq2 primary source](https://github.com/thelovelab/DESeq2/blob/devel/R/core.R) implements median ratios and re-centers their geometric mean when a reference is supplied. The installed DESeq2 1.52.0 function was inspected directly and the behavior reproduced in `test_frozen_normalization.R`. A one-sample supplied-reference call returns factor one; the corrected helper retains the training scale and passes batch-invariance and depth-scaling tests. The live development source supports the mechanism; installed runtime behavior and study corrections are verified locally rather than assumed identical to a moving branch.

White and Royston's [primary imputation study](https://pmc.ncbi.nlm.nih.gov/articles/PMC2998703/) was consulted when defining the missing-covariate threat. No multiple-imputation result is claimed or citation added to imply that the eight deterministic completion scenarios implement that method.
