#!/usr/bin/env python3

from pathlib import Path
import textwrap
import pandas as pd

PROJECT_DIR = Path("/mnt/d/thymic_apc_atlas_project")
DOC_DIR = PROJECT_DIR / "docs"
FIG_DIR = PROJECT_DIR / "results/figures/manuscript_final_panels"
DOC_DIR.mkdir(parents=True, exist_ok=True)

manuscript = r"""
# Single-cell and spatial reconstruction of human thymic stromal–immune antigen-presenting niches across early-life development

## Abstract

**Background:** The human thymus supports T-cell education through specialised epithelial, stromal and haematopoietic antigen-presenting compartments. Although thymic epithelial cells are central to classical models of positive and negative selection, the spatial organisation of broader stromal–immune antigen-presenting niches in human thymus remains incompletely resolved, particularly in early-life tissue contexts.

**Methods:** We performed an integrative analysis of public human thymus single-cell and spatial transcriptomic datasets. The Bautista human thymic stromal single-cell dataset was used as a reference for thymic epithelial and stromal architecture, while the Heimli paediatric thymus CITE-seq and Visium dataset was used to resolve immune antigen-presenting, epithelial, thymocyte and spatial niche programmes. Analyses included quality control, balanced downsampling, UMAP embedding, Leiden clustering, RNA-only marker ranking for CITE-seq interpretation, marker-programme scoring, manual conservative annotation, Visium spot-level signature scoring and spatial co-pattern analysis.

**Results:** The Bautista dataset resolved broad stromal and epithelial compartments, including TEC, fibroblast, vascular, vascular-MHCII-associated and perivascular states, supporting its use as a stromal/TEC architecture reference. The Heimli dataset resolved paediatric thymic immune and epithelial niche components, including B-cell APC, macrophage/APC, activated DC-like, pDC-like, cTEC-like, mTEC-like, thymocyte/T-cell, endothelial, fibroblast and perivascular states. Spatial projection across eight Heimli Visium sections identified recurrent MHC-II/APC-rich regions. These regions co-patterned most strongly with B-cell APC signatures, followed by mTEC-like, macrophage/APC and pDC-like programmes. In contrast, cTEC-like and T-cell/thymocyte signatures showed relative spatial separation from the strongest MHC-II/APC fields.

**Conclusions:** Public single-cell and spatial transcriptomic datasets support a multi-component model of human thymic antigen-presenting niches in which epithelial and haematopoietic APC-associated programmes are spatially organised rather than uniformly distributed. The strongest spatial signal identified an MHC-II/APC-rich axis associated with B-cell APC, mTEC-like, macrophage/APC and pDC-like signatures. These findings provide a reproducible exploratory framework for studying human thymic stromal–immune antigen-presenting architecture and generate hypotheses for future work linking thymic niche organisation to T-cell education and early-life immune competence.

## Introduction

The thymus is the central organ for T-cell development, selection and immune repertoire formation. Within the thymic microenvironment, developing thymocytes interact with epithelial, stromal and haematopoietic antigen-presenting cells that collectively shape positive selection, negative selection and central tolerance. Classical models emphasise cortical thymic epithelial cells as key mediators of positive selection and medullary thymic epithelial cells as central organisers of negative selection and tissue-restricted antigen presentation. However, human thymic antigen presentation is not mediated by epithelial cells alone. Dendritic cells, macrophages, B cells and other immune antigen-presenting populations also occupy the thymic environment and may contribute to the local organisation of tolerance and T-cell education.

Understanding how these epithelial, stromal and immune antigen-presenting compartments are organised in human thymus is particularly important in early life. The paediatric thymus is highly active, supports rapid generation of the T-cell repertoire and may influence later immune competence. Yet human thymic tissue is difficult to obtain, and spatially resolved studies remain relatively limited compared with transcriptomic profiling of dissociated cells. As a result, many questions remain about how thymic epithelial, stromal and immune antigen-presenting programmes are arranged within tissue and whether these programmes form spatially patterned niches rather than isolated cell states.

Single-cell RNA sequencing has substantially advanced the resolution of human thymic cell-state maps. These studies have identified diverse thymic epithelial, fibroblast, vascular, perivascular, thymocyte and immune populations across developmental and postnatal contexts. CITE-seq and related multimodal approaches further enable immune and stromal populations to be interpreted with both transcriptomic and protein-marker support. However, dissociated single-cell data alone cannot resolve how inferred cell states are positioned relative to one another in tissue. Spatial transcriptomic methods such as Visium provide complementary information by projecting gene-expression programmes back into tissue sections, although spot-level resolution requires careful interpretation because each spot may contain multiple cells.

Here, we used complementary public human thymus single-cell and spatial transcriptomic datasets to reconstruct stromal–immune antigen-presenting niche architecture. Rather than forcing all biological claims from a single dataset, we applied a two-arm analytical strategy. The Bautista human thymic stromal single-cell dataset was used as a reference for thymic epithelial and stromal architecture across fetal, postnatal and adult samples. The Heimli paediatric thymus CITE-seq and Visium dataset was used to resolve immune antigen-presenting, thymic B-cell, epithelial, stromal and spatially projected niche programmes. This division allowed each dataset to contribute according to its strongest experimental design.

Using this framework, we identify broad stromal and epithelial architecture in the Bautista dataset and resolve paediatric immune antigen-presenting compartments in the Heimli dataset, including B-cell APC, macrophage/APC, activated DC-like and pDC-like states. We then project single-cell-derived niche signatures into Heimli Visium thymus sections and show that MHC-II-rich spatial regions co-pattern most strongly with B-cell APC, mTEC-like, macrophage/APC and pDC-like programmes, while cTEC-like and thymocyte signatures show relative spatial separation from the strongest MHC-II/APC fields. These findings support a model in which human thymic antigen-presenting niches are spatially patterned, multi-component systems involving both epithelial and haematopoietic APC-associated programmes.

## Results

### Complementary thymus single-cell datasets support a two-arm strategy for reconstructing stromal and immune antigen-presenting niches

We first established a public-data framework for reconstructing human thymic antigen-presenting niches from complementary single-cell and spatial resources. Two datasets were prioritised because they captured distinct but biologically connected components of thymic niche architecture. The Bautista human thymus stromal single-cell dataset provided age-spanning stromal and epithelial information and was therefore treated as the primary reference for thymic epithelial and stromal architecture. The Heimli paediatric thymus CITE-seq dataset provided APC-enriched, CD45-depleted and unenriched thymic fractions, together with matched Visium spatial sections, and was therefore treated as the primary resource for immune APC, thymic B-cell and spatial niche analyses.

Exploratory UMAP analysis of the Bautista dataset resolved broad TEC, fibroblast, vascular, perivascular, vascular-MHCII-associated, cycling, capsule-associated stromal and minor contaminant compartments. Working-label annotation further separated these broad classes into biologically interpretable epithelial and stromal states, including fetal CCL25/KRT-high TEC, adult KRT5/KRT19 mTEC-like cells, immature/progenitor-like epithelial TEC, endothelial and vascular-associated stromal states, fibroblast subsets, perivascular smooth-muscle-like states and capsule-associated stromal populations. Marker-score projections supported strong epithelial, stromal and MHC-associated programmes, whereas lineage-specific immune APC and B-cell signatures were weak or inconsistent. These results confirmed that Bautista provides a robust stromal/TEC reference but should not be used as the primary source for immune APC or thymic B-cell claims.

In contrast, analysis of the Heimli paediatric thymus CITE-seq dataset recovered a broader mixture of immune, epithelial, stromal and thymocyte-associated compartments. Broad-class annotation identified B-cell APC, macrophage/APC, activated DC-like, pDC-like, TEC_cTEC-like, TEC_mTEC-like, fibroblast, endothelial, perivascular, thymocyte, T-cell and mixed low-confidence states. Fraction-level projections were internally coherent: CD45-depleted libraries contributed epithelial, stromal, endothelial, perivascular and TEC-like states; unenriched libraries were dominated by thymocyte and T-cell states; and APC-enriched libraries contributed immune/APC-associated and mixed ADT-high populations. This fraction-specific structure supported the biological plausibility of the annotation scheme and established Heimli as the primary dataset for immune APC niche analysis.

### Heimli CITE-seq resolves B-cell, myeloid and plasmacytoid antigen-presenting compartments within paediatric thymus

RNA-only marker ranking and ADT-supported annotation of the Heimli dataset identified discrete immune antigen-presenting populations within the paediatric thymus. A B-cell APC cluster was marked by MS4A1, CD79A, CD74 and multiple HLA class II genes, supporting a thymic B-cell antigen-presenting phenotype. A macrophage/monocyte-like APC cluster was marked by LYZ, CST3, CD74 and HLA-DRA/HLA-DPA1/HLA-DPB1 expression. An activated DC-like cluster was marked by CCR7, LAMP3, TNFAIP2 and BIRC3, while a pDC-like cluster was marked by LILRA4, IL3RA, GZMB, PLD4, IRF7 and JCHAIN. These populations were annotated conservatively as B-cell APC, macrophage/APC, activated DC-like and pDC-like states.

In parallel, Heimli recovered epithelial and stromal populations that connected the immune APC arm to the thymic structural framework. TEC-like clusters included a CCL25/PRSS16 cTEC-like population and KRT5/KRT19 or CALML5/S100A14 mTEC-like populations, including an antigen-presenting mTEC-like state with HLA-DRA signal. Stromal populations included fibroblast, endothelial and perivascular compartments. These compartments overlapped conceptually with the Bautista stromal architecture arm, supporting the view that Heimli captures both immune APC and epithelial/stromal niche components in paediatric thymus.

A subset of clusters was retained with low-confidence or mixed labels because their RNA and ADT patterns suggested broad immune enrichment, epithelial–MHC mixing or low lineage specificity. These mixed populations were not used as the basis for strong lineage-specific claims. Instead, interpretation focused on high-confidence or biologically coherent populations: B-cell APC, macrophage/APC, activated DC-like, pDC-like, TEC-like, endothelial, fibroblast, perivascular, thymocyte and T-cell states.

### Spatial Visium projection identifies recurrent MHC-II-rich thymic regions associated with B-cell APC, mTEC-like and myeloid APC programmes

We next projected single-cell-derived niche signatures onto eight Heimli Visium thymus sections to evaluate whether the inferred thymic APC programmes showed spatial structure. All eight sections passed first-pass spatial QC and retained sufficient tissue spots for exploratory signature analysis. Signature scoring was performed for TEC, cTEC-like, mTEC-like, MHC-II antigen-presentation, B-cell APC, macrophage/APC, activated DC-like, pDC-like, fibroblast, endothelial, perivascular and T-cell/thymocyte programmes.

At the section level, T-cell/thymocyte and MHC-II/APC signatures were consistently among the strongest mean signals, while cTEC-like programmes were recurrent across sections. Stromal signatures showed lower whole-section mean values but high local maxima, suggesting that fibroblast, endothelial and perivascular states are spatially restricted rather than uniformly distributed. These section-level patterns indicated that the tissue contains coexisting thymocyte-rich, epithelial and APC-rich regions rather than a single homogeneous antigen-presenting field.

Co-pattern analysis revealed a reproducible MHC-II/APC-centred spatial organisation. Across sections, spots with high B-cell APC and MHC-II signal were the most frequent co-high category among the analysed APC-related combinations. B-cell APC + MHC-II co-high spots represented approximately one quarter of spots on average across sections and reached the highest values in the representative S2_B1 section. mTEC-like + MHC-II, macrophage/APC + MHC-II and pDC-like + MHC-II co-high patterns were also recurrent across sections, whereas cTEC-like + MHC-II and thymocyte + APC co-high spots were comparatively infrequent.

Correlation analysis supported the co-high results. B-cell APC and MHC-II/APC signatures showed the strongest positive spatial association across all sections. mTEC-like, macrophage/APC and pDC-like signatures also showed positive association with MHC-II/APC, whereas cTEC-like and T-cell/thymocyte signatures were negatively associated with the MHC-II/APC axis. This pattern indicates that MHC-II-rich thymic regions are not simply areas of high cellularity or thymocyte density. Instead, they preferentially co-pattern with B-cell APC, mTEC-like and myeloid/pDC-like antigen-presenting programmes while being relatively separated from the strongest cTEC-like and thymocyte-associated regions.

### Integrated single-cell and spatial evidence supports a multi-component model of thymic antigen-presenting niches

The combined analysis supports a model in which human thymic antigen-presenting niches are distributed across epithelial, stromal and haematopoietic compartments. Bautista provides the stromal and TEC architecture arm, resolving epithelial, fibroblast, vascular and perivascular states across developmental and postnatal thymus. Heimli provides the paediatric immune/APC and spatial arm, resolving B-cell APC, macrophage/APC, activated DC-like, pDC-like, TEC-like, endothelial, fibroblast, perivascular and thymocyte-associated states, and projecting these programmes into spatial tissue context.

The strongest integrated spatial signal is an MHC-II-centred antigen-presenting axis that co-patterns most strongly with B-cell APC signatures, followed by mTEC-like, macrophage/APC and pDC-like programmes. In contrast, cTEC-like and thymocyte-associated signatures show relative spatial separation from the strongest MHC-II/APC fields. This suggests that thymic antigen presentation is not organised as a single homogeneous programme, but instead as a spatially patterned system involving multiple APC-associated compartments.

## Discussion

This study reconstructs human thymic stromal–immune antigen-presenting niche architecture using complementary public single-cell and spatial transcriptomic datasets. Rather than forcing all biological questions onto a single dataset, we used a two-arm analytical strategy in which Bautista defined the stromal and thymic epithelial framework, while Heimli provided the paediatric immune/APC, thymic B-cell and spatially linked niche evidence. This division proved important because the datasets differed in sampling strategy, enrichment design and biological emphasis. Bautista resolved a rich stromal/TEC landscape but had limited lineage-specific immune APC and B-cell signal, whereas Heimli captured immune APC, thymocyte, epithelial and stromal compartments and could be projected into tissue space through matched Visium sections. The resulting framework supports a multi-component model in which human thymic antigen-presenting niches are distributed across epithelial, stromal and haematopoietic APC-associated programmes.

A major finding is that paediatric thymic immune APC programmes can be resolved within the Heimli dataset as transcriptionally distinct B-cell APC, macrophage/APC, activated DC-like and pDC-like states. The B-cell APC compartment was marked by canonical B-cell and antigen-presentation features, including MS4A1, CD79A, CD74 and HLA class II genes. Myeloid APC populations showed LYZ, CST3 and MHC-II-associated expression, while pDC-like states expressed LILRA4, IL3RA, GZMB, PLD4, IRF7 and JCHAIN. These findings are consistent with the thymus containing not only epithelial antigen-presenting cells but also multiple haematopoietic APC populations capable of contributing to the antigen-presentation landscape. Importantly, we retained mixed or low-confidence clusters with conservative labels and avoided using these regions as the basis for strong lineage-specific claims.

The spatial analysis provides the strongest integrated support for a patterned antigen-presenting niche model. Across Heimli Visium sections, MHC-II/APC-rich spatial regions co-patterned most strongly with B-cell APC signatures, followed by mTEC-like, macrophage/APC and pDC-like programmes. In contrast, cTEC-like and T-cell/thymocyte signatures were negatively associated with the strongest MHC-II/APC fields. This suggests that the MHC-II-rich thymic landscape is not simply a reflection of total cellularity or thymocyte density. Instead, it appears to represent spatially organised antigen-presenting regions in which epithelial and haematopoietic APC-associated programmes converge.

The thymic B-cell/APC signal is particularly notable. In the single-cell analysis, a discrete B-cell APC cluster was recovered with strong antigen-presentation features. In the spatial analysis, B-cell APC + MHC-II was the most recurrent co-high pattern and showed the strongest positive correlation with MHC-II/APC across sections. This supports the hypothesis that thymic B cells may contribute meaningfully to local antigen-presentation niches rather than representing only rare bystander immune cells. However, because Visium spots are multicellular and the current analysis used signature scoring rather than spatial deconvolution, the result should be interpreted as co-patterning of B-cell APC and MHC-II programmes within tissue regions, not proof of direct antigen presentation by individual B cells in situ.

The mTEC-like association with MHC-II/APC also supports the biological plausibility of the spatial result. mTEC-like signatures were positively correlated with MHC-II/APC signatures across sections and were among the most frequent MHC-II co-high patterns. This is consistent with the known role of medullary epithelial antigen-presentation programmes in thymic tolerance, while also suggesting that immune APC signatures co-occupy or neighbour these antigen-presenting regions. In contrast, cTEC-like signatures were negatively associated with MHC-II/APC in the spatial analysis, suggesting that cTEC-like regions are spatially distinguishable from the strongest MHC-II-rich fields detected in these sections.

Several limitations should be considered. First, the single-cell analyses were performed on balanced downsampled objects for memory-safe exploratory processing. Therefore, composition panels should be interpreted as relative analytical summaries rather than absolute tissue abundance estimates. Second, cluster annotations are based on marker scoring, RNA-only ranked genes, ADT support and manual interpretation, and therefore remain exploratory. Third, the spatial analysis used Visium spot-level signature scoring, which cannot resolve single-cell identities within each spot and cannot establish direct cell–cell contact. Fourth, formal spatial deconvolution, neighbourhood modelling, histology-guided region annotation and external validation remain necessary before making definitive anatomical microdomain claims. Finally, although the broader motivation of this work relates to early-life immune education and vaccine responsiveness, the present analysis does not directly test vaccine-response outcomes. Such links should be treated as hypotheses for future studies.

In summary, this study identifies a spatially patterned MHC-II/APC-rich thymic antigen-presenting axis in paediatric human thymus and places it within a broader stromal–epithelial framework derived from complementary single-cell resources. The findings support a model of human thymic antigen-presenting niches as distributed, multi-cellular and spatially organised systems rather than homogeneous epithelial or immune compartments. This provides a foundation for future studies linking thymic tissue organisation to T-cell education and downstream immune responsiveness.

## Methods

### Study design and analytical overview

This study used publicly available human thymus single-cell and spatial transcriptomic datasets to reconstruct stromal, epithelial and immune antigen-presenting niche programmes in human thymic tissue. The analysis was designed as a complementary multi-dataset framework rather than a single integrated atlas forced across all samples. The Bautista human thymic stromal single-cell dataset was used as the primary reference for thymic epithelial and stromal architecture, while the Heimli paediatric thymus CITE-seq and Visium dataset was used as the primary resource for immune APC, thymic B-cell, TEC-like, stromal and spatial niche analyses.

The workflow consisted of dataset retrieval, matrix and marker audit, per-sample QC, memory-safe downsampling, single-cell normalisation, dimensionality reduction and clustering, RNA-only and ADT-separated marker interpretation, manual cluster annotation, spatial Visium section processing, spot-level signature scoring and spatial co-pattern analysis.

### Public datasets

The Bautista dataset, GSE147520, contains human thymic stromal single-cell transcriptomes across fetal, postnatal and adult samples. Seven 10x Genomics-style matrices were analysed, corresponding to fetal thymus samples at 19 and 23 weeks, postnatal samples at 6 days and 10 months, and an adult sample at 25 years. This dataset was used primarily for stromal and thymic epithelial compartment reconstruction.

The Heimli dataset, GSE207206, contains paediatric thymus CITE-seq and Visium spatial transcriptomic data. Single-cell libraries included APC-enriched, CD45-depleted and unenriched thymic fractions from five donors. PBMC libraries were retained as available raw resources but were not included in the primary thymic niche analysis. Eight Visium spatial sections were analysed: S1_A1, S1_B1, S1_C1, S1_D1, S2_A1, S2_B1, S2_C1 and S2_D1.

### Single-cell preprocessing and annotation

Raw matrix files were downloaded from GEO and inspected to confirm that matrix dimensions matched corresponding barcode and feature files. Marker-presence audits were performed before downstream analysis to confirm that genes required for epithelial, stromal, immune APC, B-cell, myeloid, pDC and antigen-presentation signatures were present.

Bautista matrices were loaded into AnnData objects, merged and filtered using at least 500 total counts, at least 200 detected genes and mitochondrial percentage below 20%. Genes detected in fewer than three cells were removed. The filtered Bautista object retained 71,488 cells. For memory-safe exploratory analysis, a balanced downsampled object was generated by selecting up to 1,500 cells per sample, producing a 9,145-cell analytical object.

Heimli thymic single-cell libraries were processed per sample. APC-enriched, CD45-depleted and unenriched thymic fractions were retained for primary analysis; PBMC libraries were excluded from the main thymic niche analysis. Cells were filtered using at least 500 total counts, at least 200 detected genes and mitochondrial percentage below 25%. Across 15 thymic libraries, 230,550 filtered cells were retained. A balanced exploratory object was generated by selecting up to 1,200 cells from each thymic library, producing an 18,000-cell object.

Downsampled objects were normalised to 10,000 counts per cell, log-transformed, embedded using PCA and UMAP, and clustered using Leiden clustering. Marker programmes were scored for TEC, cTEC, mTEC, MHC-II antigen-presentation, cDC, pDC, B-cell, macrophage/APC, T-cell/thymocyte, fibroblast, endothelial and cycling signatures. For Heimli, gene-expression and antibody-capture features were separated, and RNA-only marker ranking was used as the primary basis for cluster annotation. ADT/protein summaries were used as supportive evidence.

### Spatial preprocessing and signature scoring

Eight Heimli Visium spatial sections were processed from raw matrix, barcode, feature and spatial coordinate files. Because the files were provided in a flat GEO layout, matrices, features, barcodes and spatial position files were paired manually by section prefix. Only in-tissue spots were retained. Spots were filtered using at least 100 total counts, at least 50 detected genes and mitochondrial percentage below 30%.

Spatial signatures were defined from the single-cell analyses and included TEC, cTEC-like, mTEC-like, MHC-II/APC, B-cell APC, macrophage/APC, activated DC-like, pDC-like, fibroblast, endothelial, perivascular and T-cell/thymocyte programmes. Each spatial section was normalised to 10,000 counts per spot and log-transformed. Signature scores were calculated for each spot using available genes from each signature set.

### Spatial co-pattern analysis

Spot-level signature scores were z-scored within each section. Four broad spatial axes were calculated: epithelial/TEC, immune APC, structural stromal and T-cell/thymocyte. Dominant spatial axes and dominant signatures were assigned using the highest z-scored axis or signature when above a threshold of 0.5; spots below this threshold were labelled low/mixed.

Co-high analyses identified spots in which two signatures were simultaneously high. Co-high rules included cTEC + MHC-II, mTEC + MHC-II, B-cell APC + MHC-II, macrophage/APC + MHC-II, pDC + MHC-II, fibroblast + TEC, endothelial + APC and thymocyte + APC. For each section, co-high spot counts and percentages were calculated. Spearman correlations were also calculated between selected spatial signatures within each section.

### Software

Analyses were performed in Python using Scanpy, AnnData, pandas, NumPy, SciPy and Matplotlib in a WSL/Linux environment. Sparse-matrix-aware operations and balanced downsampling were used to avoid memory-intensive dense matrix conversion. Exact scripts and processed outputs are retained in the project repository.

## Figure legends

**Figure 1. Complementary single-cell datasets define stromal/TEC and paediatric immune-niche arms for human thymic antigen-presenting niche reconstruction.**  
**A**, UMAP of the Bautista human thymic stromal single-cell dataset annotated by broad stromal and epithelial classes, including TEC, fibroblast, vascular, vascular-MHCII-associated, perivascular, capsule-associated stromal, cycling, immune-contaminant and contaminant compartments. **B**, Relative broad-class composition across Bautista samples ordered by developmental stage. Composition values are calculated within the balanced downsampled analytical object and are not intended as absolute tissue abundance estimates. **C**, Bautista UMAP coloured by sample, showing the contribution of fetal, postnatal and adult thymus samples to the exploratory stromal/TEC atlas. **D**, UMAP of the Heimli paediatric thymus CITE-seq dataset collapsed into simplified niche classes, including immune APC, TEC/epithelial, stromal/vascular, T-cell/thymocyte, mixed/low-confidence and contaminant groups. **E**, Relative simplified-class composition across Heimli APC-enriched, CD45-depleted and unenriched fractions. Composition values are calculated within the balanced downsampled analytical object. **F**, Heimli UMAP coloured by fraction, showing fraction-specific structure across APC-enriched, CD45-depleted and unenriched thymic libraries. Together, these analyses define Bautista as the primary stromal/TEC architecture reference and Heimli as the primary paediatric immune/APC and spatially linked thymic niche dataset.

**Figure 2. Heimli single-cell marker programmes support immune APC, epithelial and annotation-confidence structure in paediatric thymus.**  
**A**, Detailed broad-class UMAP annotation of the Heimli paediatric thymus CITE-seq dataset, resolving B-cell APC, macrophage/APC, activated DC-like, pDC/APC, TEC-like, epithelial-mixed, stromal, endothelial, perivascular, thymocyte, T-cell and mixed low-confidence populations. **B–E**, UMAP projections of marker-programme scores for MHC-II/APC (**B**), B-cell APC (**C**), macrophage/APC (**D**) and pDC/APC (**E**) signatures. These scores support distinct immune antigen-presenting programmes within the paediatric thymic dataset. **F–G**, UMAP projections of cTEC-like (**F**) and mTEC-like (**G**) scores, supporting epithelial programmes that complement the immune APC compartments. **H**, Annotation-confidence map showing high-, medium- and low-confidence regions of the exploratory annotation. Low-confidence and mixed clusters were retained but interpreted conservatively. Together, the marker programmes support a multi-component thymic antigen-presenting system involving B-cell, myeloid, pDC-like and epithelial-associated compartments.

**Figure 3. Spatial co-patterning of MHC-II-rich antigen-presenting programmes in Heimli Visium thymus sections.**  
**A–E**, Representative spatial maps from Heimli Visium section S2_B1 showing within-section z-scored signature values for MHC-II/APC (**A**), B-cell APC (**B**), mTEC-like (**C**), cTEC-like (**D**) and T-cell/thymocyte (**E**) programmes. Spatial maps use a shared z-score colour scale from −3 to +3. S2_B1 was selected as a representative section because it showed strong MHC-II/APC signal and the highest B-cell APC + MHC-II co-high fraction among the analysed sections. **F**, Quantification of recurrent co-high spatial patterns across all eight Heimli Visium sections. Bars show the mean percentage of co-high spots across sections, and error bars indicate the observed section-level minimum and maximum. B-cell APC + MHC-II was the most frequent co-high pattern, followed by mTEC-like + MHC-II, macrophage/APC + MHC-II and pDC-like + MHC-II. **G**, Heatmap of section-level Spearman correlations between selected spatial signatures. MHC-II/APC showed strong positive spatial association with B-cell APC, mTEC-like, macrophage/APC and pDC-like signatures, while cTEC-like and T-cell/thymocyte signatures showed negative association with MHC-II/APC. These results support a spatially patterned MHC-II/APC-rich thymic antigen-presenting axis involving both epithelial and haematopoietic APC-associated programmes. The spatial analysis is interpreted as signature co-patterning rather than formal deconvolution or direct cell–cell contact inference.

## Data and code availability

All analyses were performed using publicly available datasets from GEO. Processed analytical objects, scripts, tables and figure files are organised within the local project directory and will be deposited in an appropriate public repository upon manuscript submission. The repository will include scripts for data retrieval, QC, single-cell processing, CITE-seq RNA/ADT marker separation, spatial signature scoring, co-pattern analysis and figure generation.

## Author contributions

R.S.M. conceived the study, designed the computational workflow, performed data analysis, interpreted the results and drafted the manuscript. Additional author contributions will be added according to project participation.

## Acknowledgements

The authors acknowledge the investigators who generated and made publicly available the human thymus single-cell and spatial transcriptomic datasets re-analysed in this study.

## Competing interests

The author declares no competing interests.

## References

References will be added and formatted according to the target journal style.
"""

manuscript = textwrap.dedent(manuscript).strip() + "\n"

out_md = DOC_DIR / "manuscript_draft_v1.md"
out_md.write_text(manuscript)

assets = [
    {
        "figure": "Figure 1",
        "status": "frozen_near_final",
        "png": str(FIG_DIR / "figure1_final_dataset_architecture.png"),
        "svg": str(FIG_DIR / "figure1_final_dataset_architecture.svg"),
        "pdf": str(FIG_DIR / "figure1_final_dataset_architecture.pdf"),
    },
    {
        "figure": "Figure 2",
        "status": "frozen_near_final",
        "png": str(FIG_DIR / "figure2_final_marker_programmes.png"),
        "svg": str(FIG_DIR / "figure2_final_marker_programmes.svg"),
        "pdf": str(FIG_DIR / "figure2_final_marker_programmes.pdf"),
    },
    {
        "figure": "Figure 3",
        "status": "frozen_near_final",
        "png": str(FIG_DIR / "figure3_final_spatial_copatterns.png"),
        "svg": str(FIG_DIR / "figure3_final_spatial_copatterns.svg"),
        "pdf": str(FIG_DIR / "figure3_final_spatial_copatterns.pdf"),
    },
]

asset_df = pd.DataFrame(assets)
asset_out = DOC_DIR / "manuscript_figure_assets_v1.tsv"
asset_df.to_csv(asset_out, sep="\t", index=False)

print("Done.")
print(f"Wrote manuscript draft: {out_md}")
print(f"Wrote figure asset table: {asset_out}")
print()
print("Figure asset status:")
print(asset_df.to_string(index=False))
