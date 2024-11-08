# Transcription Program Gradients
## Human Tumor Atlas Network (HTAN) Data Jamboree | Nov. 6-8, 2024
![Image showing overall](2024_overall_figure.png)
## HTAN Transcription Program Gradients Team
- **Yukun Tan**, Department of Bioinformatics and Computational Biology, University of Texas MD Anderson Cancer Center, Houston, TX
- **Mary Goldman**, Genomics Institute, UC Santa Cruz, Santa Cruz, CA
- **Minji Kim**, Department of Artificial Intelligence and Informatics Research, Mayo Clinic, Jacksonville FL
- **Irene Marin Goni**, Machine Learning for Computational Biology, CIMA University of Navarra, Spain / MPET Department, Mayo Clinic, Rochester, MN
- **Neel Sanghvi**, Cancer Data Science Laboratory, NCI (NIH), Bethesda, MD
- **Lon Fong**, PRIME-TR, University of Texas MD Anderson Cancer Center, Huston TX
- **Syed Abbas Bukhari**, Laboratory of Pediatric Oncology, NCI (NIH), Bethesda, MD
- **Phi Le**, Department of Medicine, University of California San Francisco, San Francisco, CA

## Overview
Local Spatial Gradient Inference (LSGI) is a tool that systematically identifies spatial locations with prominent, interpretable spatial transcriptomic gradients (STGs) from spatial transcriptomic (ST) data ([Manuscript](https://www.biorxiv.org/content/10.1101/2024.03.19.585725v1), [Github](https://github.com/qingnanl/LSGI)).

We applied LSGI to breast cancer data in HTAN, specifically exploring if transcriptional gradients vary between breast cancer subtypes.

We reanalyzed data from [Iglesia 2024](https://www.nature.com/articles/s43018-024-00773-6), which showed characteristic links in gene expression and chromatin accessibility between breast cancer subtypes. In particular they found exhausted CD8+ T cells were enriched in basal-like breast cancer, suggesting an altered means of immune dysfunction. We hope to verify this exhausted CD8+ T cells in the results from the LSGI algorithm.

![Image showing workflow](Workflow.png)

## Data Details
Visium 10X data were downloaded from the Human Tumor Atlas Network’s WUSTL atlas (HTAN WUSTL). We downloaded NN files from NN samples. This is the [direct link]() to the HTAN filtered list of samples with the freely accessible Seurat RDS (R code files). 

## Project Description
This project includes follwing R scripts:
1. `spatial_gradient_pipeline.R`
### Workflow for spatial_gradient_pipeline.R
1. **Data Preprocessing**
    - Load spatial transcriptomics data from a 10X Genomics dataset.
    - Create a Seurat object for spatial analysis, link to corresponding imaging data, and filter out unwanted features such as ribosomal and mitochondrial genes.
    - Normalize the data to make it suitable for downstream analysis.
2. **Non-negative Matrix Factorization (NMF)**
    - Apply NMF to identify spatial patterns in gene expression data.
    - This step helps reduce dimensionality, making it easier to visualize and interpret the data.
3. **Local Spatial Gradient Inference (LSGI)**
    - Extract spatial coordinates and NMF embeddings.
    - Perform LSGI analysis to infer local spatial gradients, which provides insights into the local relationships between cell populations.
4. **Functional Annotation**
    - Use the hallmark gene sets to functionally annotate the identified NMF factors.
    - Perform hypergeometric tests to enrich for biological pathways.
    - Visualize the results to understand the biological relevance of identified spatial modules.
5. **Visualization**
    - Generate plots for spatial gradients and perform distance heatmap analysis.
    - Create visualizations of the enriched gene sets to understand their biological significance.
### Input
- **10X Genomics Visium Data**: The input data should include the filtered feature-barcode matrix and the spatial image data for breast cancer samples.
- **Parameters for NMF**: Users can specify the rank for NMF analysis and the tolerance level for optimization.
### Output
- **Processed Seurat Object**: An RDS file containing the processed Seurat object (e.g.,`HT206B1-S1Fc1U2Z1B1_processed_data.rds`), which includes the spatial transcriptomics data with linked imaging and NMF embeddings.
- **LSGI Results**: An RDS file (e.g.,`HT206B1-S1Fc1U2Z1B1_lsgi.rds`) containing the results from the LSGI analysis.
- **Visualizations**: Plots of the spatial gradients, distance heatmaps, and enriched gene sets saved as PDFs.
- **Functional Annotations**: Excel files containing the hallmark enrichment results (e.g.,`MESO_fov19_Hallmark_hypeR.xlsx`).
### Requirements
- R version 4.0 or higher
- Required libraries: `Seurat`, `Matrix`, `RcppML`, `ggplot2`, `dplyr`, `LSGI`, `hypeR`, `msigdbr`
## Usage
1. Clone this repository and set the correct working directory.
2. Ensure that all required libraries are installed.
3. Run the script with the input data in the appropriate format.
4. Outputs will be saved in the current directory, including processed data, LSGI results, and visualization plots.
## Results

## Future features

## License
This project is licensed under the MIT License. See the [LICENSE](https://github.com/NCI-HTAN-Jamborees/Transcription-Program-Gradients/blob/main/LICENSE) file for more details.

## Acknowledgements
Thank you to the Human Tumor Atlas Network (HTAN), NIH, NCI, ISB-CGC and Cancer Genomics Cloud (Seven Bridges) for all support.

