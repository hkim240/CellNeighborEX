# CellContact
CellContact is a statistical framework to detect cell contact-dependent gene expression. This is a new framework to analyze spatial transcriptomics (ST) data using spatial beads composed of muliple cell types, especially for Slide-seq data with 10 µm resolution. 

The figure below shows the workflow for detecting cell contact-dependent gene expression.

<img width="1426" alt="Fig 1a" src="https://user-images.githubusercontent.com/99720939/156026752-0f3dd260-3c00-48fb-974b-76d89e7b22ea.png">

(1) Spatial transcriptomics data: Slide-seq V2 data are used. <br>
(2) Cell type deconvolution: The cell type mixture of the individual beads in Slide-seq are decomposed by RCTD. <br>
(3) Cell contact-dependent gene expresson analysis: DEGs between the heterotypic and homotypic beads are identified by CellContact. <br>
(4) Spatial validation: The expression of the DEGs are visualized via spatial mapping.<br>

The figure below shows how to create a null model of aritificial heterotypic beads.

![Suppl_Fig 1](https://user-images.githubusercontent.com/99720939/156159727-6fe693aa-ce23-4760-9506-8d6bdde4cb37.png)

The statistical significance of the DEGs is confirmed through the null model. The process is implemented in CellContact. Specifically, CellContact performs four steps as follows:

STEP1: Converting log-normalized data to z-values <br>
STEP2: Creating the null model of artificial heterotypic spots <br>
STEP3: Finding differentially expressed genes (DEGs) <br>
SETP4: Plotting heatmaps <br> 


# Citation
Hyobin Kim, Cecilia Lövkvist, Patrick C.N. Martin, Junil Kim, Kyoung Jae Won, Detecting Cell Contact-dependent Gene Expression from Spatial Transcriptomics Data, bioRxiv, 2022.
