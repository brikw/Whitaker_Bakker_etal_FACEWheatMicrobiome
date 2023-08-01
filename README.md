# "Manipulating atmospheric CO2 concentration induces shift in wheat leaf and spike microbiome and Fusarium pathogen communities over time"
### Bakker, M.G., Whitaker, B.K., McCormick, S.P., Ainsworth, E.A., Vaughan, M.M.

DOI HERE

This repository includes the R code, data files, small scripts, and a meta file to supplement the manuscript by Bakker & Whitaker et al. "Manipulating atmospheric CO2 concentration induces shift in wheat leaf and spike microbiome and Fusarium pathogen communities over time".

The field experiment was conducted as part of SOY-FACE at UIUC, with ambient and elevated CO2 treatments. 

The two files (“2017 FACE field data Biomass and DON.csv” and “2018 Data Biomass and Toxin.csv”) contain collected data and treatments information for the processing of Fusarium inoculated wheat spike samples (different samples than those processed for microbiome). Fungal biomass was constructed in 2017 from a single gene of Fusarium and single gene of wheat, using Bio_Rad, and expressed as a ratio. Fungal biomass in 2018 was constructed from multiple genes using a Fluidigm, also expressed as a ratio.

The metadata file ("Metadata.csv") includes information about each microbiome sample processed using Illumina MiSeq (different samples than those inoculated with Fusarium). The microbiome data was processed using Illumina MiSeq and DADA2. The output of DADA2 processing consisted of phyloseq objects for each year (2017 and 2018) and microbial kingdom (bacterial = b, and fungal = f). For convenience, filtered and unfiltered versions of the phyloseq objects are provided.

Detailed information about column headers in the Biomass/DON and Microbiome Metadata files can be found in the .xlsx file "FACE_Meta.xlsx".

The R/Rmd files should be run in this order: FACE-Toxin-Biomass.Rmd > FACE_sequence_processing.R > FACE-Microbes-CommDiv.Rmd > FACE-Microbes-DA.Rmd. 

The /code folder contains the batch and R scripts necessary to run small codes inside R. Data necessary to run the analyses can be found in /data folder.

Please see the manuscript for details and full reference information.
