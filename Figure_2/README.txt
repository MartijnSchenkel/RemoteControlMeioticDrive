./analysis_script.R 
This file runs the analysis required to generate Figure 2. This script is divided into two sections,
with the first generating the datasets 2A, 2B, and 2C. These track allele frequency changes across 
generations. The second section of the analysis script generates datasets 2D and 2E, which test a 
wide parameter space. 2D and 2E therefore take much longer to run. The analysis script writes all 
five datasets to the working directory (which should be set to ./data or later moved to ./data).

./figure_script.R
This file generates Figure 2 from the datasets present in ./data and writes both a .pdf and a .png 
version of the figure, labeled 2024_12_30_Figure_2.png and 2024_12_30_Figure_2.pdf to the working 
directory (which should be set to ./figure_files or later moved to ./figure_files). 
