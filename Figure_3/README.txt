./analysis_script.R 
This file runs the analysis required to generate Figure 3. This script is divided into two sections,
with the first generating the datasets 3A, 3B, and 3C. These track allele frequency changes across 
generations. The second section of the analysis script generates datasets 3D and 3E, which test a 
wide parameter space. 2D and 2E therefore take much longer to run. The analysis script writes all 
five datasets to the working directory (which should be set to ./data or later moved to ./data).

./figure_script.R
This file generates Figure 3 from the datasets present in ./data and writes both a .pdf and a .png 
version of the figure, labeled 2024_12_30_Figure_3.png and 2024_12_30_Figure_3.pdf to the working 
directory (which should be set to ./figure_files or later moved to ./figure_files). 
