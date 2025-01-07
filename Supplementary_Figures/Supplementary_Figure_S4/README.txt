./analysis_script.R 
This file runs the analysis required to generate Supplementary Figure S4. This script is divided 
into two sections, with the first generating the datasets S4A, S4B, and S4C. These track allele 
frequency changes across generations. The second section of the analysis script generates datasets 
S4D and S4E, which test a wide parameter space. S4D and S4E therefore take much longer to run. The 
analysis script writes all five datasets to the working directory (which should be set to ./data or 
later moved to ./data).

./figure_script.R
This file generates Supplementary Figure S4 from the datasets present in ./data and writes both a 
.pdf and a .png version of the figure, labeled 2024_12_30_Figure_S4.png and 2024_12_30_Figure_S4.pdf 
to the working directory (which should be set to ./figure_files or later moved to ./figure_files). 
