# Isoniazid-PopPK-Repository
This repository contains a collection of population pharmacokinetic (PopPK) models for Isoniazid (INH), integrating multiple published studies to simulate pharmacokinetic behavior under different physiological conditions and genetic phenotypes.

---
## Project Structure
- **Main Script**: `INH PopPK models.R`  
- **Output Folder**: `INH0001/` (automatically generated, contains model prediction plots)

---
## Dependencies

### R packages
Ensure the following R packages are installed before running the script:
install.packages(c("tidyverse", "mrgsolve", "rxode2", "cowplot", "ggplot2", 
                   "ggpubr", "dplyr", "PopED", "ggsci", "gridExtra"))

## Contributors
Author: Graham Ju
Affiliation: Xiangya Hospital, Central South University, Hunan Province, China
Email: jugehang@163.com or 218101055@csu.edu.cn
Date: 2025/02/23

## License
This project is open-source under the MIT License. For details, refer to the script comments.

## Acknowledgments
This repository is based on published pharmacokinetic studies. The authors of the original studies are acknowledged for their contributions to the field of pharmacokinetics.
