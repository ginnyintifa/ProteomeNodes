# ProteomeNodes
#1 Introduction 
This package achieves curve fitting and classification for proteome and phosphoproteome stimulated with different treatment dosage and under two different temperatures. 

#2 Installation and preparation 
ProteomeNodes can be downloaded and installed in R. Installation of GPD requires devtools as a prerequisite:

```{r}
install.packages("devtools")
```
Next, install ProteomeNodes by:

```{r}
library("devtools")
devtools::install_github("ginnyintifa/ProteomeNodes")
library(ProteomeNodes)
```
Input files needed are output from maxQuant

* s1_phospho_file
* s2_phospho_file
* s1_peptide_file
* s2_peptide_file


We also need the kinase-substrate relationship file which is provided in this repo. 

* ks_network_file

#3 Function 
To run the function we call

```{r}
dosage_curves_2T(
  s1_phospho_file,
  s2_phospho_file,
  s1_peptide_file,
  s2_peptide_file,
  ks_network_file,
  prob_threshold,
  dosages,
  range_threshold,
  s1_tag,
  s2_tag,
  s1_rep1_colnames,
  s1_rep2_colnames,
  s1_rep1_omit,
  s1_rep2_omit,
  s2_rep1_colnames,
  s2_rep2_colnames,
  s2_rep1_omit,
  s2_rep2_omit,
  bandwidth,
  lead_phospho_BL,
  lead_peptide_BL,
  order_s,
  both_s_BL,
  output_dir)
```

Here we further explain the choice of the the following parameters.

* bandwidth: a numeric value to adjust the smoothness of the fitted curve, default to 0.1. The larger bandwidth is the smoother the curve is. 
* range_threshold: a numeric value to categorise curves to NR and non-NR, default to 0.3. This is the difference between the maximum and minimum value of the fitted curve. 
 


