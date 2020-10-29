## 1 Introduction 
This package achieves curve fitting and classification for proteome and phosphoproteome stimulated with different treatment dosage and under two different temperatures. 

## 2  Installation 

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


## 3 Input files

Input files needed are output from maxQuant

* s1_phospho_file
* s2_phospho_file
* s1_peptide_file
* s2_peptide_file


We also need the kinase-substrate relationship file which is provided in this repo. 

* ks_network_file

## 4 Function 
To run the function we call function 

```{r}
dosage_curves_2T
```
Here shows an example of running the function with the sample file we provide.


```{r}
dosage_curves_2T(
  s1_phospho_file = "s1_phospho.tsv",
  s2_phospho_file = "s2_phospho.tsv",
  s1_peptide_file = "s1_peptide.tsv",
  s2_peptide_file = "s2_peptide.tsv",
  ks_network_file = "mergeNetwork.tsv",
  prob_threshold = 0.5,
  dosages = 10,
  range_threshold = 0.3,
  s1_tag = "37C",
  s2_tag = "52C",
  s1_rep1_colnames = "_B2",
  s1_rep2_colnames = "_B3",
  s1_rep1_omit = NA,
  s1_rep2_omit = NA,
  s2_rep1_colnames = "_B2",
  s2_rep2_colnames = "_B3",
  s2_rep1_omit = 2,
  s2_rep2_omit = 2,
  bandwidth = 0.1,
  lead_phospho_BL = T,
  lead_peptide_BL = F,
  order_s = 1,
  both_s_BL = T,
  output_dir  = "your/output/dir/")
```

Here we further explain the choice of the the following parameters.

* bandwidth: a numeric value to adjust the smoothness of the fitted curve, default to 0.1. The larger bandwidth is the smoother the curve is. 
* range_threshold: a numeric value to categorise curves to NR and non-NR, default to 0.3. This is the difference between the maximum and minimum value of the fitted curve. 
 

## 5 Output files

