## 1 Introduction 
This package achieves fold change calculation and curve fitting and classification for proteome and phosphoproteome stimulated with different treatment dosage and under two different temperatures. 


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


ProteomeNodes has two modes for two different experimental settings. If you monitor abundance changes in response to multiple doses (normally around 10 different doses), choose mode ```dosage_curves_2T```. Alternatively, if you have two treatment levels and one control conditaion, choose mode ```foldChange_2T```.

In ```dosage_curves_2T``` mode, input files needed are output from maxQuant

* s1_phospho_file
* s2_phospho_file
* s1_peptide_file
* s2_peptide_file


We also need the kinase-substrate relationship file which is provided in this repo. 

* ks_network_file

In ```foldChange_2T``` mode, in addition to the raw maxquant output files, it is necessary to provide a input specification file, such as one in "spCol_positions.txt"

```
condition	position
sp1ctrlCol	16
sp2ctrlCol	19
sp3ctrlCol	22
sp1lowCol	17
sp2lowCol	20
sp3lowCol	23
sp1highCol	18
sp2highCol	21
sp3highCol	24
```

Here, column "position" is the column number of the specified conditions in both input data files. 




## 4 Function 


Here are examples of running the two functions:

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
 
 
```{r}


foldChange_2T(
rawData_inputfilename1 = "MEHP37_F_PeptideGroups.txt",
rawData_inputfilename2 = "MEHP53_F_PeptideGroups.txt",
inputSpec_filename =  "spCol_positions.txt",
rep_number = 3,
workingDir =  "your/output/dir/",
mapDIA_flag = T,
mapDIADir = "dir/to/mapDIA-master/",
userSystem = "Linux",
defaultInput_parameter = "mapDIA/inputfile/name",
cutoff_uFc = 1.3,
cutoff_lFc = 0.7,
cutoff_mapDIA = 0.05,
outputTagS1 = "lowVSctrl",
outputTagS2 = "highVSctrl")


```
Note: mapDIA is a software used for generating statistical significance for small-sample-size comparisons. It can be downloaded from: https://sourceforge.net/projects/mapdia/files/

* rep_number: the number of replicates under each experimental condition. 
* mapDIADir: the directory where mapDIA is installed.
* defaultInput_parameter: the absolute path where mapDIA input file generated by the function will be sent to.  


## 5 Output files

In ```dosage_curves_2T``` mode, these are the final output files:


* both_ProtPhospho_curve.pdf: the figure of curves fitted for proteins and phosphorylation sites. Plots in the first and third column are protein changes at 37C and 52C respectively, whereas plots in the second and fourth column are one of phosphosite changes of corresponding proteins at 37C and 52C respectively. Blue and red colors contrast the two replicates.  

<img src="https://github.com/ginnyintifa/ProteomeNodes/blob/master/both_ProtPhospho_curve.png" align="center"/>

* not_both_ProtPhospho_curve.pdf: in the same format as the figure above, but this file has only protein or phosphosite curves. 


<img src="https://github.com/ginnyintifa/ProteomeNodes/blob/master/prot_2s.png" align = "center"/>

* prot_2s.pdf: the scatter plots of average slopes of proteins in two temperatures. Proteins beyond specified thresholds are colored and labeled. 


Detailed results for proteins and phosphosites can be found in the following .tsv files:

* 37C_protFitted.tsv
* 52C_protFitted.tsv
* 37C_phosphoFitted.tsv
* 52C_phosphoFitted.tsv
* prot_2s.tsv

In ```foldChange_2T``` mode, these are the final output files:

<img src="https://github.com/ginnyintifa/ProteomeNodes/blob/master/highVSctrl.png" align="center"/>

* highVSctrl.pdf: scatter plot of protein fold changes (high condition over control) in two temperatures. Proteins reliable and significant are colored and labeled. 
* lowVSctrl.pdf: scatter plot of protein fold changes (low condition over control) in two temperatures.

Detailed results corresponding to the two figures can be found in:

* highVSctrl.tsv
* lowVSctrl.tsv
 



 







