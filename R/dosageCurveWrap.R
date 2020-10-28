
#' Fit curves 
#' @param s1_phospho_file file name of the phosphopeptide data under scenario 1 
#' @param s2_phospho_file file name of the phosphopeptide data under scenario 2
#' @param s1_peptide_file file name of the peptide data under first scenario 1 
#' @param s2_peptide_file file name of the peptide data under second scenario 2 
#' @param ks_network_file file name of the kinase substrate network
#' @param prob_threshold  probablity threshold for phosphosites filtering 
#' @param dosages number of dosages 
#' @param s1_tag a name tag for scenario 1
#' @param s2_tag a name tag for scenario 2
#' @param s1_rep1_colnames colnames of scenario 1 replicate 1 to grep relevant columns
#' @param s1_rep2_colnames colnames of scenario 1 replicate 2 to grep relevant columns
#' @param s1_rep1_omit colnames of scenario 1 replicate 1 to omit relevant columns
#' @param s1_rep2_omit colnames of scenario 1 replicate 2 to omit relevant columns
#' @param s2_rep1_colnames colnames of scenario 2 replicate 1 to grep relevant columns
#' @param s2_rep2_colnames olnames of scenario 2 replicate 2 to grep relevant columns
#' @param s2_rep1_omit colnames of scenario 2 replicate 1 to omit relevant columns
#' @param s2_rep2_omit colnames of scenario 2 replicate 2 to omit relevant columns
#' @param bandwidth bandwidth for Gaussian kernal the larger the smoothier the curve is
#' @param range_threshold a number set to be the threshold for range of the fitted curve 
#' @param lead_phospho_BL a boolean value for filthering only leading proteins or not
#' @param lead_peptide_BL a boolean value for filthering only leading proteins or not
#' @param order_s which replicate should be arranged at the first place 
#' @param both_s_BL a boolean value, show pairs with both protein and phospho or not
#' @param output_dir output directory 
#' @import dplyr data.table magrittr 
#' @export






dosage_curves_2T = function(
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


{
  
  #  
  # s1_phospho_file = "/data/ginny/YanTing20200227/InputFiles/B2B3_37/Phospho (STY)Sites.txt"
  # s2_phospho_file =  "/data/ginny/YanTing20200227/InputFiles/B2B3_52/Phospho (STY)Sites.txt"
  # s1_peptide_file = "/data/ginny/YanTing20200227/InputFiles/B2B3_37_input/peptides.txt"
  # s2_peptide_file = "/data/ginny/YanTing20200227/InputFiles/B2B3_52_input/peptides.txt"
  # prob_threshold = 0.5
  # ks_network_file = "/data/ginny/KinaseSubstrate/mergeNetwork.tsv"
  # dosages = 10
  # range_threshold = 0.3  ### for biphasic thing
  # s1_tag = "37C"
  # s2_tag = "52C"
  # s1_rep1_colnames = "_B2"
  # s1_rep2_colnames = "_B3"
  # s1_rep1_omit = NA
  # s1_rep2_omit = NA
  # s2_rep1_colnames = "_B2"
  # s2_rep2_colnames = "_B3"
  # s2_rep1_omit = 2
  # s2_rep2_omit = 2
  # bandwidth = 0.1    ### can change to 0.1 later 
  # lead_phospho_BL = T
  # lead_peptide_BL = F
  # order_s = 1
  # both_s_BL = T
  # output_dir = "/data/ginny/YanTing20200227/OutFiles/compile0617/"
  # # 
  #  #### file cleaning 
  # 
  #### in the following function I need to add a column of peptide sequence 
  
  fileCleanPhosphoSites (original_file =  s1_phospho_file,
                         probability_threshold = prob_threshold,
                         rep1Colnames = s1_rep1_colnames,
                         rep2Colnames = s1_rep2_colnames,
                         output_name = paste0(output_dir, s1_tag, "_phospho.tsv"))
  
  fileCleanPhosphoSites (original_file =  s2_phospho_file,
                         probability_threshold = prob_threshold,
                         rep1Colnames = s2_rep1_colnames,
                         rep2Colnames = s2_rep2_colnames,
                         output_name = paste0(output_dir, s2_tag, "_phospho.tsv"))
  
  fileCleanPeptides (original_file = s1_peptide_file,
                     rep1Colnames = s1_rep1_colnames,
                     rep2Colnames = s1_rep2_colnames,
                     output_name = paste0(output_dir, s1_tag, "_peptide.tsv"))
  
  fileCleanPeptides (original_file = s2_peptide_file,
                     rep1Colnames = s2_rep1_colnames,
                     rep2Colnames = s2_rep2_colnames,
                     output_name = paste0(output_dir, s2_tag, "_peptide.tsv"))
  
  cat("files cleaned.", "\n")
  
  
  ##### get double reps
  
  
  doublePhospho_s1 = filterDoublePhospho(Longfile_name = paste0(output_dir, s1_tag, "_phospho.tsv"),
                                         rep_number = dosages,
                                         leadBL = lead_phospho_BL)
  
  
  doublePhospho_s2 = filterDoublePhospho(Longfile_name = paste0(output_dir, s2_tag, "_phospho.tsv"),
                                         rep_number = dosages,
                                         leadBL = lead_phospho_BL)
  
  
  doublePep_s1 = filterDoublePeptide(pep_filename = paste0(output_dir, s1_tag, "_peptide.tsv"),
                                     rep_number = dosages,
                                     leadBL = lead_peptide_BL)
  
  doublePep_s2 = filterDoublePeptide(pep_filename = paste0(output_dir, s2_tag, "_peptide.tsv"),
                                     rep_number = dosages,
                                     leadBL = lead_peptide_BL)
  
  cat("double rep extracted.", "\n")
  
  #### summarise peptide to proteins 
  
  sumPep_s1 = sumPeptideToProtein(peptide = doublePep_s1,
                                  leadBL = lead_peptide_BL,
                                  rep1Colnames = s1_rep1_colnames,
                                  rep2Colnames = s1_rep2_colnames,
                                  output_name = paste0(output_dir, s1_tag,"_sum_peptide.tsv"))
  
  sumPep_s2 = sumPeptideToProtein(peptide = doublePep_s2,
                                  leadBL = lead_peptide_BL,
                                  rep1Colnames = s2_rep1_colnames,
                                  rep2Colnames = s2_rep2_colnames,
                                  output_name = paste0(output_dir, s2_tag,"_sum_peptide.tsv"))
  
  
  cat("peptides summed to proteins.", "\n")
  
  
  ##### fit curves 
  
  ##### need to improve on assigning biphasic first 
  #### 1 max appears in the middle  and min apears in ends 
  #### 2 the distance between min and max is large 
  
  
  
  
  double_fitPhospho_s1 = curveFittingDouble(cleanedData = doublePhospho_s1,
                                            dosages = dosages,
                                            distanceTH = range_threshold,
                                            rep1Colnames = s1_rep1_colnames,
                                            rep2Colnames = s1_rep2_colnames,
                                            rep1Omit = s1_rep1_omit,
                                            rep2Omit = s1_rep2_omit,
                                            bandwidthSet = bandwidth)
  
  ### use the new function here 
  
  #outputName = paste0(output_dir, s1_tag, "_phosphoFitted.tsv"))
  #  outputName = "/data/ginny/YanTing20200227/OutFiles/double37Lead_curveFitting_0305.tsv")
  # 
  # pho1 = double_fitPhospho_s1%>%
  #   dplyr::filter(response == "biphasic")
  # 
  # 
  # pho2 = double_fitPhospho_s1%>%
  #   dplyr::filter(response == "hypo")
  
  
  
  ### 0511 add the hgnc symbol to the id columns 
  
  double_fitPhospho_s1_output = outputTable(phosphoSiteFit = double_fitPhospho_s1,
                                            uniprot_gn_filename = "/data/ginny/tcga_pancan/important_files/sp_0814_fasta_uniprot_gn.tsv",
                                            kinaseSubstrate_filename = ks_network_file,
                                            outputName = paste0(output_dir, s1_tag, "_phosphoFitted.tsv"))
  
  
  
  
  double_fitPhospho_s2 = curveFittingDouble(cleanedData = doublePhospho_s2,
                                            dosages = dosages,
                                            distanceTH = range_threshold,
                                            rep1Colnames = s2_rep1_colnames,
                                            rep2Colnames = s2_rep2_colnames,
                                            rep1Omit = s2_rep1_omit,
                                            rep2Omit = s2_rep2_omit,
                                            bandwidthSet = bandwidth)
  
  double_fitPhospho_s2_output = outputTable(phosphoSiteFit = double_fitPhospho_s2,
                                            uniprot_gn_filename = "/data/ginny/tcga_pancan/important_files/sp_0814_fasta_uniprot_gn.tsv",
                                            kinaseSubstrate_filename = ks_network_file,
                                            outputName = paste0(output_dir, s2_tag, "_phosphoFitted.tsv"))
  
  
  
  
  ### may need to implement different threshold for proteins 
  ### becuase they are more steady 
  ### can discuss with YT 
  
  
  double_fitPep_s1 =  curveFittingDouble(cleanedData = sumPep_s1,
                                         dosages = dosages,
                                         distanceTH = range_threshold,
                                         rep1Colnames = s1_rep1_colnames,
                                         rep2Colnames = s1_rep2_colnames,
                                         rep1Omit = s1_rep1_omit,
                                         rep2Omit = s1_rep2_omit,
                                         bandwidthSet = bandwidth)
  
  #  outputName = paste0(output_dir, s1_tag,"_protFitted.tsv"))
  ### work on this after lunch 
  
  
  
  double_fitPep_s1_output = outputTable_protein(protFit = double_fitPep_s1,
                                                uniprot_gn_filename = "/data/ginny/tcga_pancan/important_files/sp_0814_fasta_uniprot_gn.tsv",
                                                outputName = paste0(output_dir, s1_tag, "_protFitted.tsv"))
  
  
  ### a new function for protein output table 
  
  
  
  double_fitPep_s2 =  curveFittingDouble(cleanedData = sumPep_s2,
                                         dosages = dosages,
                                         distanceTH = range_threshold,
                                         rep1Colnames = s2_rep1_colnames,
                                         rep2Colnames = s2_rep2_colnames,
                                         rep1Omit = s2_rep1_omit,
                                         rep2Omit = s2_rep2_omit,
                                         bandwidthSet = bandwidth)
  
  
  
  double_fitPep_s2_output = outputTable_protein(protFit = double_fitPep_s2,
                                                uniprot_gn_filename = "/data/ginny/tcga_pancan/important_files/sp_0814_fasta_uniprot_gn.tsv",
                                                outputName = paste0(output_dir, s2_tag, "_protFitted.tsv"))
  
  
  
  cat("curves fitted.", "\n")
  
  ##### plot final output 
  
  
  
  
  arrPlot(phosphoSite1_fit = double_fitPhospho_s1,
          phosphoSite2_fit = double_fitPhospho_s2,
          pepSum1_fit = double_fitPep_s1,
          pepSum2_fit = double_fitPep_s2,
          sc1_tag = s1_tag,
          sc2_tag = s2_tag,
          orderS = order_s,
          sc1Rep1Omit = s1_rep1_omit,
          sc1Rep2Omit = s1_rep2_omit,
          sc2Rep1Omit = s2_rep1_omit,
          sc2Rep2Omit = s2_rep2_omit, #c(2,3)
          dosages = dosages,
          output_name = paste0(output_dir, "both_ProtPhospho_curve.pdf"))
  
  
  
  arrPlot_notBoth(phosphoSite1_fit = double_fitPhospho_s1,
                  phosphoSite2_fit = double_fitPhospho_s2,
                  pepSum1_fit = double_fitPep_s1,
                  pepSum2_fit = double_fitPep_s2,
                  sc1_tag = s1_tag,
                  sc2_tag = s2_tag,
                  bothSBL = both_s_BL,
                  orderS = order_s,
                  sc1Rep1Omit = s1_rep1_omit,
                  sc1Rep2Omit = s1_rep2_omit,
                  sc2Rep1Omit = s2_rep1_omit,
                  sc2Rep2Omit = s2_rep2_omit, #c(2,3)
                  dosages = dosages,
                  output_name = paste0(output_dir, "not_both_ProtPhospho_curve.pdf"))
  
  
  #### generate scatter plot with two temperatures 
  
  ### a function for scatterplot 
  
  ### for proteins 
  
  
  
  join_prot = joinScatterOutput_curve(
    protFit1 = double_fitPep_s1_output,
    protFit2 = double_fitPep_s2_output,
    slope_upper = 0.5,
    slope_lower = -0.5,
    outputName = paste0(output_dir, "prot_2s.tsv"),
    pdfName = paste0(output_dir, "prot_2s.pdf"))
  
  #### 0429 use ggplot to beautify the plot of prot_2s 
  
  
  cat("plotted.", "\n")
  cat("MISSION COMPLETED", "\n")
  
  
}

