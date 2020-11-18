


# base functions ----------------------------------------------------------


gauss.kernel.smooth = function(xx, yy, k.width, outlier=TRUE) {
  
  
  # xx = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  # yy = c(1.0000000,1.0706709,1.1524108,1.0434486,0.8909648,0.6184499,
  #        0.4114971,0.4031042,0.3749818,0.3804068)
   #k.width = 0.07
  # xx = sq1
  # xx = sq1/dosages
  # yy = norm_data1
  # k.width = 0.1
  # outlier = T

  yest = yy
  n = length(yy)
  medy = median(yy, na.rm=TRUE)
  mady = mad(yy, na.rm=TRUE)
  oid = abs(yy - medy) / mady > 9 ### relax the outlier threshold from 5 to 8 (times distance of mady)
  oid[is.na(oid)] = FALSE
  if(sum(oid)<=1)
  {
    yy[oid] = NA    #### catch outlier only if there is one outlier 
    
  }
  for(i in 1:n) {
    wt = (xx - xx[i]) / k.width
    wt = dnorm(wt, 0, 1)
    wt[is.na(yy)] = NA
    yest[i] = sum(wt * yy, na.rm=TRUE) / sum(wt, na.rm=TRUE)
  }
  
  
  yest
  
  
}
# 
#   plot(xx,yy, pch = 16, col = "red")
#   lines(xx, yest, type  = "b")
# ### function to calculate average slope 

cal_slope = function(xx, yy)
{
  seqL = length(yy)
  
  yy_lag = yy[-1]
  xx_lag = xx[-1]
  yy_diffs = yy_lag - yy[-seqL]
  xx_diffs = xx_lag - xx[-seqL]
  slopes = yy_diffs/xx_diffs
  
  av_slopes = mean(slopes)
  
  return(av_slopes)
  
}


#sign of the average slope of each interval
 
slope_sign = function(yy1, yy2, seq1, seq2)
{
  
  ### what is your expected results?
  
  # 
  # 
  # yy1 = this_curve1
  # yy2 = this_curve2
  # seq1 = sq1
  # seq2 = sq2

  sl1 = seq1[length(seq1)]
  sl2 = seq2[length(seq2)]
   
  sl = max(sl1, sl2)

  yy1_lag = yy1[-1]
  yy2_lag = yy2[-1]
  
  seq1_lag = seq1[-1]
  seq2_lag = seq2[-1]
  
 
  yy1_diffs = yy1_lag - yy1[-length(seq1)]
  yy2_diffs = yy2_lag - yy2[-length(seq2)]
  
  diff1 = rep(pi,(sl-1))
  diff1[seq1_lag-1] = yy1_diffs
  
  diff2 = rep(pi,(sl-1))
  diff2[seq2_lag-1] = yy2_diffs
  
  
  diffAv = diff1 + diff2
  
  which_miss = which(diffAv == 2*pi)
  
  if(length(which_miss)>0)
  {
    diff = diffAv[-which_miss]
    
  }else{
    diff = diffAv
  }

  

  each_slope = rep("-", length(diff))
  
  each_slope[which(diff>0)] = "+"
  
  
  return(each_slope)
  
  
  
}




slope_sign_single = function(yy1,  seq1)
{
   # yy1 = this_curve1
   # seq1 = sq1

  sl1 = seq1[length(seq1)]

  yy1_lag = yy1[-1]

  seq1_lag = seq1[-1]
  
  yy1_diffs = yy1_lag - yy1[-length(seq1)]

  diff1 = rep(pi,(sl1-1))
  diff1[seq1_lag-1] = yy1_diffs
  
  each_slope = rep("-", length(diff1))
  
  each_slope[which(diff1>0)] = "+"
  
  return(each_slope)
  
  
  
}


zero_length = function(x)
{
  which_zero = which(x ==0)
  return(length(which_zero))
}



# clean phospho site file -------------------------------------------------



### add the peptide sequence 


fileCleanPhosphoSites = function(original_file,
                                 rep1Colnames,
                                 rep2Colnames,
                                 probability_threshold,
                                 output_name)
{
  
  
  
  # 
  #  original_file =  s1_phospho_file
  # probability_threshold = prob_threshold
  # rep1Colnames = s1_rep1_colnames
  # rep2Colnames = s1_rep2_colnames
  # output_name = paste0(output_dir, s1_tag, "_phospho.tsv")
  # # 
  
  maxQuant = fread(original_file,
                   stringsAsFactors = F)
  
  
  old_colnames = colnames(maxQuant)
  
  new_colnames = unlist(lapply(1:length(old_colnames), function(p) {
    
    sep_names = unlist(strsplit(old_colnames[p], split = " ", fixed = T))
    new_name = paste(sep_names, collapse = "_")
    
    return(new_name)
    
    
  }))
  
  colnames(maxQuant) = new_colnames
  
  
  
  proteins_notIncluded = c("CON_")
  leadings_notIncluded = c("REV_")
  
  retain_data = maxQuant%>%
    dplyr::filter(!grepl(proteins_notIncluded, Proteins))%>%
    dplyr::filter(!grepl(leadings_notIncluded, Leading_proteins))
  
  ### next split the proteins 
  
  proteins = retain_data$Proteins
  
  
  #### change this part to include the header name B1/B2/B3
  
  
  data_col_rep1 = grep(paste0("Reporter_intensity_\\d+",rep1Colnames), colnames(retain_data))
  data_col_rep1_names =  grep(paste0("Reporter_intensity_\\d+",rep1Colnames), colnames(retain_data), value = T)
  
  
  data_col_rep2 = grep(paste0("Reporter_intensity_\\d+",rep2Colnames), colnames(retain_data))
  data_col_rep2_names =  grep(paste0("Reporter_intensity_\\d+",rep2Colnames), colnames(retain_data), value = T)
  
  data_col_names = c(data_col_rep1_names, data_col_rep2_names)
  data_col = c(data_col_rep1, data_col_rep2)
  
  unique_names = unique(unlist(lapply(1:length(data_col_names), function(q){
    this_name_sep = unlist(strsplit(data_col_names[q],split = "_" ))
    
    new_name = paste(this_name_sep[1:4], collapse = "_")
    return(new_name)
    
  })))
  
  
  detail_info = rbindlist(lapply(1:length(proteins), function(x) {


    disect_proteins = unlist(strsplit(proteins[x], split = ";"))
    
    this_positions =  as.integer(unlist(strsplit(retain_data$Positions_within_proteins[x], split = ";")))
    
    this_frame = data.frame(proteins = disect_proteins, 
                            leadProtein = retain_data$Protein[x],
                            positionsOnProteins = this_positions,
                            residue = retain_data$Amino_acid[x],
                            seq = retain_data$Sequence_window[x],
                            probability = retain_data$Localization_prob[x],
                            stringsAsFactors = F)
    
    
    this_data = retain_data[x, data_col]
    
    seq = c(1:20)
    first_seq= 3*(seq-1)+1
    first_pep = this_data[first_seq]
    
    
    second_pep = this_data[(first_seq+1)]
    third_pep = this_data[(first_seq+2)]
    
    names(first_pep) = unique_names
    names(second_pep) = unique_names
    names(third_pep) = unique_names
    
    long_data = rbind(rbind(first_pep, second_pep), third_pep)
    
    
    full_frame = rbindlist(lapply(1:nrow(this_frame), function(p){
      header = this_frame[p,]
      this_expand = data.frame(header[rep(1,3),], pepS = c(1,2,3), long_data)
      row.names(this_expand) = NULL
      return(this_expand)    
    }))
    
    
    
    if(x%%100 ==0)
      cat(x, "\n")
    
    return(full_frame)
    
  }))
  
  
  
  # 
  # detail_info = detail_info%>%
  #   dplyr::arrange(proteins, positionsOnProteins)
  # 
  # 
  # write.table(detail_info, "/data/ginny/YanTingPhospho/PhosphoSTYsitesDetail_long_test.tsv",
  #             sep = "\t", row.names = F, quote = F)
  
  detail_fil = detail_info%>%
    dplyr::arrange(proteins, positionsOnProteins)%>%
    dplyr::filter(probability>=probability_threshold)
  
  
  write.table(detail_fil, output_name,
              sep = "\t", row.names = F, quote = F)
  
  
  
  
}





# clean peptide file ------------------------------------------------------


fileCleanPeptides = function(original_file,
                             rep1Colnames,
                             rep2Colnames,
                             output_name)
{
  # 
  # original_file = "/data/ginny/YanTing20200227/InputFiles/B2B3_37_input/peptides.txt"
  # output_name = "/data/ginny/YanTing20200227/OutFiles/peptide37_Long.tsv"
  # 
  # 
  # 
  pep = fread(original_file,
              stringsAsFactors = F)
  
  
  
  old_colnames = colnames(pep)
  
  new_colnames = unlist(lapply(1:length(old_colnames), function(p) {
    
    sep_names = unlist(strsplit(old_colnames[p], split = " ", fixed = T))
    new_name = paste(sep_names, collapse = "_")
    
    return(new_name)
    
    
  }))
  
  colnames(pep) = new_colnames
  
  
  
  proteins_notIncluded = "CON_|REV_"
  
  retain_data = pep%>%
    dplyr::filter(!grepl(proteins_notIncluded, Proteins))%>%
    dplyr::filter(!grepl(proteins_notIncluded, Leading_razor_protein))
  
  ### next split the proteins 
  
  
  ### next split the proteins 
  
  proteins = retain_data$Proteins
  
  
  
  data_col_rep1 = grep(paste0("Reporter_intensity_\\d+",rep1Colnames), colnames(retain_data))
  data_col_rep1_names =  grep(paste0("Reporter_intensity_\\d+",rep1Colnames), colnames(retain_data), value = T)
  
  
  data_col_rep2 = grep(paste0("Reporter_intensity_\\d+",rep2Colnames), colnames(retain_data))
  data_col_rep2_names =  grep(paste0("Reporter_intensity_\\d+",rep2Colnames), colnames(retain_data), value = T)
  
  data_col_names = c(data_col_rep1_names, data_col_rep2_names)
  data_col = c(data_col_rep1, data_col_rep2)
  
  
  
  unique_names = unique(unlist(lapply(1:length(data_col_names), function(q){
    this_name_sep = unlist(strsplit(data_col_names[q],split = "_" ))
    
    new_name = paste(this_name_sep[1:4], collapse = "_")
    return(new_name)
    
  })))
  
  ####
  
  
  
  detail_info = rbindlist(lapply(1:length(proteins), function(x) {
    
    disect_proteins = unlist(strsplit(proteins[x], split = ";"))
    
    # this_positions =  as.integer(unlist(strsplit(pep$Positions_within_proteins[x], split = ";")))
    
    
    this_frame = data.frame(sequence = retain_data$Sequence[x],
                            leadProtein = retain_data$Leading_razor_protein[x],
                            startPosition = retain_data$Start_position[x],
                            endPosition = retain_data$End_position[x],
                            retain_data[x, data_col],
                            stringsAsFactors = F)
    
    full_frame = data.frame(proteins = disect_proteins, this_frame[rep(1, length(disect_proteins)),], stringsAsFactors = F)
    
    if(x%%100 ==0)
      cat(x, "\n")
    
    return(full_frame)
    
  }))
  
  
  
  
  detail_info = detail_info%>%
    dplyr::arrange(proteins, startPosition)
  
  write.table(detail_info, output_name,
              sep = "\t", row.names = F, quote = F)
  
  
}






# Filter and get double replicate phosphosites ----------------------------


filterDoublePhospho = function(Longfile_name,
                               rep_number = 10,
                               leadBL = T)
{
  # 
  # Longfile_name = "/data/ginny/YanTing20200227/OutFiles/phospho37_FilterLong.tsv"
  # rep_number = 10
  
  longf = fread(Longfile_name,
                stringsAsFactors = F)
  ### first remove rows with all zeros
  
  
  data_col = grep("Reporter_", colnames(longf), value = T)
  
  dataf = longf %>%
    dplyr::select(one_of(data_col))
  
  rowSum = apply(dataf, 1, sum)
  
  
  
  ### get rid of not necessary rows 
  
  cleanLongf = longf[which(rowSum>0),]
  
  
  clean_dataf = cleanLongf%>%
    dplyr::select(one_of(data_col))
  
  
  
  check_row_first = apply(clean_dataf[,1:rep_number], 1, zero_length)
  check_row_second = apply(clean_dataf[,(rep_number+1):(2*rep_number)], 1, zero_length)
  
  single_first_row = which(check_row_first == 0 & check_row_second == rep_number)
  single_second_row = which(check_row_first == rep_number  & check_row_second == 0)
  
  double_row = which(check_row_first == 0 & check_row_second == 0)
  
  single_firstf = cbind(cleanLongf[,1:6], clean_dataf[,1:rep_number])[single_first_row,]
  single_secondf = cbind(cleanLongf[,1:6], clean_dataf[,(rep_number+1):(2*rep_number)])[single_second_row,]
  doublef = cleanLongf[double_row,]
  
  
  
  if(leadBL == T)
  {
    
    doublef = doublef%>%
      dplyr::filter(proteins == leadProtein)
    
    
  }
  
  return(doublef)
  
  
}



#### fit curves for phosphosites with double replicates only save leading protein phosphosites 

# fit curves --------------------------------------------------------------



#which(doublePhospho_s1$proteins == "sp|Q9UK76-3|JUPI1_HUMAN" & doublePhospho_s1$pepS == 1 &doublePhospho_s1$positionsOnProteins == 85)

classify_curve = function(curve,
                          slope,
                          SlopeSeq,
                          disT)
{
  
  max_av = max(curve)
  min_av = min(curve)
  max_wav = which(curve == max_av)
  min_wav = which(curve == min_av)
  dis_av = max_av - min_av
  
  acl = length(curve)
  min_av_front = min(curve[c(1,2,3)])
  min_av_end = min(curve[c(acl-2,acl-1,acl)])
  
  
  this_response = "NR"
  b1 = F
  b1_2 = F
  b2 = F
  b3 = F
  b4 = F
  b5 = F
  b6 = F
  
  if(max_wav %in% c(4,5,6,7) & min_wav %in% c(1,2,3,8,9,10) )
    b1 = T
  
  if(min_wav %in% c(4,5,6,7) & max_wav %in% c(1,2,3,8,9,10) )
    b1_2 = T
  
  
  if(dis_av > disT)
    b2 = T
  
  
  if(grepl("---",SlopeSeq) & grepl("\\+\\+\\+", SlopeSeq))
    b3 = T
  
  ### change 6 to 5 consectives 
  
  if(grepl("-----",SlopeSeq))
    b4 = T
  
  if(grepl("\\+\\+\\+\\+\\+",SlopeSeq))
    b5 = T
  
  
  ### add a b6 to make sure the two ends does not have too much distance 
  
  if(abs(min_av_end - min_av_front)<0.45*dis_av)   ### pay attention to this 0.35 thing
    b6 = T
  
  
  if((b1== T | b1_2 == T) & b2 == T & b3 == T & b6 == T)
    this_response = "biphasic"
  # 
  
  if(b5==T & b2 == T & this_response!="biphasic"&slope>0)
    this_response = "hyper"
  
  if(b4==T & b2 == T& this_response!="biphasic"& slope<0)
    
    this_response = "hypo"
  
  
  #### a second round of adjustment 
  
  if(this_response == "NR" & slope >=0.3)
    this_response = "hyper"
  
  if(this_response == "NR" & slope <= -0.3)
    this_response = "hypo"
  
  
  if(this_response == "hyper" & slope <0.3)
    this_response = "NR"
  
  
  if(this_response == "hypo" & slope >-0.3)
    this_response = "NR"
  
  return(this_response)
  
}






curveFittingDouble = function(cleanedData,
                              dosages,
                              distanceTH,
                              rep1Colnames, 
                              rep2Colnames,
                              rep1Omit,
                              rep2Omit,
                              bandwidthSet)
  
  
  
{
  cat("starts curve fitting.", "\n")
  
  # cleanedData = doublePhospho_s1
  # dosages = dosages
  # distanceTH = range_threshold
  # rep1Colnames = s1_rep1_colnames
  # rep2Colnames = s1_rep2_colnames
  # rep1Omit = s1_rep1_omit
  # rep2Omit = s1_rep2_omit
  # bandwidthSet = bandwidth
  # # 
  
  
  data1_col = grep(rep1Colnames, colnames(cleanedData))
  data2_col = grep(rep2Colnames, colnames(cleanedData))
  info_col = setdiff(seq(1:ncol(cleanedData)), c(data1_col, data2_col))
  
  sq1 = seq(1:dosages)
  sq2 = seq(1:dosages)
  
  if(!is.na(rep1Omit))
  {
    data1_col = data1_col[-rep1Omit]
    
    sq1 = sq1[-rep1Omit]
  }
  
  if(!is.na(rep2Omit))
  {
    
    data2_col = data2_col[-rep2Omit]
    sq2 = sq2[-rep2Omit]
  }
  
  sl1 = length(sq1)
  sl2 = length(sq2)
  
  
  smooth_slope = rbindlist(lapply(1:nrow(cleanedData), function(i) {
    #smooth_slope = rbindlist(lapply(3245:3247, function(i) {
    
    
    #  i = 2171
    
    if(i%%1000 == 0)
      cat(i, "\n")
    this_data1 = as.numeric(cleanedData[i,data1_col])
    this_data2 = as.numeric(cleanedData[i,data2_col])
    
    ### normalise the data by dividing the first dosage 
    if(this_data1[1]!=0)
    {
      norm_data1 = this_data1/this_data1[1]
    }else
    {
      norm_data1 = rep(0,sl1)
    }
    
    if(this_data2[1]!=0)
    {
      norm_data2 = this_data2/this_data2[1]
    }else
    {
      norm_data2 = rep(0,sl2)
    }
    ## think about this 
    
    
    set.seed(123)
    this_curve1 = gauss.kernel.smooth(xx = sq1/dosages,
                                      yy = norm_data1,
                                      k.width = bandwidthSet, 
                                      outlier = T)
    
    
    curve1_avSlope = cal_slope(xx = sq1/dosages,
                               yy = this_curve1)
    
    
    
    
    set.seed(123)
    this_curve2 = gauss.kernel.smooth(xx = sq2/dosages,
                                      yy = norm_data2,
                                      k.width = bandwidthSet, 
                                      outlier = T)
    
    curve2_avSlope = cal_slope(xx = sq2/dosages,
                               yy = this_curve2)
    
    
    curve1_intervalSlope = slope_sign_single(yy1 = this_curve1,
                                             seq1 = sq1)
    
    slopeSeq1 = paste(curve1_intervalSlope, collapse = "")
    
    curve2_intervalSlope = slope_sign_single(yy1 = this_curve2,
                                             seq1 = sq2)
    
    slopeSeq2 = paste(curve2_intervalSlope, collapse = "")
    
    
    av_intervalSlope = slope_sign(yy1 = this_curve1,
                                  yy2 = this_curve2,
                                  seq1 = sq1,
                                  seq2 = sq2)
    
    avSlopeSeq = paste(av_intervalSlope, collapse = "")
    #### if the two replicates both have data 
    #### I will calculate a average version 
    
    final_slope = mean(c(curve1_avSlope, curve2_avSlope))
    
    cat_original = paste(c(norm_data1,norm_data2), collapse = "_")
    cat_smooth = paste(c(this_curve1,this_curve2), collapse = "_")
    
    # 
    # sl = max(sq1[length(sq1)],sq2[length(sq2)])
    # 
    # curve1 = rep(NA, sl)
    # curve2 = rep(NA, sl)
    # curve1[sq1] = this_curve1
    # curve2[sq2] = this_curve2
    # 
    # miss_curve1 = which(is.na(curve1))
    # miss_curve2 = which(is.na(curve2))
    # 
    # curve1[miss_curve1] = curve2[miss_curve1]
    # curve2[miss_curve2] = curve1[miss_curve2]
    # 
    # ### remove NA 
    # 
    # miss_curve1 = which(is.na(curve1))
    # miss_curve2 = which(is.na(curve2))
    # 
    # if(length(miss_curve1)>0)
    # {
    #   curve1 = curve1[-miss_curve1]
    #   
    # }
    # 
    # if(length(miss_curve2)>0)
    # {
    #   curve2 = curve2[-miss_curve2]
    #   
    # }
    # 
    # 
    # av_curve = curve1/2 + curve2/2
    # 
    
    #cat(i, this_curve1[1], slopeSeq1, "\n")
    response1 = classify_curve(curve = this_curve1, 
                               slope = curve1_avSlope,
                               disT = distanceTH,
                               SlopeSeq = slopeSeq1 )
    
    #cat(i, response1, "\n")
    
    #cat(i, this_curve2[1], slopeSeq2 ,"\n")
    
    response2 = classify_curve(curve = this_curve2, 
                               slope = curve2_avSlope,
                               disT = distanceTH,
                               SlopeSeq = slopeSeq2 )
    #cat(i, response2, "\n")
    
    #cat("##############################","\n")
    
    this_response = "NR"
    if(response1 == response2)
    {
      this_response = response1
    }
    
    
    sign_df = data.frame(t(av_intervalSlope), stringsAsFactors = F)
    colnames(sign_df) = paste("interval", seq(1:length(av_intervalSlope)), sep = "_")
    
    
    #### need to add other columns 
    ### kinase and so 
    #### in the final output 
    
    ## add in a separate function tomorrow 
    
    
    this_df = data.frame(cleanedData[i, info_col], 
                         original = cat_original,
                         smooth = cat_smooth,
                         firstSlope = curve1_avSlope,
                         secondSlope = curve2_avSlope,
                         firstClass = response1,
                         secondClass = response2,
                         response = this_response,
                         avSlope = final_slope, 
                         stringsAsFactors = F)
    
    
    finish_df = cbind(this_df, sign_df)
    
    
    #}
    
    
    return(finish_df)
    
    
    
  }))
  
  
  arr_smooth_slope = smooth_slope%>%
    dplyr::arrange(avSlope)
  
  
  
  # # drop_cols = c("original", "smooth")
  # # 
  # # output_dataf = arr_smooth_slope%>%
  # #   dplyr::select(-one_of(drop_cols))
  # # 
  # # write.table(output_dataf, outputName,
  # #             quote = F, sep = "\t", row.names = F)
  # 
  return (arr_smooth_slope)
  #return(smooth_slope)
  
}








curveFittingDouble_old = function(cleanedData,
                              dosages,
                              distanceTH,
                              rep1Colnames, 
                              rep2Colnames,
                              rep1Omit,
                              rep2Omit,
                              bandwidthSet)
  
  
  
 {
  # cleanedData = doublePhospho_s1
  # dosages = dosages
  # distanceTH = range_threshold
  # rep1Colnames = s1_rep1_colnames
  # rep2Colnames = s1_rep2_colnames
  # rep1Omit = s1_rep1_omit
  # rep2Omit = s1_rep2_omit
  # bandwidthSet = bandwidth
  # 
  
  data1_col = grep(rep1Colnames, colnames(cleanedData))
  data2_col = grep(rep2Colnames, colnames(cleanedData))
  info_col = setdiff(seq(1:ncol(cleanedData)), c(data1_col, data2_col))
  
  sq1 = seq(1:dosages)
  sq2 = seq(1:dosages)
  
  if(!is.na(rep1Omit))
  {
    data1_col = data1_col[-rep1Omit]
    
    sq1 = sq1[-rep1Omit]
  }
  
  if(!is.na(rep2Omit))
  {
    
    data2_col = data2_col[-rep2Omit]
    sq2 = sq2[-rep2Omit]
  }
  
  sl1 = length(sq1)
  sl2 = length(sq2)
  
  
  smooth_slope = rbindlist(lapply(1:nrow(cleanedData), function(i) {

#   i = 7766


 if(i%%1000 == 0)
      cat(i, "\n")
    this_data1 = as.numeric(cleanedData[i,data1_col])
    this_data2 = as.numeric(cleanedData[i,data2_col])
    
    ### normalise the data by dividing the first dosage 
    if(this_data1[1]!=0)
    {
      norm_data1 = this_data1/this_data1[1]
    }else
    {
      norm_data1 = rep(0,sl1)
    }
    
    if(this_data2[1]!=0)
    {
      norm_data2 = this_data2/this_data2[1]
    }else
    {
      norm_data2 = rep(0,sl2)
    }
    # 
    #   norm_data1 = this_data1
    #   norm_data2 = this_data2 
    #   
    #   
    
    ### think about this 
    
    
    set.seed(123)
    this_curve1 = gauss.kernel.smooth(xx = sq1/dosages,
                                      yy = norm_data1,
                                      k.width = bandwidthSet, 
                                      outlier = T)
    
    
    curve1_avSlope = cal_slope(xx = sq1/dosages,
                               yy = this_curve1)
    
    
    
    
    set.seed(123)
    this_curve2 = gauss.kernel.smooth(xx = sq2/dosages,
                                      yy = norm_data2,
                                      k.width = bandwidthSet, 
                                      outlier = T)
    
    curve2_avSlope = cal_slope(xx = sq2/dosages,
                               yy = this_curve2)
    
    
    
    av_intervalSlope = slope_sign(yy1 = this_curve1,
                                  yy2 = this_curve2,
                                  seq1 = sq1,
                                  seq2 = sq2)
    
    avSlopeSeq = paste(av_intervalSlope, collapse = "")
    #### if the two replicates both have data 
    #### I will calculate a average version 
    
    final_slope = mean(c(curve1_avSlope, curve2_avSlope))
    
    cat_original = paste(c(norm_data1,norm_data2), collapse = "_")
    cat_smooth = paste(c(this_curve1,this_curve2), collapse = "_")
    
    
    
    
    #### ok this part design a way to evaluate whether it is biphasic or not 
    
    ### a new function to do it 
    
    #### from the fitted curve, calculate the distance between max and min 
    
    #### think about diff 1/2 length 
    ### not quite right 
    

    sl = max(sq1[length(sq1)],sq2[length(sq2)])

    curve1 = rep(NA, sl)
    curve2 = rep(NA, sl)
    curve1[sq1] = this_curve1
    curve2[sq2] = this_curve2
    
    miss_curve1 = which(is.na(curve1))
    miss_curve2 = which(is.na(curve2))
    
    curve1[miss_curve1] = curve2[miss_curve1]
    curve2[miss_curve2] = curve1[miss_curve2]
    
    ### remove NA 
    
    miss_curve1 = which(is.na(curve1))
    miss_curve2 = which(is.na(curve2))
    
    if(length(miss_curve1)>0)
    {
      curve1 = curve1[-miss_curve1]
      
    }
    
    if(length(miss_curve2)>0)
    {
      curve2 = curve2[-miss_curve2]
      
    }
    

    av_curve = curve1/2 + curve2/2
    
    
    #### the above remain unchanged 
    
    max_av = max(av_curve)
    min_av = min(av_curve)
    max_wav = which(av_curve == max_av)
    min_wav = which(av_curve == min_av)
    dis_av = max_av - min_av
    
    acl = length(av_curve)
    min_av_front = min(av_curve[c(1,2,3)])
    min_av_end = min(av_curve[c(acl-2,acl-1,acl)])
    
    
    # 
    # sp = paste(rep("+", (length(common_seq)-1)), collapse = "")
    # sn = paste(rep("-", (length(common_seq)-1)), collapse = "")
    # 
    # 
    # 
    
    this_response = "NR"
    b1 = F
    b2 = F
    b3 = F
    b4 = F
    b5 = F
    b6 = F
    
    if(max_wav %in% c(4,5,6,7) & min_wav %in% c(1,2,3,8,9,10) )
      b1 = T
    if(dis_av > distanceTH)
      b2 = T
    
    
    if(grepl("---",avSlopeSeq) & grepl("\\+\\+\\+", avSlopeSeq))
      b3 = T
    
    ### change 6 to 5 consectives 
    
    if(grepl("-----",avSlopeSeq))
      b4 = T
    
    if(grepl("\\+\\+\\+\\+\\+",avSlopeSeq))
      b5 = T
    
    
    
    ### add a b6 to make sure the two ends does not have too much distance 
    
    if(abs(min_av_end - min_av_front)<0.3*dis_av)   ### pay attention to this 0.35 thing
      b6 = T
    
    
    
    
    if(b1== T & b2 == T & b3 == T & b6 == T)
      this_response = "biphasic"

    
    # if(b5==T & b2 == T & this_response!="biphasic" & final_slope>0)
    #   this_response = "hyper"
    # 
    # if(b4==T & b2 == T& this_response!="biphasic"& final_slope<0)
    #   this_response = "hypo"
    # 
    
    if(b5==T & b2 == T & this_response!="biphasic"&final_slope>0)
      this_response = "hyper"
    
    if(b4==T & b2 == T& this_response!="biphasic"& final_slope<0)
      
      this_response = "hypo"
    
    
    #### a second round of adjustment 
    
    if(this_response == "NR" & final_slope >=0.3)
      this_response = "hyper"
    
    if(this_response == "NR" & final_slope <= -0.3)
      this_response = "hypo"
    
    
    if(this_response == "hyper" & final_slope <0.3)
      this_response = "NR"
    
    
    if(this_response == "hypo" & final_slope >-0.3)
      this_response = "NR"
    
    
    
    #### for hyper and hypo ask for slope!!!!
    
    
    sign_df = data.frame(t(av_intervalSlope), stringsAsFactors = F)
    colnames(sign_df) = paste("interval", seq(1:length(av_intervalSlope)), sep = "_")
    
    
    #### need to add other columns 
    ### kinase and so 
    #### in the final output 
    
    ## add in a separate function tomorrow 
    
    
    this_df = data.frame(cleanedData[i, info_col], 
                         original = cat_original,
                         smooth = cat_smooth,
                         firstSlope = curve1_avSlope,
                         secondSlope = curve2_avSlope,
                         response = this_response,
                         avSlope = final_slope, 
                         stringsAsFactors = F)
    
    
    finish_df = cbind(this_df, sign_df)
    
    
   #}


    return(finish_df)



  }))


  arr_smooth_slope = smooth_slope%>%
    dplyr::arrange(avSlope)
  
  # drop_cols = c("original", "smooth")
  # 
  # output_dataf = arr_smooth_slope%>%
  #   dplyr::select(-one_of(drop_cols))
  # 
  # write.table(output_dataf, outputName,
  #             quote = F, sep = "\t", row.names = F)
  
  return (arr_smooth_slope)
  
  
}







# output a table  --------------------------------------


## this is for phoshposite 

outputTable = function(phosphoSiteFit,
                       uniprot_gn_filename,
                       kinaseSubstrate_filename,
                       outputName)
{

  # # 
  # 
  # phosphoSiteFit = double_fitPhospho_s1
  # uniprot_gn_filename = "/data/ginny/tcga_pancan/important_files/sp_0814_fasta_uniprot_gn.tsv"
  # kinaseSubstrate_filename = ks_network_file
  # outputName = paste0(output_dir, s1_tag, "_phosphoFitted.tsv")

  #### what are the info needed here 
  drop_cols = c("original", "smooth")
  get_phospho = phosphoSiteFit%>%
    dplyr::select(-one_of(drop_cols))
  
  
  
  prot = get_phospho$proteins 
  
  uni_prot = rbindlist(lapply(1:length(prot), function(x) {
    
    this_acc = unlist(strsplit(prot[x], split = "|", fixed = T))[2]
    this_ph = unlist(strsplit(prot[x], split = "|", fixed = T))[3]
    this_uni = unlist(strsplit(this_acc, split = "-", fixed = T))[1]
     
    this_pn = unlist(strsplit(this_ph,split = "_", fixed = T))[1]
    
    this_df = data.frame(uniprot = this_uni, protName = this_pn, stringsAsFactors = F)
    
    return(this_df)
    
  }))
  


  phoFit = cbind(uni_prot, get_phospho)
  
  
  
  uniprot_gn = fread(uniprot_gn_filename, stringsAsFactors = F)
  
    ks = fread(kinaseSubstrate_filename, stringsAsFactors = F)
  
  
  drop_cols = c("source1","source2","source3")
  
  
  
  
  phoFit_gn = phoFit%>%
    dplyr::left_join(uniprot_gn, by = c("uniprot" = "accession"))%>%
    dplyr::mutate(posRes = paste0(residue, positionsOnProteins))%>%
    dplyr::left_join(ks, by = c("geneName" = "substrate", "posRes" = "substrate_site"))%>%
    dplyr::select(-one_of(drop_cols))%>%
    dplyr::group_by_at(vars(uniprot:posRes))%>%
    dplyr::summarise(allK = paste(kinase,collapse = "_"))%>%
    dplyr::ungroup()%>%
    dplyr::arrange(avSlope)%>%
    dplyr::select(uniprot, protName, geneName, posRes, allK, everything())
  
  
  ### I think ready to be output now 
  write.table(phoFit_gn, outputName, sep = "\t", quote = F, row.names = F)
  
  return(phoFit_gn)
  
}
  

# output a table for protein results --------------------------------------



outputTable_protein = function(protFit,
                       uniprot_gn_filename,
                       outputName)
{
  
  # 
  #  protFit = double_fitPep_s1
  # uniprot_gn_filename = "/data/ginny/tcga_pancan/important_files/sp_0814_fasta_uniprot_gn.tsv"
  # outputName = paste0(output_dir, s1_tag, "_protFitted.tsv")
  # # 
  #### what are the info needed here 
  drop_cols = c("original", "smooth")
  get_prot = protFit%>%
    dplyr::select(-one_of(drop_cols))
  
  
  
  prot = get_prot$proteins 
  
  uni_prot = rbindlist(lapply(1:length(prot), function(x) {
  

    this_acc = unlist(strsplit(prot[x], split = "|", fixed = T))[2]
    this_ph = unlist(strsplit(prot[x], split = "|", fixed = T))[3]
    this_uni = unlist(strsplit(this_acc, split = "-", fixed = T))[1]
    this_iso = unlist(strsplit(this_acc, split = "-", fixed = T))[2]
    
    this_pn = unlist(strsplit(this_ph,split = "_", fixed = T))[1]
    
    

      

    
    this_df = data.frame(uniprot = this_uni, protName = this_pn,
                         isoform = this_iso,
                         stringsAsFactors = F)
    
    return(this_df)
    
  }))
  
  
  
  pFit = cbind(uni_prot, get_prot)
  
  
  
  
  
  # 
   uniprot_gn = fread(uniprot_gn_filename, stringsAsFactors = F)
  # 

   
   #### 0511 one more extra step, get the isoform from  the protein and make it another symbol 
   

  
  
  # 
  pFit_ar = pFit%>%
    dplyr::left_join(uniprot_gn, by = c("uniprot" = "accession"))%>%
    dplyr::arrange(avSlope)%>%
    dplyr::mutate(symbol = case_when(is.na(isoform) ~ geneName,
                                     T ~ paste0(geneName,"-", isoform)))%>%
    dplyr::select(uniprot, protName, geneName,isoform, symbol, everything())
    
  # 
  
  ### I think ready to be output now 
  write.table(pFit_ar, outputName, sep = "\t", quote = F, row.names = F)
  
  return(pFit_ar)
  
}






# get peptide level doubles only ------------------------------------------


filterDoublePeptide = function(pep_filename,
                               rep_number = 10,
                               leadBL = F)
  
{
  
  
  
  
  
  pepf = fread(pep_filename,
               stringsAsFactors = F)
  
  
  data_col = grep("Reporter_", colnames(pepf), value = T)
  
  dataf = pepf %>%
    dplyr::select(one_of(data_col))
  
  rowSum = apply(dataf, 1, sum)
  
  
  
  ### get rid of not necessary rows 
  
  cleanPepf = pepf[which(rowSum>0),]
  
  
  
  
  
  #### from cleanLongf I want to separate double replicate and single replicate 
  
  
  clean_dataf = cleanPepf%>%
    dplyr::select(one_of(data_col))
  
  
  
  check_row_first = apply(clean_dataf[,1:rep_number], 1, zero_length)
  check_row_second = apply(clean_dataf[,(rep_number+1):(2*rep_number)], 1, zero_length)
  
  single_first_row = which(check_row_first == 0 & check_row_second == rep_number)
  single_second_row = which(check_row_first == rep_number  & check_row_second == 0)
  
  double_row = which(check_row_first == 0 & check_row_second == 0)
  
  single_firstf = cbind(cleanPepf[,1:5], clean_dataf[,1:rep_number])[single_first_row,]
  single_secondf = cbind(cleanPepf[,1:5], clean_dataf[,(rep_number+1):(2*rep_number)])[single_second_row,]
  doublef = cleanPepf[double_row,]
  
  if(leadBL == T)
  {
    
    doublef = doublef%>%
      dplyr::filter(proteins == leadProtein)
    
    
  }
  
  return(doublef)
  
  
  
}




# sum peptide to protein --------------------------------------------------


sumPeptideToProtein = function(peptide,
                               leadBL = F,
                               rep1Colnames,
                               rep2Colnames,
                               output_name)
{
  #peptide_filename = "/data/ginny/YanTing20200227/OutFiles/peptide37.tsv"
  
  # 
  #   peptide = doublePep_s1
  #   leadBL = lead_peptide_BL
  #   output_name = paste0(output_dir, s1_tag,"_sum_peptide.tsv")
  #   
  
  
  if(leadBL == T)
  {
    peptide = peptide %>%
      dplyr::filter(proteins == leadProtein)
  }
  
  
  #### I will have a bug if user is not using B2 AND B3
  first_var = quo(paste0("Reporter_intensity_1", rep1Colnames))
  last_var = quo(paste0("Reporter_intensity_10", rep2Colnames))
  
  
  
  
  agg_pep = peptide%>%
    dplyr::group_by(proteins)%>%
    dplyr::summarise_at(vars(!!first_var:!!last_var), sum)%>%
    dplyr::ungroup()
  
  
  write.table(agg_pep, output_name, quote = F, row.names =  F, sep = "\t")
  return(agg_pep)
  
}






# arrange and plot data ---------------------------------------------------


arrPlot = function(phosphoSite1_fit,
                   phosphoSite2_fit,
                   pepSum1_fit,
                   pepSum2_fit,
                   sc1_tag,
                   sc2_tag,
                   #  bothSBL = T,
                   orderS = 1,
                   sc1Rep1Omit = NA,
                   sc1Rep2Omit = NA,
                   sc2Rep1Omit = 2,
                   sc2Rep2Omit = 2,
                   dosages = 10,
                   output_name)



{
  # # 
  # phosphoSite1_fit = double37Lead_fit03
  # phosphoSite2_fit = double52Lead_fit03
  # pepSum1_fit = double37Pep_fit03
  # pepSum2_fit = double52Pep_fit03
  # sc1_tag = "37C"
  # sc2_tag = "52C"
  # bothSBL = T
  # orderS = 1
  # sc1Rep1Omit = NA
  # sc1Rep2Omit = NA
  # sc2Rep1Omit = 2
  # sc2Rep2Omit = 2  #c(2,3)
  # dosages = 10
  # output_name = "/data/ginny/YanTing20200227/OutFiles/double_arrPlots_04.pdf"
  # # 
  
  ### simply plot 
  
  
  ## proteins with both temperatures 
  
  orderProt = pepSum1_fit
  followProt = pepSum2_fit
  
  orderPhospho = phosphoSite1_fit
  followPhospho = phosphoSite2_fit
  
  order_tag = sc1_tag
  follow_tag = sc2_tag
  orderRep1Omit = sc1Rep1Omit
  orderRep2Omit = sc1Rep2Omit
  followRep1Omit = sc2Rep1Omit
  followRep2Omit = sc2Rep2Omit
  
  
  if(orderS == 2)
  {
    
    orderProt = pepSum2_fit
    followProt = pepSum1_fit
    
    orderPhospho = phosphoSite2_fit
    followPhospho = phosphoSite1_fit
    
    order_tag = sc2_tag
    follow_tag = sc1_tag
    
    orderRep1Omit = sc2Rep1Omit
    orderRep2Omit = sc2Rep2Omit
    followRep1Omit = sc1Rep1Omit
    followRep2Omit = sc1Rep2Omit
    
    
    
  }
  
  
  
  
  # 
  # prot_order = intersect(orderPhospho$proteins, orderProt$proteins)
  # prot_follow = intersect(followPhospho$proteins, followProt$proteins)
  
  prot_order =  orderProt$proteins
  prot_follow = followProt$proteins
  
  
  
  protInt = intersect(prot_order, prot_follow)
  
  #  if(bothSBL == T)
  # {
  prot_order = protInt
  prot_follow = protInt
  
  #}
  
  
  
  #### order according to prot1
  
  
  orderProt_order = orderProt%>%
    dplyr::filter(proteins %in% prot_order)%>%
    dplyr::arrange(avSlope)
  
  
  #### for each of these protein find the peptides and plot them 
  pdf(output_name, useDingbats = F,width = 12, height = 12)
  
  par(mfrow = c(4,4))
  
  
  for(i in 1:nrow(orderProt_order))
    
    
  {
    
    
    this_protein = orderProt_order$proteins[i]
    
    #### order according to peptides in the first temperature 
    
    
    this_orderPhospho = orderPhospho%>%
      dplyr::filter(proteins == this_protein)
    
    
    this_followPhospho = followPhospho%>%
      dplyr::filter(proteins == this_protein)
    
    this_orderProt = orderProt%>%
      dplyr::filter(proteins == this_protein)
    
    this_followProt = followProt%>%
      dplyr::filter(proteins == this_protein)
    
    
    ### now let me only consider the scenario where only intersection proteins are considered
    
    
    orderProt_original = as.numeric(unlist(strsplit(this_orderProt$original[1], split = "_")))
    orderProt_smooth = as.numeric(unlist(strsplit(this_orderProt$smooth[1], split = "_")))
    
    ### separate the two replicate's data 
    rep_seq = seq(1:dosages)
    
    order_rep1 = rep_seq
    order_rep2 = rep_seq
    
    if(is.na(orderRep1Omit)==F)
    {
      order_rep1 = rep_seq[-orderRep1Omit]
      
    }
    
    if(is.na(orderRep2Omit) == F)
    {
      order_rep2 = rep_seq[-orderRep2Omit]+dosages
      
    }
    
    or1 = seq(1:length(order_rep1))
    or2 = c((1+length(order_rep1)):(length(order_rep1)+length(order_rep2)))
    
    orderProt_data1 = orderProt_original[or1]
    orderProt_curve1 = orderProt_smooth[or1]
    
    orderProt_data2 = orderProt_original[or2]
    orderProt_curve2 = orderProt_smooth[or2]
    
    
    
    ############### 
    followProt_original = as.numeric(unlist(strsplit(this_followProt$original[1], split = "_")))
    followProt_smooth = as.numeric(unlist(strsplit(this_followProt$smooth[1], split = "_")))
    
    follow_rep1 = rep_seq
    follow_rep2 = rep_seq
    
    if(is.na(followRep1Omit)==F)
    {
      follow_rep1 = rep_seq[-followRep1Omit]
      
    }
    
    if(is.na(followRep2Omit) == F)
    {
      follow_rep2 = rep_seq[-followRep2Omit]
      
    }
    fr1 = seq(1:length(follow_rep1))
    fr2 = c((1+length(follow_rep1)):(length(follow_rep1)+length(follow_rep2)))
    
    followProt_data1 = followProt_original[fr1]
    followProt_curve1 = followProt_smooth[fr1]
    
    followProt_data2 = followProt_original[fr2]
    followProt_curve2 = followProt_smooth[fr2]
    
    
    #######
    
    orderProt_xlab = paste0("Doses    Blue Rep1_", round(this_orderProt$firstSlope[1], digits = 2),
                            "    Red Rep2_",
                            round(this_orderProt$secondSlope[1], digits = 2))
    
    followProt_xlab = paste0("Doses    Blue Rep1_", round(this_followProt$firstSlope[1], digits = 2),
                             "    Red Rep2_",
                             round(this_followProt$secondSlope[1], digits = 2))
    
    
    orderProt_max = max(c(orderProt_original, orderProt_smooth))
    followProt_max = max(c(followProt_original, followProt_smooth))
    
    
    
    #### same row for the same phosphoSite 
    #### think about no rows scenario 
    
    ### order  this_orderPhospho
    
    # order_this_orderPhospho = this_orderPhospho%>%
    #   dplyr::arrange(avSlope)
    # ### think about how to identify phosphosites 
    
    ### take care of the ones appear in both temperatures
    join_order_follow = this_orderPhospho%>%
      dplyr::inner_join(this_followPhospho, by = c("proteins",
                                                   "leadProtein",
                                                   "positionsOnProteins",
                                                   "residue"))%>%
      dplyr::arrange(avSlope.x)
    
    
    if(nrow(join_order_follow)>0)
    {
      for(k in 1:nrow(join_order_follow))
      {
        
        ### make the name 
        
        cat(i, k , "\n")
        this_name = paste(join_order_follow$proteins[k], 
                          join_order_follow$positionsOnProteins[k],
                          join_order_follow$residue[k],
                          sep = "-")
        
        order_name = paste0(this_name, "-S", join_order_follow$pepS.x[k])
        follow_name = paste0(this_name, "-S", join_order_follow$pepS.y[k])
        
        
        
        ##### plot them
        
        order_original = as.numeric(unlist(strsplit(join_order_follow$original.x[k], split = "_")))
        order_smooth = as.numeric(unlist(strsplit(join_order_follow$smooth.x[k], split = "_")))
        
        
        
        follow_original = as.numeric(unlist(strsplit(join_order_follow$original.y[k], split = "_")))
        follow_smooth = as.numeric(unlist(strsplit(join_order_follow$smooth.y[k], split = "_")))
        
        
        order_data1 = order_original[or1]
        order_curve1 = order_smooth[or1]
        
        order_data2 = order_original[or2]
        order_curve2 = order_smooth[or2]
        
        
        
        follow_data1 = follow_original[fr1]
        follow_curve1 = follow_smooth[fr1]
        
        follow_data2 = follow_original[fr2]
        follow_curve2 = follow_smooth[fr2]
        
        # so1 = length(order_data1)
        # so2 = length(order_data2)
        # sf1 = length(follow_data1)
        # sf2 = length(follow_data2)
        # 
        
        
        order_min = min(c(order_original, order_smooth))
        order_max = max(c(order_original, order_smooth))
        
        follow_min = min(c(follow_original, follow_smooth))
        follow_max = max(c(follow_original, follow_smooth))
        
        
        order_xlab = paste0("Doses    Blue Rep1_", round(join_order_follow$firstSlope.x[k], digits = 2),
                            "    Red Rep2_",
                            round(join_order_follow$secondSlope.x[k], digits = 2))
        
        follow_xlab = paste0("Doses    Blue Rep1_", round(join_order_follow$firstSlope.y[k], digits = 2),
                             "    Red Rep2_",
                             round(join_order_follow$secondSlope.y[k], digits = 2))
        
        
        ## for the protein 
        plot(NULL, xlab = orderProt_xlab, ylab = paste0("prot abundance ", order_tag), 
             xlim = c(1,dosages), ylim = c(0, orderProt_max),
             main = paste0(this_protein,"_",round(this_orderProt$avSlope[1], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(order_rep1, orderProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(order_rep1, orderProt_curve1, type = "l", col = "blue")
        lines(order_rep2, orderProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(order_rep2, orderProt_curve2, type = "l", col = "red", pch = 16)
        
        
        
        plot(NULL, xlab = order_xlab, ylab = paste0("ptm abundance ", order_tag), 
             xlim = c(1,dosages), ylim = c(0, order_max),
             main = paste0(order_name,"_",round(join_order_follow$avSlope.x[k], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(order_rep1, order_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(order_rep1, order_curve1, type = "l", col = "blue")
        lines(order_rep2, order_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(order_rep2, order_curve2, type = "l", col = "red")
        
        
        
        ## for the protein 
        
        plot(NULL, xlab = followProt_xlab, ylab = paste0("prot abundance ", follow_tag), 
             xlim = c(1,dosages), ylim = c(0, followProt_max),
             main = paste0(this_protein,"_",round(this_followProt$avSlope[1], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(follow_rep1, followProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(follow_rep1, followProt_curve1, type = "l", col = "blue")
        lines(follow_rep2, followProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(follow_rep2, followProt_curve2, type = "l", col = "red", pch = 16)
        
        
        
        
        ### for the follow site 
        
        plot(NULL, xlab = follow_xlab, ylab = paste0("ptm abundance ", follow_tag), 
             xlim = c(1,dosages), ylim = c(0, follow_max),
             main = paste0(follow_name,"_",round(join_order_follow$avSlope.y[k], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(follow_rep1, follow_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(follow_rep1, follow_curve1, type = "l", col = "blue")
        lines(follow_rep2, follow_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(follow_rep2, follow_curve2, type = "l", col = "red")
        
        
        
        
      }
      
      
    }
    
    
    
    #dev.off()
    
  }
  
  
  dev.off()
  
  
  
  
}




# arrange and plot but with either or none temperatures -------------------


arrPlot_notBoth = function(phosphoSite1_fit,
                           phosphoSite2_fit,
                           pepSum1_fit,
                           pepSum2_fit,
                           sc1_tag,
                           sc2_tag,
                           bothSBL = T,
                           orderS = 1,
                           sc1Rep1Omit = NA,
                           sc1Rep2Omit = NA,
                           sc2Rep1Omit = 2,
                           sc2Rep2Omit = 2,
                           dosages = 10,
                           output_name)



{
  #   # 
  # phosphoSite1_fit = double37Lead_fit03
  # phosphoSite2_fit = double52Lead_fit03
  # pepSum1_fit = double37Pep_fit03
  # pepSum2_fit = double52Pep_fit03
  # sc1_tag = "37C"
  # sc2_tag = "52C"
  # bothSBL = T
  # orderS = 1
  # sc1Rep1Omit = NA
  # sc1Rep2Omit = NA
  # sc2Rep1Omit = 2
  # sc2Rep2Omit = 2  #c(2,3)
  # dosages = 10
  # output_name = "/data/ginny/YanTing20200227/OutFiles/double_notBoth_arrPlots_04.pdf"
  #   # 
  #   
  ### simply plot 
  
  
  ## proteins with both temperatures 
  
  orderProt = pepSum1_fit
  followProt = pepSum2_fit
  
  orderPhospho = phosphoSite1_fit
  followPhospho = phosphoSite2_fit
  
  order_tag = sc1_tag
  follow_tag = sc2_tag
  orderRep1Omit = sc1Rep1Omit
  orderRep2Omit = sc1Rep2Omit
  followRep1Omit = sc2Rep1Omit
  followRep2Omit = sc2Rep2Omit
  
  
  if(orderS == 2)
  {
    
    orderProt = pepSum2_fit
    followProt = pepSum1_fit
    
    orderPhospho = phosphoSite2_fit
    followPhospho = phosphoSite1_fit
    
    order_tag = sc2_tag
    follow_tag = sc1_tag
    
    orderRep1Omit = sc2Rep1Omit
    orderRep2Omit = sc2Rep2Omit
    followRep1Omit = sc1Rep1Omit
    followRep2Omit = sc1Rep2Omit
    
    
    
  }
  
  
  
  
  # 
  # prot_order = intersect(orderPhospho$proteins, orderProt$proteins)
  # prot_follow = intersect(followPhospho$proteins, followProt$proteins)
  
  
  prot_order = orderProt$proteins
  prot_follow = followProt$proteins
  
  
  
  
  
  protInt = intersect(prot_order, prot_follow)
  
  if(bothSBL == T)
  {
    prot_order = protInt
    prot_follow = protInt
    
  }
  
  
  
  #### order according to prot1
  
  
  orderProt_order = orderProt%>%
    dplyr::filter(proteins %in% prot_order)%>%
    dplyr::arrange(avSlope)
  
  
  #### for each of these protein find the peptides and plot them 
  pdf(output_name, useDingbats = F,width = 12, height = 12)
  
  par(mfrow = c(4,4))
  
  
  for(i in 1:nrow(orderProt_order))
    
    
  {
    # i =40
    
    this_protein = orderProt_order$proteins[i]
    
    #### order according to peptides in the first temperature 
    
    
    this_orderPhospho = orderPhospho%>%
      dplyr::filter(proteins == this_protein)
    
    
    this_followPhospho = followPhospho%>%
      dplyr::filter(proteins == this_protein)
    
    this_orderProt = orderProt%>%
      dplyr::filter(proteins == this_protein)
    
    
    
    this_followProt = followProt%>%
      dplyr::filter(proteins == this_protein)
    
    
    ### now let me only consider the scenario where only intersection proteins are considered
    
    
    orderProt_original = as.numeric(unlist(strsplit(this_orderProt$original[1], split = "_")))
    orderProt_smooth = as.numeric(unlist(strsplit(this_orderProt$smooth[1], split = "_")))
    
    ### separate the two replicate's data 
    rep_seq = seq(1:dosages)
    
    order_rep1 = rep_seq
    order_rep2 = rep_seq
    
    if(is.na(orderRep1Omit)==F)
    {
      order_rep1 = rep_seq[-orderRep1Omit]
      
    }
    
    if(is.na(orderRep2Omit) == F)
    {
      order_rep2 = rep_seq[-orderRep2Omit]+dosages
      
    }
    
    or1 = seq(1:length(order_rep1))
    or2 = c((1+length(order_rep1)):(length(order_rep1)+length(order_rep2)))
    
    orderProt_data1 = orderProt_original[or1]
    orderProt_curve1 = orderProt_smooth[or1]
    
    orderProt_data2 = orderProt_original[or2]
    orderProt_curve2 = orderProt_smooth[or2]
    
    
    orderProt_xlab = paste0("Doses    Blue Rep1_", round(this_orderProt$firstSlope[1], digits = 2),
                            "    Red Rep2_",
                            round(this_orderProt$secondSlope[1], digits = 2))
    
    orderProt_max = max(c(orderProt_original, orderProt_smooth))
    
    
    
    ############### 
    
    followProt_flag = T
    
    if(nrow(this_followProt)>0)
      
    {
      followProt_original = as.numeric(unlist(strsplit(this_followProt$original[1], split = "_")))
      followProt_smooth = as.numeric(unlist(strsplit(this_followProt$smooth[1], split = "_")))
      
      follow_rep1 = rep_seq
      follow_rep2 = rep_seq
      
      if(is.na(followRep1Omit)==F)
      {
        follow_rep1 = rep_seq[-followRep1Omit]
        
      }
      
      if(is.na(followRep2Omit) == F)
      {
        follow_rep2 = rep_seq[-followRep2Omit]
        
      }
      fr1 = seq(1:length(follow_rep1))
      fr2 = c((1+length(follow_rep1)):(length(follow_rep1)+length(follow_rep2)))
      
      followProt_data1 = followProt_original[fr1]
      followProt_curve1 = followProt_smooth[fr1]
      
      followProt_data2 = followProt_original[fr2]
      followProt_curve2 = followProt_smooth[fr2]
      
      
      followProt_xlab = paste0("Doses    Blue Rep1_", round(this_followProt$firstSlope[1], digits = 2),
                               "    Red Rep2_",
                               round(this_followProt$secondSlope[1], digits = 2))
      
      followProt_max = max(c(followProt_original, followProt_smooth))
      
      
    }else{
      
      followProt_flag = F
      
      # followProt_xlab = ""
      # followProt_main = ""
      # followProt_max = orderProt_max
      # 
    }
    
    
    
    #######
    
    
    #### same row for the same phosphoSite 
    #### think about no rows scenario 
    
    ### order  this_orderPhospho
    
    # order_this_orderPhospho = this_orderPhospho%>%
    #   dplyr::arrange(avSlope)
    # ### think about how to identify phosphosites 
    
    ### take care of the ones appear in both temperatures
    join_order_follow = this_orderPhospho%>%
      dplyr::full_join(this_followPhospho, by = c("proteins",
                                                  "leadProtein",
                                                  "positionsOnProteins",
                                                  "residue"))%>%
      dplyr::arrange(avSlope.x)
    
    #### get any NA and plot the NA ones 
    
    # 
    # this_orderPhosphoNA = which(is.na(join_order_follow$probability.x))
    # this_followPhosphoNA = which(is.na(join_order_follow$probability.y))
    # 
    this_onlyFollow = join_order_follow%>%
      dplyr::filter(is.na(probability.x))
    this_onlyOrder = join_order_follow%>%
      dplyr::filter(is.na(probability.y))
    
    
    ### plot Order first 
    
    if(nrow(this_onlyOrder)>0)
    {
      for(k in 1:nrow(this_onlyOrder))
      {
        
        ### make the name 
        
        cat(i, k , "\n")
        order_name = paste(this_onlyOrder$proteins[k], 
                           this_onlyOrder$positionsOnProteins[k],
                           this_onlyOrder$residue[k],
                           this_onlyOrder$pepS.x[k],
                           sep = "-")
        
        
        ##### plot them
        
        order_original = as.numeric(unlist(strsplit(this_onlyOrder$original.x[k], split = "_")))
        order_smooth = as.numeric(unlist(strsplit(this_onlyOrder$smooth.x[k], split = "_")))
        
        order_data1 = order_original[or1]
        order_curve1 = order_smooth[or1]
        
        order_data2 = order_original[or2]
        order_curve2 = order_smooth[or2]
        
        order_min = min(c(order_original, order_smooth))
        order_max = max(c(order_original, order_smooth))
        
        order_xlab = paste0("Doses    Blue Rep1_", round(this_onlyOrder$firstSlope.x[k], digits = 2),
                            "    Red Rep2_",
                            round(this_onlyOrder$secondSlope.x[k], digits = 2))
        
        
        ## for the protein 
        plot(NULL, xlab = orderProt_xlab, ylab = paste0("prot abundance ", order_tag), 
             xlim = c(1,dosages), ylim = c(0, orderProt_max),
             main = paste0(this_protein,"_",round(this_orderProt$avSlope[1], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(order_rep1, orderProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(order_rep1, orderProt_curve1, type = "l", col = "blue")
        lines(order_rep2, orderProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(order_rep2, orderProt_curve2, type = "l", col = "red", pch = 16)
        
        
        
        plot(NULL, xlab = order_xlab, ylab = paste0("ptm abundance ", order_tag), 
             xlim = c(1,dosages), ylim = c(0, order_max),
             main = paste0(order_name,"_",round(this_onlyOrder$avSlope.x[k], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(order_rep1, order_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(order_rep1, order_curve1, type = "l", col = "blue")
        lines(order_rep2, order_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(order_rep2, order_curve2, type = "l", col = "red")
        
        ## for the protein 
        
        
        
        if(followProt_flag == T)
        {
          plot(NULL, xlab = followProt_xlab, ylab = paste0("prot abundance ", follow_tag), 
               xlim = c(1,dosages), ylim = c(0, followProt_max),
               main = paste0(this_protein,"_",round(this_followProt$avSlope[1], digits = 2)),
               cex.main = 0.8, cex.lab = 0.8)
          
          lines(follow_rep1, followProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
          lines(follow_rep1, followProt_curve1, type = "l", col = "blue")
          lines(follow_rep2, followProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
          lines(follow_rep2, followProt_curve2, type = "l", col = "red", pch = 16)
        }else{
          
          
          plot(NULL, xlab = "", ylab = paste0("prot abundance ", follow_tag), 
               xlim = c(1,dosages), ylim = c(0, orderProt_max),
               main = "",
               cex.main = 0.8, cex.lab = 0.8)
          
        }
        
        
        ### for the follow site 
        
        plot(NULL, xlab = "", ylab = paste0("ptm abundance ", follow_tag), 
             xlim = c(1,dosages), ylim = c(0, order_max),
             main = "",
             cex.main = 0.8, cex.lab = 0.8)
        
        
        
        
        
      }
      
      
    }
    
    
    #### for follow only 
    
    if(nrow(this_onlyFollow)>0)
    {
      for(k in 1:nrow(this_onlyFollow))
      {
        
        ### make the name 
        
        cat(i, k , "\n")
        follow_name = paste(this_onlyFollow$proteins[k], 
                            this_onlyFollow$positionsOnProteins[k],
                            this_onlyFollow$residue[k],
                            this_onlyFollow$pepS.y[k],
                            sep = "-")
        
        
        ##### plot them
        
        follow_original = as.numeric(unlist(strsplit(this_onlyFollow$original.y[k], split = "_")))
        follow_smooth = as.numeric(unlist(strsplit(this_onlyFollow$smooth.y[k], split = "_")))
        
        follow_data1 = follow_original[fr1]
        follow_curve1 = follow_smooth[fr1]
        
        follow_data2 = follow_original[fr2]
        follow_curve2 = follow_smooth[fr2]
        
        follow_min = min(c(follow_original, follow_smooth))
        follow_max = max(c(follow_original, follow_smooth))
        
        follow_xlab = paste0("Doses    Blue Rep1_", round(this_onlyFollow$firstSlope.y[k], digits = 2),
                             "    Red Rep2_",
                             round(this_onlyFollow$secondSlope.y[k], digits = 2))
        
        
        
        ## for the protein 
        plot(NULL, xlab = orderProt_xlab, ylab = paste0("prot abundance ", order_tag), 
             xlim = c(1,dosages), ylim = c(0, orderProt_max),
             main = paste0(this_protein,"_",round(this_orderProt$avSlope[1], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(order_rep1, orderProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(order_rep1, orderProt_curve1, type = "l", col = "blue")
        lines(order_rep2, orderProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(order_rep2, orderProt_curve2, type = "l", col = "red", pch = 16)
        
        plot(NULL, xlab = "", ylab = paste0("ptm abundance ", order_tag), 
             xlim = c(1,dosages), ylim = c(0, follow_max),
             main = "",
             cex.main = 0.8, cex.lab = 0.8)
        
        
        ## for the protein 
        
        if(followProt_flag == T)
        {
          plot(NULL, xlab = followProt_xlab, ylab = paste0("prot abundance ", follow_tag), 
               xlim = c(1,dosages), ylim = c(0, followProt_max),
               main = paste0(this_protein,"_",round(this_followProt$avSlope[1], digits = 2)),
               cex.main = 0.8, cex.lab = 0.8)
          
          lines(follow_rep1, followProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
          lines(follow_rep1, followProt_curve1, type = "l", col = "blue")
          lines(follow_rep2, followProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
          lines(follow_rep2, followProt_curve2, type = "l", col = "red", pch = 16)
          
        }else{
          
          
          plot(NULL, xlab = "", ylab = paste0("prot abundance ", follow_tag), 
               xlim = c(1,dosages), ylim = c(0, orderProt_max),
               main = "",
               cex.main = 0.8, cex.lab = 0.8)
          
          
          
        }
        
        
        ### for the follow site 
        
        plot(NULL, xlab = follow_xlab, ylab = paste0("ptm abundance ", follow_tag), 
             xlim = c(1,dosages), ylim = c(0, follow_max),
             main = paste0(follow_name,"_",round(this_onlyFollow$avSlope.y[k], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(follow_rep1, follow_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(follow_rep1, follow_curve1, type = "l", col = "blue")
        lines(follow_rep2, follow_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(follow_rep2, follow_curve2, type = "l", col = "red")
        
        
        
        
        
      }
      
      
    }
    
    ### no phosphorylation sites for either
    if((!nrow(this_onlyOrder)>0)&(!nrow(this_onlyFollow)>0))
    {
      cat(i, "\n")
      
      ## for the protein 
      plot(NULL, xlab = orderProt_xlab, ylab = paste0("prot abundance ", order_tag), 
           xlim = c(1,dosages), ylim = c(0, orderProt_max),
           main = paste0(this_protein,"_",round(this_orderProt$avSlope[1], digits = 2)),
           cex.main = 0.8, cex.lab = 0.8)
      
      lines(order_rep1, orderProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
      lines(order_rep1, orderProt_curve1, type = "l", col = "blue")
      lines(order_rep2, orderProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
      lines(order_rep2, orderProt_curve2, type = "l", col = "red", pch = 16)
      
      
      plot(NULL, xlab = "", ylab = paste0("ptm abundance ", order_tag), 
           xlim = c(1,dosages), ylim = c(0, orderProt_max),
           main = "",
           cex.main = 0.8, cex.lab = 0.8)
      
      
      ## for the protein 
      
      if(followProt_flag == T)
      {
        plot(NULL, xlab = followProt_xlab, ylab = paste0("prot abundance ", follow_tag), 
             xlim = c(1,dosages), ylim = c(0, followProt_max),
             main = paste0(this_protein,"_",round(this_followProt$avSlope[1], digits = 2)),
             cex.main = 0.8, cex.lab = 0.8)
        
        lines(follow_rep1, followProt_data1, cex = 0.8, type = "p", col = "blue", pch = 16)
        lines(follow_rep1, followProt_curve1, type = "l", col = "blue")
        lines(follow_rep2, followProt_data2, cex = 0.8, type = "p", col = "red", pch = 16)
        lines(follow_rep2, followProt_curve2, type = "l", col = "red", pch = 16)
        
      }else{
        
        plot(NULL, xlab = "", ylab = paste0("prot abundance ", follow_tag), 
             xlim = c(1,dosages), ylim = c(0, orderProt_max),
             main = "",
             cex.main = 0.8, cex.lab = 0.8)
        
        
      }
      
      
      ### for the follow site 
      
      plot(NULL, xlab = "", ylab = paste0("ptm abundance ", follow_tag), 
           xlim = c(1,dosages), ylim = c(0, followProt_max),
           main = "",
           cex.main = 0.8, cex.lab = 0.8)
      
      
    }
    
    
    #dev.off()
    
  }
  
  
  dev.off()
  
  
  
  
}





#### scatter plot for two proteins 




joinScatterOutput_curve = function(protFit1,
                             protFit2,
                             slope_upper = 0.5,
                             slope_lower = -0.5,
                             outputName,
                             pdfName)

{
  # ### use slop as fold change 
 
  # 
  # protFit1 = double_fitPep_s1_output
  # protFit2 = double_fitPep_s2_output
  # slope_upper = 0.5
  # slope_lower = -0.5
  # outputName = paste0(output_dir, "prot_2s.tsv")
  # pdfName = paste0(output_dir, "prot_2s.pdf") 
  # 
  
  
  
   # 
  join_prot = protFit1%>%
    dplyr::left_join(protFit2, by = "proteins")%>%
    dplyr::mutate(uniprot = uniprot.x, protName = protName.x, geneName = geneName.x, symbol = symbol.x)%>%
    dplyr::select(uniprot, protName, geneName, symbol,proteins,response.x, response.y, avSlope.x, avSlope.y)%>%
    na.omit()

  write.table(join_prot, outputName, quote = F, row.names = F, sep = "\t")

  
  
  join_prot_plot = join_prot%>%
    dplyr::group_by(uniprot, protName, geneName,response.x, response.y, avSlope.x, avSlope.y)%>%
    dplyr::arrange(symbol)%>%
    dplyr::slice(1)%>%
    dplyr::mutate(thermal = case_when(avSlope.y >= slope_upper | avSlope.y <= slope_lower ~ "hit",
                                      T ~ "nonhit"))%>%
    dplyr::mutate(abundance = case_when(avSlope.x >= slope_upper | avSlope.x <= slope_lower ~ "hit",
                                        T ~ "nonhit"))%>%
    dplyr::mutate (mark = case_when(thermal == "hit" & abundance == "hit" ~ "both",
                                    thermal == "hit" & abundance == "nonhit" ~ "thermal",
                                    thermal == "nonhit" & abundance == "hit" ~ "abundance",
                                    T ~ "none"))
    
    
    
    ggplot(join_prot_plot, aes(avSlope.x, avSlope.y))+
    geom_point(aes(color = mark),
               size = 1,
               alpha = 0.6) +
    scale_colour_manual(values=c("darkgreen", "darkorange", "grey", "navyblue"), name = "", labels = c("abundance", "both", "unaffected", "thermal"))+
    geom_hline(yintercept = slope_upper, colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) + 
    geom_hline(yintercept = slope_lower, colour = "black", linetype = "dashed", size = 0.3,alpha = 0.5) +
    geom_vline(xintercept = slope_upper, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) + 
    geom_vline(xintercept = slope_lower, colour = "black", linetype="dashed", size = 0.3,alpha = 0.5) +
    xlim(min(join_prot$avSlope.x*1.01), max(join_prot$avSlope.x*1.01))+
    ylim(min(join_prot$avSlope.y*1.01), max(join_prot$avSlope.y*1.01)) +
    geom_label_repel(data = join_prot_plot[which(join_prot_plot$mark!="none"),],
                     aes(label = symbol,color = mark),
                     size = 2,
                     #alpha = 0.8,
                     label.size = 0.1,
                     label.r = 0.1,
                     segment.size = 0.3,
                     segment.alpha = 0.3,
                     force = 2,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"),
                     seed = 123)+
      theme_bw()+
    theme(text = element_text(size=15)) +
    labs(x="abundance \n(average slope at 37C)", y = "thermal stability \n(average slope at 52C)")
  
  
  ggsave(pdfName)
  
  
    
  ###  plot scatter 
  ### have the ones marked with gene name 
  
  ### low/control
  # pdf(pdfName, useDingbats = F)
  # 
  # plot(join_prot$avSlope.x, join_prot$avSlope.y,
  #      type = "p", pch = 16, cex = 0.3,  #### size of the dots and shape
  #      xlim = c(-1.5,1.5), ### range of the plot
  #      ylim = c(-1.5,1.5),
  #      xlab = "avSlope_37C",
  #      ylab = "avSlope_52C",
  #      main = "average slope")
  # abline(h= slope_upper, lty = 3)
  # abline(h= slope_lower, lty = 3)
  # abline(v= slope_upper, lty = 3)
  # abline(v= slope_lower, lty = 3)
  # abline(h = 0, lty = 2)
  # abline(v = 0, lty = 2)
  # 
  # 
  # ### add marks to the fc ones 
  # 
  # join_hit = join_prot%>%
  #   dplyr::filter(avSlope.x>slope_upper| avSlope.x<slope_lower | avSlope.y>slope_upper| avSlope.y<slope_lower)
  # 
  # 
  # 
  # text(jitter(join_hit$avSlope.x), jitter(join_hit$avSlope.y), labels = join_hit$symbol, pos = 4, cex = 0.25)
  # 
  # 
  # 
  # dev.off()
  # 
  # 
  # 
  #### I gonna change this in a second after I modifiy the ggplot output in the other file 
  
  
  
  return(join_prot)
  
  
  
}


























