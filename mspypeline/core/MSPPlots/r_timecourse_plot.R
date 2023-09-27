r_time_course_FC = function(df, 
                          genelist,
                          genelist_name,
                          plot_conditions,
                          ctrl_condition,
                          logscale = TRUE,
                          plot_errorbar = FALSE,
                          plot_title,
                          savedir,
                          match_time_norm,
                          align_yaxis,
                          ifFDR = FALSE,
                          df_to_use,
                          selected_normalizer)
{# R script for plotting time-course FC data. Modified from Elisa Holstein's script (S18_E27_TimeCourse on 2/21/2023)
 # Parameters, all are received from Python (core -> MSPPlots -> BasePlotter.py -> plot_r_timecourse):
  # ---------------------------------
  # df: dataframe of protein intensities, passed from python. The intensities (raw, iBAQ, or LFQ) were log2-transformed
  #     Indexes of df are gene names
  # genelist: list of genes whose protein levels will be plotted
  # plot_conditions: list of all conditions to be plotted
  # ctrl_condition: condition to normalize data against
  # plot_errorbar: Option to whether to plot error bars instead of dots
  # plot_title: name of the plot pdf file
  # savedir: the directory where the plots will be saved
  # match_time_norm: (TRUE or FALSE) whether the normalization will be done by each time point of the ctrl_condition
  # align_yaxis: (TRUE or FALSE) whether the y axes will be aligned across all plots
  # ifFDR: should p-values should be adjusted (Benjamini-Hochberg), to be passed to plot_significance()
  #For naming the output file:
  # df_to_use: type of intensities plotted (i.e. raw, iBAQ, or LFQ), this is only for labeling the plots
  # selected_normalizer: which normalizer was used
  # -----------------------------
  
  # 0: Set up
  library(gtools)
  library(ggrepel)
  library(egg)
  library(dplyr)
  library(data.table)
  #library(ggsci)
  library(ggplot2)
  library(tidyverse)
  library(stringr)
  #library(MBQN)
  #library(preprocessCore)
  library(gridExtra)
  library(matrixStats)
  library(limma)
  #library(EnhancedVolcano)
  print("Plotting in progress :D The R script used to make this plot was adapted from Elisa Holstein's work")
  
  # 1: remove all irrelevant genes and conditions ---------
  df_gene <- subset(df, row.names(df) %in% genelist)  # This would also remove NA gene names
  vector <- colnames(df_gene)
  if (ctrl_condition %in% plot_conditions)
    all_conditions <- plot_conditions
  else
    all_conditions <- append(plot_conditions, ctrl_condition)
  ctrl_condition = paste0(ctrl_condition, "_")
  all_conditions = unlist(lapply(all_conditions, function(x) paste0(x,"_")))
  df_gene_cond <- df_gene %>%
    dplyr::select(starts_with(all_conditions))
  # get the timepoint vector (the second to last position in column names)
  time_vector <- c()
  for (sample in colnames(df_gene_cond))
  {t = tail(unlist(strsplit(sample, "_")), 2)[1]
    time_vector = append(time_vector, t)
  }
  time_vector = unique(time_vector)
  time_unit = ""
  while ((is.na(as.numeric(t))) && (t != ""))
  {
    n = nchar(t)
    time_unit = paste(substr(t,n,n), time_unit, sep = "")
    t = substr(t,1, n-1)
  }
  if (time_unit == "")
    {time_unit = "unit not found"
    time_vector = sort(as.numeric(time_vector))
    time_vector = as.character(time_vector)}
   else
    {time_vector = gsub(time_unit, "", time_vector)
    time_vector = sort(as.numeric(time_vector))
    time_vector = lapply(time_vector, function(x) paste0(x,time_unit))}
  # 2 Normalization -----------------------------------
  # Mean and standard dev of log2 intensities are stored in 2 dataframes: all_mean_norm and all_sd_norm
  # after processing, columns of both dataframes are named by the genes, plus two columns indicating time and condition
  
  # start an empty dataframe
  df_all_norm <- data.frame("time"= c(), "condition"= c(), "gene"= c(), "logFC"= c())
  # log-transform the dataframe if it has not been transformed yet
  if (logscale == FALSE)
    df_gene_cond = log(df_gene_cond, base = 2)
  # if match
  
  # compute the fold change for each time point and add the result to df_all_norm
  starting_timepoint = time_vector[1]
  df_t = df_gene_cond %>%
    dplyr::select(contains(paste0("_", starting_timepoint, "_"))) %>%
    dplyr::select(starts_with(ctrl_condition))
  ctrl_col = rowMeans(df_t, na.rm = TRUE)
  for (t in time_vector)
  {
  timepoint = paste("_",t, "_", sep="")
  df_t = df_gene_cond %>%
    dplyr::select(contains(timepoint))
  if (match_time_norm)
    ctrl_col = rowMeans(df_t %>% dplyr::select(starts_with(ctrl_condition)), na.rm = TRUE)
  df_norm <- df_t - ctrl_col # Normalization: subtract columns from the mean of the ctrl one
  for (condition in all_conditions)
  {df_per_condition = df_norm %>%
    dplyr::select(starts_with(condition))
  # give warning if the condition was not measured at the timepoint
  if (ncol(df_per_condition) == 0)
  {print(paste("Warning: cannot find data for ", condition, " at timepoint ", t, sep=""))
  # if at the starting timepoint: set protein level to be similar to that of ctrl_condition (i.e. fold change to be all 0)
  if (t == time_vector[1])
  {
    print(paste("The protein level of this condition was set to that of the normalizing sample (aka control): ", ctrl_condition, sep= ""))
    logFC = rep(0, nrow(df_per_condition))
    df_per_condition_to_all = data.frame("time"= rep(t,length(logFC)),
                                  "condition"= rep(condition, length(logFC)), 
                                  "gene"= rep(row.names(df_per_condition)),
                                  "logFC"= logFC)
    df_all_norm = rbind(df_all_norm, df_per_condition_to_all)
  }}
  else
    {logFC = unlist(df_per_condition)
  
  df_per_condition_to_all = data.frame("time"= rep(t,length(logFC)),
                                "condition"= rep(condition, length(logFC)), 
                                "gene"= rep(row.names(df_norm),ncol(df_per_condition)),
                                "logFC"= logFC)
  df_all_norm = rbind(df_all_norm, df_per_condition_to_all)}
  }
  }
  if (time_unit != "unit not found")
    df_all_norm$time <- gsub(time_unit, "", df_all_norm$time) # remove all non-digits inthe column time
  df_all_norm$time <- as.numeric(df_all_norm$time)
  df_all_norm = na.omit(df_all_norm) # remove all NAs
  
  # 3 plot with ggplot2 -------------------------------
  plot_FC <- ggplot(df_all_norm, aes(x=time, y=logFC, color=condition)) 
  if (plot_errorbar)
    plot_FC = plot_FC +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  linewidth=1) +
    stat_summary(fun = "mean", geom="point",  size=2)+
    stat_summary(fun.data = "mean_se", geom = "errorbar", linewidth=0.5, 
                 width = 1,
                 fun.args = list(mult=1))
  else
    plot_FC = plot_FC + geom_point(aes(fill=condition), shape=16, size=2) +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  linewidth=1)
  # adjust plot appearance
  plot_FC <-  plot_FC + 
    theme_article() +
    labs(title=paste("Fold change of", df_to_use, "intensities with respect to", ctrl_condition),
         subtitle = paste("Gene list:", genelist_name),
         x =paste("time (",time_unit,")",sep=""), y = "log2(FC)") +
    theme(plot.title = element_text(hjust = 0.5, color="black", size=18, face="bold"), 
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text=element_text(size=13, color = "black"),
          axis.line = element_line(colour="black"),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position = "top",
          # strip.text.x = element_blank(), # if you want to leave out x headers for each plot
          strip.text.x = element_text(size = 14, face = "bold")) + # remove axis labels of facet plots
      scale_x_continuous(limits = c(min(unique(df_all_norm$time)), max(unique(df_all_norm$time))))
  # align the y axis if align_yaxis = TRUE
  if (align_yaxis)
    plot_FC = plot_FC + 
    scale_y_continuous(limits = c(min(df_all_norm$logFC), max(df_all_norm$logFC)))
  # 4 return the plot object ---------------------------
  plot_title = paste(plot_title, ".pdf", sep="")
  plot_height = (length(unique(df_all_norm$gene))%/%4 + 1)*2 +3
  plotdir = paste(savedir, "/plots", sep="")  # directory to save plot
  # Handle exception for when the directory is too long
  if (nchar(paste0(plotdir,'/', plot_title)) > 250)
  {print("Warning: directory to the time-course plot exceeds 250 characters")
    plot_title1 = paste0(genelist_name,' ',df_to_use, '.pdf')
    if (nchar(paste0(plotdir,'/', plot_title1)) > 250)
      print("Cannot save the plot. Please copy the data file to a shorter directory and try again")
    else
      {print(paste("The time-course plot will be saved under the name '",plot_title1,"'"))
      print("To make sure no plot is overwritten, please copy the data file to a shorter directory and try again")
      plot_title = plot_title1
    }
  }
  ggsave(plot_title, plot_FC, path = plotdir, width = 13, height = plot_height, 
         limitsize = FALSE,units = 'in', dpi = 300)
  
  # 5 pass the dataframe to plot_significance for statistical testing
  plot_significance(df_all_norm, savedir, genelist_name, match_time_norm, df_to_use, selected_normalizer, ifFDR)
  print("Plotting done! Check the prompt for any warning")
}

r_time_course_intensity = function(df, 
                            genelist, 
                            genelist_name,
                            plot_conditions,
                            logscale = TRUE,
                            plot_errorbar = FALSE,
                            plot_title,
                            savedir,
                            align_yaxis,
                            ifFDR = FALSE,
                            df_to_use,
                            selected_normalizer)
{# R script for plotting time-course protein intensities (no fold change). Modified from Elisa Holstein's script (S18_E27_TimeCourse on 2/21/2023)
  # Parameters, all are received from Python (core -> MSPPlots -> BasePlotter.py -> plot_r_timecourse):
  # ---------------------------------
  # df: dataframe of protein intensities, passed from python. The intensities (raw, iBAQ, or LFQ) were log2-transformed
  #     Indexes of df are gene names
  # genelist: list of genes whose protein levels will be plotted
  # plot_conditions: list of all conditions to be plotted
  # plot_errorbar: Option to whether to plot error bars instead of dots
  # plot_title: name of the plot pdf file
  # plotdir: the directory where the plots will be saved
  # align_yaxis: (TRUE or FALSE) whether the y axes will be aligned across all plots
  # ifFDR: should p-values be adjusted for multiple comparisons (Benjamini-Hochberg), to be passed to plot_significance()
  #For naming the output file:
  # genelist_name: name of the genelist
  # df_to_use: type of intensities plotted (i.e. raw, iBAQ, or LFQ), this is only for labeling the plots
  # selected_normalizer: which normalizer was used
  # -----------------------------
  
  # 0: Set up
  library(gtools)
  #library(ggrepel)
  library(egg)
  library(dplyr)
  library(data.table)
  #library(ggsci)
  library(ggplot2)
  library(tidyverse)
  library(stringr)
  #library(MBQN)
  #library(preprocessCore)
  library(gridExtra)
  library(matrixStats)
  library(limma)
  #library(EnhancedVolcano)
  print("Plotting in progress :D The R script used to make this plot was adapted from Elisa Holstein's work")
  
  # 1: remove all irrelevant genes and conditions ---------
  df_gene <- subset(df, row.names(df) %in% genelist)  # This would also remove NA gene names
  plot_conditions = unlist(lapply(plot_conditions, function(x) paste0(x,"_")))
  df_gene_cond <- df_gene %>%
    dplyr::select(starts_with(plot_conditions))
  # get the timepoint vector (the second to last position in column names)
  time_vector <- c()
  for (sample in colnames(df_gene_cond))
  {t = tail(unlist(strsplit(sample, "_")), 2)[1]
  time_vector = append(time_vector, t)
  }
  time_vector = unique(time_vector)
  
  # get the time unit
  time_unit = ""
  while ((is.na(as.numeric(t))) && (t != ""))
  {
    n = nchar(t)
    time_unit = paste(substr(t,n,n), time_unit, sep = "")
    t = substr(t,1, n-1)
  }
  if (time_unit == "")
    time_unit = "unit not found"
  
  
  # 2 Normalization -----------------------------------
  # Mean and standard dev of log2 intensities are stored in 2 dataframes: all_mean_norm and all_sd_norm
  # after processing, columns of both dataframes are named by the genes, plus two columns indicating time and condition
  
  # start an empty dataframe
  df_all <- data.frame("time"= c(), "condition"= c(), "gene"= c(), "intensity"= c())
  # log-transform the dataframe if it has not been transformed yet
  if (logscale == FALSE)
    df_gene_cond = log(df_gene_cond, base = 2)
  # compute the fold change for each time point and add the result to df_all_norm
  for (t in time_vector)
  {
    timepoint = paste("_",t, "_", sep="")
    df_t = df_gene_cond %>%
      dplyr::select(contains(timepoint))
    for (condition in plot_conditions)
    {df_per_condition = df_t %>%
        dplyr::select(starts_with(condition))
    # give warning if the condition was not measured at the timepoint
     if (ncol(df_per_condition) == 0)
      print(paste("Warning: cannot find data for ", condition, " at timepoint ", t, sep=""))
    else
    { # otherwise, add data to df_all
    intensity = unlist(df_per_condition)
    df_per_condition_to_all = data.frame("time"= rep(t,length(intensity)),
                                  "condition"= rep(condition, length(intensity)), 
                                  "gene"= rep(row.names(df_t),ncol(df_per_condition)),
                                  "intensity"= intensity)
    df_all = rbind(df_all, df_per_condition_to_all)}
    }
  }
  if (time_unit != "unit not found")
    df_all$time <- gsub(time_unit, "", df_all$time) # remove all non-digits inthe column time
  df_all$time <- as.numeric(df_all$time)
  df_all = na.omit(df_all) # remove all NAs
  
  # 3 plot with ggplot2 -------------------------------
  plot_intensity <- ggplot(df_all, aes(x=time, y=intensity, color=condition)) 
  if (plot_errorbar)
    plot_intensity = plot_intensity +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  linewidth=1) +
    stat_summary(fun = "mean", geom="point",  size=2)+
    stat_summary(fun.data = "mean_se", geom = "errorbar", linewidth=0.5, 
                 width = 1,
                 fun.args = list(mult=1))
  else
    plot_intensity = plot_intensity + geom_point(aes(fill=condition), shape=16, size=2) +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  linewidth=1)
  # adjust plot appearance
  plot_intensity <-  plot_intensity + 
    theme_article() +
    labs(title=paste(df_to_use, "intensities of selected genes"),
         subtitle = paste("Gene list:", genelist_name),
         x =paste("time (",time_unit,")",sep=""), y = "log2(intensity)") +
    theme(plot.title = element_text(hjust = 0.5, color="black", size=18, face="bold"), 
          axis.title.x = element_text(size=14, face="bold"),
          axis.title.y = element_text(size=14, face="bold"),
          axis.text=element_text(size=13, color = "black"),
          axis.line = element_line(colour="black"),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.position = "top",
          # strip.text.x = element_blank(), # if you want to leave out x headers for each plot
          strip.text.x = element_text(size = 14, face = "bold")) + # remove axis labels of facet plots
  scale_x_continuous(limits = c(min(unique(df_all$time)), max(unique(df_all$time))))
  # align the y axis if align_yaxis = TRUE
  if (align_yaxis)
    plot_intensity = plot_intensity + 
    scale_y_continuous(limits = c(min(df_all$intensity), max(df_all$intensity)))
  
  # 4 return the plot object ---------------------------
  plot_title = paste(plot_title, ".pdf", sep="")
  plot_height = (length(unique(df_all$gene))%/%4+1)*2 +3
  plotdir = paste(savedir, "/plots", sep="")  # directory to save plot
  # Handle exception for when the directory is too long
  if (nchar(paste0(plotdir, '/',plot_title)) > 250)
  {print("Warning: directory to the time-course plot exceeds 250 characters")
    plot_title1 = paste0(genelist_name,' ',df_to_use, '.pdf')
    if (nchar(paste0(plotdir, '/', plot_title1)) > 250)
      print("Cannot save the plot. Please copy the data file to a shorter directory and try again")
    else
    {print(paste("The time-course plot will be saved under the name '",plot_title1,"'"))
      print("To make sure no plot is overwritten, please copy the data file to a shorter directory and try again")
      plot_title = plot_title1
    }
  }
  ggsave(plot_title, plot_intensity, path = plotdir, width = 13, height = plot_height, 
         limitsize = FALSE,units = 'in', dpi = 300)
  
  # 5 pass the dataframe to plot_significance for statistical testing
  plot_significance(df_all, savedir, genelist_name, match_time_norm = FALSE, df_to_use, selected_normalizer, ifFDR)  # match_time_norm set to FALSE since no normalization was done
  print("Plotting done! Check the prompt for any warning")
}


plot_significance = function(df, savedir, genelist_name, match_time_norm, df_to_use, selected_normalizer, ifFDR = FALSE)
{# perform statistical testing of selected conditions and genes. Return heatmaps and a csv file showing the p-value
  # of the testing for all genes
  # For each gene: use a two-way ANOVA with time and condition (aka treatment/cell line/...) as factors,
  # to see if there is a significant difference among conditions. Then use Tukey Honestly-Significant-Difference
  # (TukeyHSD) for pair-wise comparisons. p-values from the TukeyHSD tests of condition are plotted and exported
  # as heatmaps and csv file, respectively
  
  # Parameters----------------------------------------------------
  # df: processed dataframes of fold change/intensities. Passed from r_timecourse_FC or r_timecourse_intensity
  # savedir: need to add either "/plots" or "/csv_significance" to become directory for saving results
  # match_time_norm: whether data was normalized per time point
  # ifFDR: should the p-values be corrected for multiple testings (Benjamini-Hochberg)
  # For naming the output file:
  # df_to_use: type of intensities plotted (i.e. raw, iBAQ, or LFQ), this is only for labeling the plots
  # selected_normalizer: which normalizer was used
  # ----------------------------------------------------
  
  # 0 ------------------------------------------------------
  # Rename df columns so data processing is easier
  library(scales)
  library(ggtext)  # for coloring axes' text
  colnames(df) = c("time", "condition", "gene", "intensity")
  if (length(unique(df$condition)) == 1)
    {print('Statistical testing is skipped since only one condition was selected')
    return()}
  if (match_time_norm)
    print("Warning: normalization per timepoint was chosen. Statistical testings might not reflect true significant differences between conditions")
  df$time = sapply(df$time, function(x) paste(x))
  all_conditions = sort(unique(df$condition))
  condition_code = c()
  for (i in 1:length(all_conditions))
    condition_code[all_conditions[i]] = paste(i)
  df$condition = gsub("-", "@", df$condition)
  # make a modified condition list to fill missing values
  all_conditions_mod = gsub("-", "@", all_conditions)
  condition_code_mod = c()
  for (i in 1:length(all_conditions_mod))
    condition_code_mod[all_conditions_mod[i]] = paste(i)
  # 1  ---------------------------------------------------------------
  all_sign_df = data.frame("sample1" = c(), "sample2" = c(), "gene" = c(), "p_value" = c())
  all_genes = unique(df$gene)
  for (g in all_genes)
  {# Two-way ANOVA testing per gene
    df_gene = df %>%
      subset(gene == g) %>%
      group_by(time, condition) %>%
      filter(n()>1)
    if ((length(unique(df_gene$time))<2) || (length(unique(df_gene$condition))<2))
    {
      print(paste(genelist_name, ": skip statistical testing for gene", g, "due to insufficient observations"))
      df_gene = data.frame("sample1" = all_conditions_mod[-1],
                           "sample2" = rep(all_conditions_mod[1], length(all_conditions_mod)-1),
                           "gene" = rep(g, length(all_conditions_mod)-1),
                           "p_value" = rep(NA, length(all_conditions_mod)-1))
    } else
    {two.way = aov(intensity ~ time + condition, data = df_gene)
    tukey = TukeyHSD(two.way)
    tukey = tukey$condition
    df_gene = data.frame("sample1" = sapply(row.names(tukey), function(x) unlist(strsplit(x, "-"))[1]),
                         "sample2" = sapply(row.names(tukey), function(x) unlist(strsplit(x, "-"))[2]),
                         "gene" = rep(g, length(row.names(tukey))),
                         "p_value" = tukey[,4])}
    # add rows with NA p-value if not all comparisons were made so that color-coding the axes are correct
    sample1_list = unique(df_gene$sample1)
    sample2_list = unique(df_gene$sample2)
    if (length(unique(sample2_list))< (length(all_conditions_mod)-1))
      for (sample in all_conditions_mod[-length(all_conditions_mod)])
        if (!(sample %in% sample2_list))
          df_gene = rbind(df_gene, data.frame("sample1"= sample1_list, "sample2"= rep(sample, length(sample1_list)),
                                            "gene"= rep(g, length(sample1_list)),
                                            "p_value"=rep(NA, length(sample1_list))))
    # add df_gene to the dataframe containing everything
    all_sign_df = rbind(all_sign_df, df_gene)
  }
  all_sign_df$sample1 = gsub("@", "-", all_sign_df$sample1)
  all_sign_df$sample2 = gsub("@", "-", all_sign_df$sample2)
  all_sign_df$label_sample1 = sapply(all_sign_df$sample1, function(x) condition_code[x])
  all_sign_df$label_sample2 = sapply(all_sign_df$sample2, function(x) condition_code[x])
  # If ifFDR is TRUE: correct the p-values for multiple hypothesis testings
  if (ifFDR)
    all_sign_df$p_value = p.adjust(all_sign_df$p_value, method = "BH")
  all_sign_df$sign_level <- cut(all_sign_df$p_value, breaks=c(-1,0.00001, 0.0001, 0.001, 0.01, 0.05, 1), 
                         label=c("*****","****" ,"***", "**", "*", ""))  # Create column of significance labels
  
  # 2 save all_sign_df to a csv file -------------------------------------------------
  # assemble file name
  file_title = paste("significance ",genelist_name, "samples ")
  for (sample in all_conditions)
    file_title = paste(file_title, sample, ", ", sep="")
  if (match_time_norm)
    file_title = paste(file_title, "per time point,", sep="")
  if (selected_normalizer != "None")
    file_title = paste(file_title, selected_normalizer)
  file_title = paste(file_title, df_to_use, sep="")
  if (ifFDR)
    file_title = paste0(file_title, "_adjusted_pval")
  # save all_sign_df as a csv file
  csvdir = paste(savedir, "/csv_significance/", file_title,".csv", sep= "")
  # handle exception when directory is too long
  if (nchar(csvdir)>250)
  {print("Warning: directory to csv file of p-values exceeds 250 characters.")
    # set the filename to be shorter, including only "significance" and genelist_name
    csvdir1 = paste(savedir, "/csv_significance/significance ", genelist_name,' ',df_to_use,".csv", sep= "")
    if (nchar(csvdir)>250)
      print("Cannot save the csv file. Please copy the data file to a shorter directory and try again")
    else
    {print(paste("csv file of significance will be saved under the name: significance ",genelist_name,' ',df_to_use ,".csv", sep= ""))
      print("To make sure no file is overwritten, please copy the data file to a shorter directory and try again")
      csvdir = csvdir1}
  }
  write.csv(all_sign_df, file = csvdir, row.names = FALSE)
  
  # 3 plot the heatmap ------------------------------------------------------------------
  # add 10^-7 if p-value = 0 so it appears red on the plot
  all_sign_df[which(all_sign_df$p_value == 0),"p_value"] = 1e-7
  # sample legend
  sample_legend = "Sample labels: "
  for (i in 1:length(attributes(condition_code)$names))
    sample_legend = paste(sample_legend, condition_code[i], "=", attributes(condition_code)$names[i],", ", sep="")
  # significance legend
  sign_legend = "p-value: ***** < 0.00001 < **** < 0.0001 < *** < 0.001 < ** < 0.01 < * < 0.05"                               
  # legend of the test used
  test_legend = "Statistical testing done by two-way ANOVA for each gene followed by Tukey Honestly-Significant-Difference"
  
  # color hex code
  color_code = hue_pal()(length(all_conditions))
  # add color code of sample1 (y axis) to all_sign_df
  label_sample1_color = rep(NA, nrow(all_sign_df))
  label_sample2_color = rep(NA, nrow(all_sign_df))
  sample_labels = c()
  for (condition in 1:length(condition_code))
  {p1 = which(all_sign_df$label_sample1 == condition_code[condition])
    p2 = which(all_sign_df$label_sample2 == condition_code[condition])
    color_text = paste0("<span style = 'color:", color_code[condition], "'>", condition, "</span>")  # encoding the color of the text by ggtext
    label_sample1_color[p1] = color_text
    label_sample2_color[p2] = color_text
    sample_labels = append(sample_labels, color_text)
  }
  all_sign_df$label_sample1 = label_sample1_color
  all_sign_df$label_sample2 = label_sample2_color
  # Rearrange the samples in the heatmap
  all_sign_df$label_sample1 = as.factor(all_sign_df$label_sample1)
  all_sign_df$label_sample1 = factor(all_sign_df$label_sample1, levels = sample_labels[-1])
  all_sign_df$label_sample2 = as.factor(all_sign_df$label_sample2)
  all_sign_df$label_sample2 = factor(all_sign_df$label_sample2, levels = sample_labels[-length(sample_labels)])
  
  # Plot!
  plot_sign = ggplot(all_sign_df, aes(x= label_sample2, y= label_sample1))+
    geom_tile(aes(fill= p_value), color= "white")+
    geom_text(aes(label=sign_level), color="black", size=5) + 
    facet_wrap(~ gene, ncol=4, scales = "free_x")
  
  # adjust appearance
  plot_sign = plot_sign +
    theme_article() +
    labs(title=paste("p-values from pair-wise comparisons of selected conditions,", df_to_use,"intensities"),
         subtitle = paste("Gene list:", genelist_name),
         y = "",
         x = "Sample labels (explanation at the top of the plot)") +
    theme(plot.title = element_text(hjust = 0.5, color="black", size=18, face="bold"), 
          axis.title.x = element_text(size=15, face= "bold"),
          axis.title.y = element_text(size=20, face= "bold"),
          axis.text.y = element_markdown(size=15, face = 'bold'),
          axis.text.x = element_markdown(size=15, face = 'bold'),
          axis.line = element_line(colour="black"),
          legend.text = element_text(size=14),
          legend.title = element_blank(),
          legend.key.width = unit(6, "cm"),
          legend.position = "top",
          # strip.text.x = element_blank(), # if you want to leave out x headers for each plot
          strip.text.x = element_text(size = 14, face = "bold")) + # remove axis labels of facet plots
    scale_x_discrete(breaks = sample_labels[-length(sample_labels)], expand = c( 0, 0))+
    scale_y_discrete(breaks = sample_labels[-1], expand = c(0,0))+
    scale_fill_gradient(high = "#FFFFE0", low = "red", na.value = "white", trans = "log",
                        breaks = c(0.00001, 0.0001, 0.001, 0.01, 0.05),
                        guide = guide_colorbar(ticks.colour = "black"))
  
  # add the legend
  plot_sign = grid.arrange(plot_sign, top = sample_legend)                        
  plot_sign = grid.arrange(plot_sign, top = sign_legend)
  plot_sign = grid.arrange(plot_sign, top = test_legend)
  # 3 save results -----------------------------------------------------------------
  plot_name = paste(file_title, ".pdf", sep="")
  # save heatmaps as a pdf file
  plotdir = paste(savedir, "/plots", sep= "")
  # handle exception when directory is too long
  if (nchar(paste0(plotdir,"/",plot_name))>250)
  {print(paste("Warning: directory to the plot of significant level (heatmap): ", plotdir,"/",plot_name,", exceeds 250 characters.", sep=""))
    # compile a shorter plot name, including only "significance" and genelist_name
    plot_name1 = paste0("significance ",genelist_name,".pdf")
    # check if the new directory is short enough
    if (nchar(paste0(plotdir,"/",plot_name1))>250)
      print("Cannot save the plot. Please copy the data file to a shorter directory and try again")
    else
    {print(paste("The significance plot will be saved under the name '",plot_name1,"'"))
        print("To make sure no plot is overwritten, please copy the data file to a shorter directory and try again")
        plot_name = plot_name1}
    }
  
  plot_height = (length(unique(all_sign_df$gene))%/%4+1)*3 + 3
  ggsave(filename=plot_name, 
         plot=plot_sign, path = plotdir, width = 12, height = plot_height, 
         limitsize = FALSE,units = 'in')
}
