r_time_course_FC = function(df, 
                          genelist, 
                          plot_conditions,
                          ctrl_condition,
                          logscale = TRUE,
                          plot_errorbar = FALSE,
                          plot_title,
                          plotdir,
                          match_time_norm,
                          align_yaxis,
                          df_to_use)
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
  # plotdir: the directory where the plots will be saved
  # match_time_norm: (TRUE or FALSE) whether the normalization will be done by each time point of the ctrl_condition
  # align_yaxis: (TRUE or FALSE) whether the y axes will be aligned across all plots
  # df_to_use: type of intensities plotted (i.e. raw, iBAQ, or LFQ), this is only for labeling the plots
  # -----------------------------
  
  # 0: Set up
  library(gtools)
  library(ggrepel)
  library(egg)
  library(dplyr)
  library(data.table)
  library(ggsci)
  library(ggplot2)
  library(tidyverse)
  library(stringr)
  library(MBQN)
  library(preprocessCore)
  library(gridExtra)
  library(matrixStats)
  library(limma)
  library(EnhancedVolcano)
  print("Plotting in progress :D The R script used to make this plot was adapted from Elisa Holstein's work")
  
  # 1: remove all irrelevant genes and conditions ---------
  df_gene <- subset(df, row.names(df) %in% genelist)  # This would also remove NA gene names
  vector <- colnames(df_gene)
  if (ctrl_condition %in% plot_conditions)
    all_conditions <- plot_conditions
  else
    all_conditions <- append(plot_conditions, ctrl_condition)
  df_gene_cond <- df_gene %>%
    dplyr::select(contains(all_conditions))
  # get the timepoint vector (the second to last position in column names)
  time_vector <- c()
  for (sample in colnames(df_gene_cond))
  {t = tail(unlist(strsplit(sample, "_")), 2)[1]
    time_vector = append(time_vector, t)
  }
  time_vector = sort(unique(time_vector))
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
  df_all_norm <- data.frame("time"= c(), "condition"= c(), "gene"= c(), "logFC"= c())
  # log-transform the dataframe if it has not been transformed yet
  if (logscale == FALSE)
    df_gene_cond = log(df_gene_cond, base = 2)
  # if match
  
  # compute the fold change for each time point and add the result to df_all_norm
  starting_timepoint = time_vector[1]
  df_t = df_gene_cond[grepl(paste("_", starting_timepoint, "_", sep=""), colnames(df_gene_cond))]
  ctrl_col = rowMeans(df_t[grepl(ctrl_condition,colnames(df_t))], na.rm = TRUE)
  for (t in time_vector)
  {
  timepoint = paste("_",t, "_", sep="")
  df_t = df_gene_cond[grepl(timepoint, colnames(df_gene_cond))]
  #print(df_t)
  if (match_time_norm)
    ctrl_col = rowMeans(df_t[grepl(ctrl_condition,colnames(df_t))], na.rm = TRUE)
  df_norm <- df_t - ctrl_col # Normalization: subtract columns from the mean of the ctrl one
  #print(df_norm)
  for (condition in all_conditions)
  {samples_per_condition = grepl(condition,colnames(df_norm))
  
  # give warning if the condition was not measured at the timepoint
  if (sum(samples_per_condition) == 0)
  {print(paste("Warning: cannot find data for ", condition, " at timepoint ", t, sep=""))
  
  # if at the starting timepoint: set protein level to be similar to that of ctrl_condition (i.e. fold change to be all 0)
  if (t == time_vector[1])
  {
    print(paste("The protein level of this condition was set to that of the normalizing sample (aka control): ", ctrl_condition, sep= ""))
    logFC = rep(0, nrow(df_norm))
    df_per_condition = data.frame("time"= rep(t,length(logFC)),
                                  "condition"= rep(condition, length(logFC)), 
                                  "gene"= rep(row.names(df_norm)),
                                  "logFC"= logFC)
    df_all_norm = rbind(df_all_norm, df_per_condition)
  }}
  else
    {logFC = unlist(df_norm[samples_per_condition])
  
  df_per_condition = data.frame("time"= rep(t,length(logFC)),
                                "condition"= rep(condition, length(logFC)), 
                                "gene"= rep(row.names(df_norm),sum(samples_per_condition)),
                                "logFC"= logFC)
  df_all_norm = rbind(df_all_norm, df_per_condition)}
  }
  }
  df_all_norm$time <- gsub("\\D", "", df_all_norm$time) # remove all non-digits inthe column time
  df_all_norm$time <- as.numeric(df_all_norm$time)
  df_all_norm = na.omit(df_all_norm) # remove all NAs
  
  # 3 plot with ggplot2 -------------------------------
  plot_FC <- ggplot(df_all_norm, aes(x=time, y=logFC, color=condition)) 
  if (plot_errorbar)
    plot_FC = plot_FC +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  size=1) +
    stat_summary(fun = "mean", geom="point",  size=2)+
    stat_summary(fun.data = "mean_se", geom = "errorbar", size=1, 
                 width = length(unique(df_all_norm$time))/8,
                 fun.args = list(mult=1))
  else
    plot_FC = plot_FC + geom_point(aes(fill=condition), shape=16, size=2) +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  size=1)
  # adjust plot appearance
  plot_FC <-  plot_FC + 
    theme_article() +
    labs(title=paste("Fold change of", df_to_use, "intensities with respect to", ctrl_condition),
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
    #scale_colour_manual(values = blue4) +
    scale_x_continuous(breaks = unique(df_all_norm$time), 
                       limits = c(min(unique(df_all_norm$time)), max(unique(df_all_norm$time))))
  # align the y axis if align_yaxis = TRUE
  if (align_yaxis)
    plot_FC = plot_FC + 
    scale_y_continuous(limits = c(min(df_all_norm$logFC), max(df_all_norm$logFC)))
  # 4 return the plot object ---------------------------
  plot_title = paste(plot_title, ".pdf", sep="")
  plot_height = (length(unique(df_all_norm$gene))%/%4 + 1)*2 +2.5
  ggsave(plot_title, plot_FC, path = plotdir, width = 13, height = plot_height, 
         limitsize = FALSE,units = 'in')
  print("Plotting done! Check the prompt for any warning")
}

r_time_course_intensity = function(df, 
                            genelist, 
                            plot_conditions,
                            logscale = TRUE,
                            plot_errorbar = FALSE,
                            plot_title,
                            plotdir,
                            align_yaxis,
                            df_to_use)
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
  # df_to_use: type of intensities plotted (i.e. raw, iBAQ, or LFQ), this is only for labeling the plots
  # -----------------------------
  
  # 0: Set up
  library(gtools)
  library(ggrepel)
  library(egg)
  library(dplyr)
  library(data.table)
  library(ggsci)
  library(ggplot2)
  library(tidyverse)
  library(stringr)
  library(MBQN)
  library(preprocessCore)
  library(gridExtra)
  library(matrixStats)
  library(limma)
  library(EnhancedVolcano)
  print("Plotting in progress :D The R script used to make this plot was adapted from Elisa Holstein's work")
  
  # 1: remove all irrelevant genes and conditions ---------
  df_gene <- subset(df, row.names(df) %in% genelist)  # This would also remove NA gene names
  df_gene_cond <- df_gene %>%
    dplyr::select(contains(plot_conditions))
  # get the timepoint vector (the second to last position in column names)
  time_vector <- c()
  for (sample in colnames(df_gene_cond))
  {t = tail(unlist(strsplit(sample, "_")), 2)[1]
  time_vector = append(time_vector, t)
  }
  time_vector = sort(unique(time_vector))
  
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
    df_t = df_gene_cond[grepl(timepoint, colnames(df_gene_cond))]
    
    for (condition in plot_conditions)
    {samples_per_condition = grepl(condition,colnames(df_t))
     
    # give warning if the condition was not measured at the timepoint
     if (sum(samples_per_condition) == 0)
      print(paste("Warning: cannot find data for ", condition, " at timepoint ", t, sep=""))
    else
    { # otherwise, add data to df_all
    intensity = unlist(df_t[samples_per_condition])
    df_per_condition = data.frame("time"= rep(t,length(intensity)),
                                  "condition"= rep(condition, length(intensity)), 
                                  "gene"= rep(row.names(df_t),sum(samples_per_condition)),
                                  "intensity"= intensity)
    df_all = rbind(df_all, df_per_condition)}
    }
  }
  df_all$time <- gsub("\\D", "", df_all$time) # remove all non-digits inthe column time
  df_all$time <- as.numeric(df_all$time)
  df_all = na.omit(df_all) # remove all NAs
  
  # 3 plot with ggplot2 -------------------------------
  plot_intensity <- ggplot(df_all, aes(x=time, y=intensity, color=condition)) 
  if (plot_errorbar)
    plot_intensity = plot_intensity +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  size=1) +
    stat_summary(fun = "mean", geom="point",  size=2)+
    stat_summary(fun.data = "mean_se", geom = "errorbar", size=1, 
                 width = length(unique(df_all$time))/8,
                 fun.args = list(mult=1))
  else
    plot_intensity = plot_intensity + geom_point(aes(fill=condition), shape=16, size=2) +
    facet_wrap( ~ gene, ncol=4, scales = "free") +
    stat_summary(fun = "mean", geom="line",  size=1)
  # adjust plot appearance
  plot_intensity <-  plot_intensity + 
    theme_article() +
    labs(title=paste(df_to_use, "intensities of selected genes"),
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
    #scale_colour_manual(values = blue4) +
    scale_x_continuous(breaks = unique(df_all$time), 
                       limits = c(min(unique(df_all$time)), max(unique(df_all$time))))
  # align the y axis if align_yaxis = TRUE
  if (align_yaxis)
    plot_intensity = plot_intensity + 
    scale_y_continuous(limits = c(min(df_all$intensity), max(df_all$intensity)))
  
  # 4 return the plot object ---------------------------
  plot_title = paste(plot_title, ".pdf", sep="")
  plot_height = (length(unique(df_all$gene))%/%4+1)*2 +2.5
  ggsave(plot_title, plot_intensity, path = plotdir, width = 13, height = plot_height, 
         limitsize = FALSE,units = 'in')
  
}
