#########################
  #METABOANALYST CODE
#########################

########  SET UP  ##############

#libraries
library(MetaboAnalystR)
library(fitdistrplus) #got error message first time running line 28 bc didn't have this package
library(RJSONIO) #same here 

#CUSTOM ADDUCTS
add.vec <- c("M+FA-H [1-]","M-H [1-]","2M-H [1-]","M-H+O [1-]","M(C13)-H [1-]",
             "2M+FA-H [1-]","M-3H [3-]","M-2H [2-]","M+ACN-H [1-]",
             "M+HCOO [1-]","M+CH3COO [1-]","M-H2O-H [1-]","M [1+]","M+H [1+]",
             "M+2H [2+]","M+3H [3+]","M+H2O+H [1+]","M-H2O+H [1+]",
             "M(C13)+H [1+]","M(C13)+2H [2+]","M(C13)+3H [3+]","M-NH3+H [1+]",
             "M+ACN+H [1+]","M+ACN+2H [2+]","M+2ACN+2H [2+]","M+3ACN+2H [2+]",
             "M+NH4 [1+]","M+H+NH4 [2+]","2M+H [1+]","2M+ACN+H [1+]")


##################################

#4 PFAS MIXTURE 

#create object for storing data
mSet1<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet1<- SetPeakFormat(mSet1, "rmp") #ranking by p-value 
mSet1<-UpdateInstrumentParameters(mSet1, 5.0, "mixed", "no"); #not enforcing primary ions (includes adducts we dont have)
mSet1<-Read.PeakListData(mSet1, "C:/users/efleury/Downloads/QG_comp_full_METABOLOMICS_4PFAS.csv");
mSet1<-SanityCheckMummichogData(mSet1)


#map selected adducts to current data 
mSet1<-Setup.AdductData(mSet1, add.vec); #custom adducts 
mSet1<-PerformAdductMapping(mSet1, "mixed") #positive and negative adducts

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet1<-SetPeakEnrichMethod(mSet1, "mum", "v2")
mSet1<-SetMummichogPval(mSet1, .05) #parameter sweep: 0.05, 0.1, 0.2--selected 0.05 

#the next step takes three or four minutes to run 
mSet1<-PerformPSEA(mSet1, "hsa_mfn", "current", 3 , 10000) 


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_4pfas<- as.data.frame(mSet1$mummi.resmat)

#store as csv
write.csv(mummi_results_4pfas, "mummi_4pfas_res.csv")

###Individual PFAS

#PFOA
#create object for storing data
mSet2<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet2<- SetPeakFormat(mSet2, "rmp")
mSet2<-UpdateInstrumentParameters(mSet2, 5.0, "mixed", "no");
mSet2<-Read.PeakListData(mSet2, "C:/users/efleury/Downloads/PFOA_METABOANALYST.csv");
mSet2<-SanityCheckMummichogData(mSet2)

#map selected adducts to current data 
mSet2<-Setup.AdductData(mSet2, add.vec);
mSet2<-PerformAdductMapping(mSet2, "mixed")

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet2<-SetPeakEnrichMethod(mSet2, "mum", "v2")
mSet2<-SetMummichogPval(mSet2, .05) #pval

#the next step takes three or four minutes to run 
mSet2<-PerformPSEA(mSet2, "hsa_mfn", "current", 3 , 10000)


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_pfoa<- as.data.frame(mSet2$mummi.resmat)

#store as csv
write.csv(mummi_results_pfoa, "mummi_pfoa_res.csv")

##########PFOS

#create object for storing data
mSet3<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet3<- SetPeakFormat(mSet3, "rmp")
mSet3<-UpdateInstrumentParameters(mSet3, 5.0, "mixed", "no");
mSet3<-Read.PeakListData(mSet3, "C:/users/efleury/Downloads/PFOS_METABOANALYST.csv");
mSet3<-SanityCheckMummichogData(mSet3)

#map selected adducts to current data 
mSet3<-Setup.AdductData(mSet3, add.vec);
mSet3<-PerformAdductMapping(mSet3, "mixed")

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet3<-SetPeakEnrichMethod(mSet3, "mum", "v2")
mSet3<-SetMummichogPval(mSet3, .05) #pval

#the next step takes three or four minutes to run 
mSet3<-PerformPSEA(mSet3, "hsa_mfn", "current", 3 , 10000)


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_pfos<- as.data.frame(mSet3$mummi.resmat)

#store as csv
write.csv(mummi_results_pfos, "mummi_pfos_res.csv")

##########PFNA


#create object for storing data
mSet4<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet4<- SetPeakFormat(mSet4, "rmp")
mSet4<-UpdateInstrumentParameters(mSet4, 5.0, "mixed", "no");
mSet4<-Read.PeakListData(mSet4, "C:/users/efleury/Downloads/PFNA_METABOANALYST.csv");
mSet4<-SanityCheckMummichogData(mSet4)

#map selected adducts to current data 
mSet4<-Setup.AdductData(mSet4, add.vec);
mSet4<-PerformAdductMapping(mSet4, "mixed")

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet4<-SetPeakEnrichMethod(mSet4, "mum", "v2")
mSet4<-SetMummichogPval(mSet4, .05) #pval

#the next step takes three or four minutes to run 
mSet4<-PerformPSEA(mSet4, "hsa_mfn", "current", 3 , 10000)


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_pfna<- as.data.frame(mSet4$mummi.resmat)

#store as csv
write.csv(mummi_results_pfna, "mummi_pfna_res.csv")

##########PFHxS


#create object for storing data
mSet5<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet5<- SetPeakFormat(mSet5, "rmp")
mSet5<-UpdateInstrumentParameters(mSet5, 5.0, "mixed", "no");
mSet5<-Read.PeakListData(mSet5, "C:/users/efleury/Downloads/PFHxS_METABOANALYST.csv");
mSet5<-SanityCheckMummichogData(mSet5)

#map selected adducts to current data 
mSet5<-Setup.AdductData(mSet5, add.vec);
mSet5<-PerformAdductMapping(mSet5, "mixed")

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet5<-SetPeakEnrichMethod(mSet5, "mum", "v2")
mSet5<-SetMummichogPval(mSet5, .05) #pval

#the next step takes three or four minutes to run 
mSet5<-PerformPSEA(mSet5, "hsa_mfn", "current", 3 , 10000)


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_pfhxs<- as.data.frame(mSet5$mummi.resmat)

#store as csv
write.csv(mummi_results_pfhxs, "mummi_pfhxs_res.csv")

#####MeFOSAA

#create object for storing data
mSet6<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet6<- SetPeakFormat(mSet6, "rmp")
mSet6<-UpdateInstrumentParameters(mSet6, 5.0, "mixed", "no");
mSet6<-Read.PeakListData(mSet6, "C:/users/efleury/Downloads/ME_PFOSA_ACOH_METABOANALYST.csv");
mSet6<-SanityCheckMummichogData(mSet6)

#map selected adducts to current data 
mSet6<-Setup.AdductData(mSet6, add.vec);
mSet6<-PerformAdductMapping(mSet6, "mixed")

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet6<-SetPeakEnrichMethod(mSet6, "mum", "v2")
mSet6<-SetMummichogPval(mSet6, .05) #pval

#the next step takes three or four minutes to run 
mSet6<-PerformPSEA(mSet6, "hsa_mfn", "current", 3 , 10000)


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_mefosaa<- as.data.frame(mSet6$mummi.resmat)

#store as csv
write.csv(mummi_results_mefosaa, "mummi_mefosaa_res.csv")

######### 5 PFAS MIXTURE

#create object for storing data
mSet7<-InitDataObjects("mass_all", "mummichog", FALSE)

#set peak format 
mSet7<- SetPeakFormat(mSet7, "rmp")
mSet7<-UpdateInstrumentParameters(mSet7, 5.0, "mixed", "no");
mSet7<-Read.PeakListData(mSet7, "C:/users/efleury/Downloads/QG_comp_full_METABOLOMICS_4PFAS_and_ME_PFOSA.csv");
mSet7<-SanityCheckMummichogData(mSet7)

#map selected adducts to current data 
mSet7<-Setup.AdductData(mSet7, add.vec);
mSet7<-PerformAdductMapping(mSet7, "mixed")

#perform mummichog algorithm using selected adducts, using version 2 of the mummichog algorithm
mSet7<-SetPeakEnrichMethod(mSet7, "mum", "v2")
mSet7<-SetMummichogPval(mSet7, .05) #pval

#the next step takes three or four minutes to run 
mSet7<-PerformPSEA(mSet7, "hsa_mfn", "current", 3 , 10000)


#below line of code gives you the plot that you see on the web browswer
#mSet<- PlotPeaks2Paths(mSet), "peaks_to_paths_0_", "png", 72, width=NA)

#store the results as a data frame so you can plot with ggplot 
mummi_results_five<- as.data.frame(mSet7$mummi.resmat)

#store as csv
write.csv(mummi_results_five, "mummi_five_res.csv")




##########PLOT###########
# MORE DETAILED CODE FOUND IN NEW_PLOTS.RMD THIS IS INCLUDED HERE TO CHECK THE 
# PARAMETERS

library(ggplot2)
library(tidyverse)
library(dplyr)
library(tibble)

#read in data if not running first part 
#results<- read.csv("mummi_4pfas_res.csv")

#get data in correct format 

#function to filter for rows with gamma less than 0.05

filter_gamma_less_than_05<- function(df){
  df_filtered<-df%>%
    filter(Gamma < 0.05)
  return(df_filtered)
}

mummi_for_plot<- filter_gamma_less_than_05(results)

mummi_for_plot<- rownames_to_column(mummi_for_plot, "pathway_name")#convert rownames to a column

#PLOT: x= log10(pval), y=pathway name, size= significant hits/expected hits 

ggplot(mummi_for_plot, aes(x = -log10(Gamma), y = pathway_name, size = Hits.sig / Expected)) +
  geom_point() + 
  scale_size_continuous(guide = guide_legend(title = "Enrichment Factor")) +
  labs(x = "-log10(p-value)", 
       y = "Metabolic Pathways", 
       title = "Metabolic Pathways Enrichment Analysis", 
       subtitle = "Data from both HILIC positive and C18 negative modes") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", 
        plot.title = element_text( face = "bold", size = 14, hjust = 0.5), 
        plot.subtitle = element_text( size = 12, hjust = 0.5),
        axis.title.x = element_text( size = 12),
        axis.title.y = element_text( size = 12),
        legend.title = element_text( size = 10),
        legend.text = element_text( size = 10),
        legend.box = "vertical") 



