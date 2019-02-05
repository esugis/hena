# This script plots the statistics of the interactions in the integrated dataset
# The plot will be save to the folder /results/integration/

load("/results/integration/integrated_int.RData")

library(ggplot2)

integrated_int$data_source <- factor(integrated_int$data_source, levels = c("PBA","IAH", "ADIA","SIA", #PPI
"ADN", "ABA_CA1", "ABA_CA2","ABA_CA3","ABA_CA4","ABA_DG","ABA_Sptn","ABA_subiculum", #co-expression
 "TGEN", "HBTRC","ADNI_VER", # epistasis +IGRI
"ADNI_ADAS11_slope" , "ADNI_ADAS13_slope",                
"ADNI_CDRSB_slope" , "ADNI_EcogPtDivatt_slope", "ADNI_EcogPtLang_slope",            
"ADNI_EcogPtMem_slope","ADNI_EcogPtPlan_slope","ADNI_EcogPtTotal_slope",           
"ADNI_EcogPtVisspat_slope","ADNI_EcogSPDivatt_slope","ADNI_EcogSPLang_slope",            
"ADNI_EcogSPMem_slope", "ADNI_EcogSPOrgan_slope", "ADNI_EcogSPPlan_slope",            
"ADNI_EcogSPTotal_slope", "ADNI_EcogSPVisspat_slope","ADNI_FAQ_slope",                   
"ADNI_MMSE_slope" , "ADNI_MOCA_slope","ADNI_RAVLT_forgetting_slope",      
"ADNI_RAVLT_immediate_slope" ,"ADNI_RAVLT_learning_slope", "ADNI_RAVLT_perc_forgetting_slope", 
"ADNI_latest_ADAS11", "ADNI_latest_ADAS13", "ADNI_latest_CDRSB_1E-8",           
"ADNI_latest_EcogPtDivatt","ADNI_latest_EcogPtLang" , "ADNI_latest_EcogPtMem",            
"ADNI_latest_EcogPtOrgan" , "ADNI_latest_EcogPtPlan" , "ADNI_latest_EcogPtTotal",          
"ADNI_latest_EcogPtVisspat","ADNI_latest_EcogSPDivatt" , "ADNI_latest_EcogSPLang",           
"ADNI_latest_EcogSPMem", "ADNI_latest_EcogSPOrgan" ,"ADNI_latest_EcogSPPlan",           
"ADNI_latest_EcogSPTotal", "ADNI_latest_EcogSPVisspat" ,  "ADNI_latest_FAQ_1E-8",             
"ADNI_latest_MMSE", "ADNI_latest_MOCA" , "ADNI_latest_RAVLT_forgetting" ,    
"ADNI_latest_RAVLT_immediate", "ADNI_latest_RAVLT_learning", "ADNI_latest_RAVLT_perc_forgetting"))
                                                                                         

integrated_int_plot <- integrated_int[,4:5]


pdf("results/integration/integrated_int_stats_stacked_hist.pdf")
ggplot(data = integrated_int_plot, aes(x = data_source,fill = interaction_type)) +
  geom_histogram(alpha=0.7, stat="count") +
  stat_count( geom="text", colour="white", size=2, angle=90,
              aes(label=..count.., group=interaction_type, y=0.5*(..count..)))+
  theme_bw()+
  labs(x = "Data source")+
  labs(y= "Number of Interactions")+
  theme(axis.text.y= element_blank(),
        axis.text.x = element_text(angle=90, hjust=1,vjust=-0.001,size = 6),
        #axis.text.y = element_text(angle=90),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.y = element_blank())+
  guides(fill=guide_legend(title="Interaction type"))+ scale_y_continuous(trans='log2')

dev.off()


# Plot table containing #nodes #interactions for each datasource
data_summary<-integrated_int[1, -c(4:5)]
colnames(data_summary)<-c("ds_ID", "nodes", "interactions")
data_summary<-unname(data_summary)
dsor<-unique(integrated_int$data_source)
for(i in 1:length(dsor)){
    data_s<-as.character(dsor[i])
    ds<-integrated_int[integrated_int$data_source%in%data_s,]
    int<-dim(ds)[1]
    node<-length(unique(c(ds$ensg.A,ds$ensg.B)))
    data_summary<-rbind(data_summary, c(data_s, node, int))
}
data_summary<- data_summary[-1,]
colnames(data_summary)<-c("ds_ID", "nodes", "interactions")
write.table(data_summary, file="results/integration/integrated_int_stats.txt", sep="\t", quote=F, row.names=F)
