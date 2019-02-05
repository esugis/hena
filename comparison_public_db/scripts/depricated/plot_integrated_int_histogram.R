# Create dataset
#integrated_int <- data.frame(
#  int_type =  c( rep('coexpression', 40), rep('PPI', 30), rep('epistasis', 14), rep('IGRI', 6)) ,
#  data_source = c( rep('ADN', 10), rep('CA1_ABA', 10), rep('CA2_ABA', 10), rep('DG_ABA', 10), rep("IAH",10), rep("ADIA",10),rep("SIA",10),rep('ADNI', 14), rep('ADNI', 6) )
#)
library(dplyr)
print("Loading integrated dataset")

load(file="../results/integration/integrated_int.RData")

data <- integrated_int%>%group_by(interaction_type, data_source) %>% count

data <- data.frame(data, stringsAsFactors = T)

data$interaction_type <- as.factor(data$interaction_type)
#
# plot circulat hist
# library
library(tidyverse)

# Set a number of 'empty bar' to add at the end of each group
empty_bar=2
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$interaction_type), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$interaction_type=rep(levels(data$interaction_type), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(interaction_type)
data$id=seq(1, nrow(data))

# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)


# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=log2(n), fill=interaction_type)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=log2(n)+10, label=paste(data_source, n, sep="  "), hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) 

p

ggsave("integrated_int_stats_simple_log2.pdf", plot = p, device = "pdf", path = "results/comparisons/",
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)



### Plot with lines

# Create dataset
#data_raw <- data.frame(
#  interaction_type =  c( rep('coexpression', 40), rep('PPI', 30), rep('epistasis', 14), rep('IGRI', 6)) ,
#  data_source = c( rep('ADN', 10), rep('CA1_ABA', 10), rep('CA2_ABA', 10), rep('DG_ABA', 10), rep("IAH",10), rep("ADIA",10),rep("SIA",10),rep('ADNI', 14), rep('ADNI', 6) ) 
#)

library(dplyr)

#data <- integrated_int %>% 
#  group_by(interaction_type, data_source) %>%
#  count

#data <- data.frame(data, stringsAsFactors = F)



# Set a number of 'empty bar' to add at the end of each group
empty_bar=1
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$interaction_type), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$interaction_type=rep(levels(data$interaction_type), each=empty_bar)
data=rbind(data, to_add)
data=data %>% arrange(interaction_type)
data$id=seq(1, nrow(data))

# Get the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)


base_data=data %>% 
  group_by(interaction_type) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))
#!!!! Change start and  end upon nesessity, it will draw aproper group line

base_data[1,4]<- 9
base_data[4,4]<- 113

# Make the plot
p = ggplot(data, aes(x=as.factor(id), y=log2(n), fill=interaction_type)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=log2(n), fill=interaction_type), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=log2(n)+10, label= paste(data_source, n, sep="  "), hjust=hjust), color="black", fontface="bold",alpha=0.6, size=1.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
 # geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -10, label=interaction_type), hjust=c(1,1,0,0), colour = "black", alpha=0.8, size= 1.8, fontface="bold", inherit.aes = FALSE)

p

ggsave("integrated_int_stats_groups_log2.pdf", plot = p, device = "pdf", path = "results/comparisons/",
       scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       dpi = 300, limitsize = TRUE)


# Stacked histogram
load("results/integration/integrated_int.RData")

library(ggplot2)

ds <- unique(integrated_int$data_source)


integrated_int_plot <- integrated_int[,4:5]
rm(integrated_int)
integrated_int_plot$data_source <- factor(integrated_int_plot$data_source, levels = c("PBA","IAH", "ADIA","SIA", #PPI
"ADN", "ABA_CA1", "ABA_CA2","ABA_CA3","ABA_CA4","ABA_DG","ABA_Sptn","ABA_subiculum", #co-expression
"TGEN", "HBTRC", "ADNI_VER", # epistasis +IGRI
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


pdf("integrated_int_stats_stacked_hist.pdf")
ggplot(data = integrated_int_plot, aes(x = data_source,fill = interaction_type)) +
geom_histogram(alpha=0.7, stat="count") +
stat_count( geom="text", colour="white", size=2, angle=90,
aes(label=..count.., group=interaction_type, y=0.5*(..count..)))+
theme_bw()+
labs(x = "Data source")+
labs(y= "Number of Interactions")+
theme(axis.text.y= element_blank(),
axis.text.x = element_text(angle=90, hjust=1,size = 6),
#axis.text.y = element_text(angle=90),
panel.grid.major.x = element_blank(),
panel.grid.minor.x = element_blank(),
axis.ticks.y = element_blank())+
guides(fill=guide_legend(title="Interaction type"))+ scale_y_continuous(trans='log2')

dev.off()

