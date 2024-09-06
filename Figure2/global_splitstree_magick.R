setwd("/Users/ben/Desktop/baseclear/networks")

##remotes::install_github("KlausVigo/phangorn",ref="2753573")
#remotes::install_github("KlausVigo/phangorn")
##remotes::install_github("KlausVigo/tanggle",ref="f6c93c8")
#remotes::install_github("KlausVigo/tanggle")
#install.packages("magick")

library(tanggle)
library(phangorn)
library(ggplot2)
library(ggrepel)
library(magick)

dutch_net <- read.nexus.networx("dutch_combined_filtered.splitstree.nex",splits=T)

# read in tiplabel metadata
tip_metadata <- read.delim("baseclear_metadata.csv", sep=";", header=TRUE,check.names=FALSE,stringsAsFactor=F)
tip_metadata <- tip_metadata[,c("Sample","Haplotype")]
#tip_metadata has some trouble
#first some have lowercase V
tip_metadata$Sample <- gsub("v","V",tip_metadata$Sample)
#then there is V273-05 in Nnet but not in tip_metadata, and V196-05 in tip_metadata but not Nnet
#tip_metadata$Sample[tip_metadata$Sample=="V196-05"] <- "V273-05"

p1 <- ggsplitnet(dutch_net,lwd=0.2) + labs(title="A) Dutch Dataset")
ggsave(filename = 'plot1_dutch.png', device = 'png',width=7,height=7)

colors <- c("wt R"="purple1","46"="dodgerblue3","92"="dodgerblue3","34"="dodgerblue3","wt"="orange2")
shapes <- c("wt R"=18,"46"=15,"92"=17,"34"=16,"wt"=18)

p2 <- ggsplitnet(dutch_net,lwd=0) %<+% tip_metadata + 
  geom_tippoint(aes(color=Haplotype,shape=Haplotype),size=3.5,alpha=0.55) +
  #geom_text_repel(aes(label=label),size=1,max.overlaps=100,max.time=10,box.padding=0.1)+
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  labs(title="A) Dutch Dataset")+
  theme(legend.position = c(0.1,0.85),
        rect=element_rect(fill=NA),
        panel.background = element_rect(fill=NA))
ggsave(filename = 'plot2_dutch.png', device = 'png', bg = 'transparent',width=7,height=7)

plot1 <- image_read('plot1_dutch.png')
plot2 <- image_read('plot2_dutch.png')

img <- c(plot1, plot2)

t <- image_mosaic(img)
image_browse(t)
image_write(t,"dutch_combined.png",quality=100)

#Now we want to read in the global samples
splits <- read.nexus.splits("global_combined_filtered.min4.stree6")
global_net <- as.networx(splits)

# read in tiplabel metadata
tip_metadata <- read.delim("/Users/ben/Downloads/baseclear_metadata.csv", sep=";", header=TRUE,check.names=FALSE,stringsAsFactor=F)
tip_metadata <- tip_metadata[,c("Sample","Haplotype")]
#tip_metadata has some trouble
#first some have lowercase V
tip_metadata$Sample <- gsub("v","V",tip_metadata$Sample)
#then there is V273-05 in Nnet but not in tip_metadata, and V196-05 in tip_metadata but not Nnet
#tip_metadata$Sample[tip_metadata$Sample=="V196-05"] <- "V273-05"

p1 <- ggsplitnet(global_net,lwd=0.05) +
  labs(title="B) Global Dataset")
  
ggsave(filename = 'plot1_global.png', device = 'png',width=7,height=7)

colors <- c("wt R"="purple1","46"="dodgerblue3","92"="dodgerblue3","34"="dodgerblue3","wt"="orange2")
shapes <- c("wt R"=18,"46"=15,"92"=17,"34"=16,"wt"=18)

p2 <- ggsplitnet(global_net,lwd=0) %<+% tip_metadata + 
  geom_tippoint(aes(color=Haplotype,shape=Haplotype),size=3,alpha=0.55) +
  scale_color_manual(values=colors)+
  scale_shape_manual(values=shapes)+
  labs(title="B) Global Dataset")+
  theme(legend.position = "none",
        rect=element_rect(fill=NA),
        panel.background = element_rect(fill=NA))
ggsave(filename = 'plot2_global.png', device = 'png', bg = 'transparent',width=7,height=7)

plot1 <- image_read('plot1_global.png')
plot2 <- image_read('plot2_global.png')

img <- c(plot1, plot2)

global_magick <- image_mosaic(img)
image_browse(global_magick)
image_write(global_magick,"combined_global.png",quality=100)


#now want to read in the subsetted plots----
#first the 50% sampling
tip_metadata <- read.delim("/Users/ben/Downloads/global_dataset_metadata.csv", sep="\t", header=TRUE,check.names=FALSE,stringsAsFactor=F)
tip_metadata <- tip_metadata[,c("accession_id","azole_resistance_sensitive")]
colors <- c("Resistant"="dodgerblue3","Susceptible"="orange2")

net_50 <- as.networx(read.nexus.splits("global_fifty_percent.min4.stree6"))

p1 <- ggsplitnet(net_50,lwd=0.05)+  labs(title="C) 50% Triazole Resistant")
ggsave(filename = 'fifty_percent_plot1.png', device = 'png',width=5,height=4)

p2 <- ggsplitnet(net_50,lwd=0) %<+% tip_metadata + 
  geom_tippoint(aes(color=azole_resistance_sensitive),size=3,alpha=0.55) +
  scale_color_manual(values=colors)+
  labs(title="C) 50% Triazole Resistant")+
  theme(legend.position = "none",
        rect=element_rect(fill=NA),
        panel.background = element_rect(fill=NA))
ggsave(filename = 'fifty_percent_plot2.png', device = 'png', bg = 'transparent',width=5,height=4)

plot1_50 <- image_read('fifty_percent_plot1.png')
plot2_50 <- image_read('fifty_percent_plot2.png')

t_50 <- image_mosaic(c(plot1_50, plot2_50))
image_write(t_50,"global_fifty_percent.png",quality=100)


#Now we read in the 15% subsample----

net_15 <- as.networx(read.nexus.splits("global_fifteen_percent.min4.stree6"))

p1 <- ggsplitnet(net_15,lwd=0.05) + labs(title="D) 15% Triazole Resistant")
ggsave(filename = 'fifteen_percent_plot1.png', device = 'png',width=5,height=4)

p2 <- ggsplitnet(net_15,lwd=0) %<+% tip_metadata + 
  geom_tippoint(aes(color=azole_resistance_sensitive),size=3,alpha=0.55) +
  scale_color_manual(values=colors)+
  labs(title="D) 15% Triazole Resistant")+
  theme(legend.position = "none",
        rect=element_rect(fill=NA),
        panel.background = element_rect(fill=NA))
ggsave(filename = 'fifteen_percent_plot2.png', device = 'png', bg = 'transparent',width=5,height=4)

plot1_15 <- image_read('fifteen_percent_plot1.png')
plot2_15 <- image_read('fifteen_percent_plot2.png')

t_15 <- image_mosaic(c(plot1_15, plot2_15))
image_write(t_15,"global_fifteen_percent.png",quality=100)

#now want to read in the 5% subsetted data----

net_5 <- as.networx(read.nexus.splits("global_five_percent.min4.stree6"))

p1 <- ggsplitnet(net_5,lwd=0.05) +   labs(title="E) 5% Triazole Resistant")
ggsave(filename = 'five_percent_plot1.png', device = 'png',width=5,height=4)

p2 <- ggsplitnet(net_5,lwd=0) %<+% tip_metadata + 
  geom_tippoint(aes(color=azole_resistance_sensitive),size=3,alpha=0.55) +
  scale_color_manual(values=colors)+
  labs(title="E) 5% Triazole Resistant")+
  theme(legend.position = c(0.15,0.85),
        legend.title = element_blank(),
        rect=element_rect(fill=NA),
        panel.background = element_rect(fill=NA))
ggsave(filename = 'five_percent_plot2.png', device = 'png', bg = 'transparent',width=5,height=4)

plot1_5 <- image_read('five_percent_plot1.png')
plot2_5 <- image_read('five_percent_plot2.png')

t_5 <- image_mosaic(c(plot1_5, plot2_5))
image_write(t_5,"global_five_percent.png",quality=100)
