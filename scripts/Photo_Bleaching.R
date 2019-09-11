#Title: Bleaching estimates from pictures
#Project: Porites Nutrition June 2019
#Author: HM Putnam
#Edited by: KH Wong
#Date Last Modified: 20190710
#See Readme file for details

rm(list=ls()) #clears workspace 

if ("vegan" %in% rownames(installed.packages()) == 'FALSE') install.packages('vegan') 
if ("ggpubr" %in% rownames(installed.packages()) == 'FALSE') install.packages('ggpubr') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("lsmeans" %in% rownames(installed.packages()) == 'FALSE') install.packages('lsmeans') 
if ("multcompView" %in% rownames(installed.packages()) == 'FALSE') install.packages('multcompView') 


#Read in required libraries
##### Include Versions of libraries
library("vegan")
library("ggpubr")
library("gridExtra")
library("plyr") 
library("lsmeans")
library("multcompView")
library("ggplot2")


#Required Data files

Data <- read.csv("data/2019/Color_Score/C/olor_Score.csv", header=T, sep=",", na.string="NA") #read in data file
Sample.Info <- read.csv("data/2019/Meta_Data/Master_Coral_Sheet.csv", header=T, sep=",", na.string="NA") #read in data file 
#Sample.Info <- Sample.Info[,c(7,5,9)]
#Tank.Info <- read.csv("Data/Tank_to_Treatment.csv", header=T, sep=",", na.string="NA") #read in data file 

Data <-na.omit(Data)

Data$Red.Norm.Coral <- Data$Red.Coral/Data$Red.Standard #normalize to color standard
Data$Green.Norm.Coral <- Data$Green.Coral/Data$Green.Standard #normalize to color standard
Data$Blue.Norm.Coral <- Data$Blue.Coral/Data$Blue.Standard #normalize to color standard

blch.scor <- as.matrix(cbind(Data$Red.Norm.Coral,Data$Green.Norm.Coral,Data$Blue.Norm.Coral)) #create matrix
rownames(blch.scor) <- Data$Coral.ID.TP #name columns in dataframe

dist <- vegdist(blch.scor, method="euclidean") #calculate distance matrix of color scores

PCA.color <- princomp(dist) #run principal components Analysis
summary(PCA.color) # view variance explained by PCs

Blch <- as.data.frame(PCA.color$scores[,1]) #extract PC1
#Blch$Coral.ID.TP <- rownames(blch.scor)
#Blch <- merge(Blch, Data, by="PLUG.ID")
Blch  <- cbind(Blch, Data$Fragment.ID) #make a dataframe of PC1 and experiment factors
colnames(Blch) <- c("Bleaching.Score", "Fragment.ID")

#Blch$SpGroup <- paste(Blch$Treatment, Blch$Species)
#Blch$Bleaching.Score <- Blch$`PCA.color$scores[, 1]`

Blch$Bleaching.Score <- -Blch$Bleaching.Score

# write.table(Blch,"~/MyProjects/Holobiont_Integration/RAnalysis/Output/Bleaching_Score.csv",sep=",", row.names=FALSE)

write.table(Blch, file = "data/2019/Color_Score/Bleaching_Score.csv", append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

#Attach metadata

Blch.meta <- merge(Blch, Sample.Info, by = "Fragment.ID")

Blch.meta2 <- Blch.meta %>%
  filter(Status == "good")

All.Means <- ddply(Blch.meta2, c('Origin', 'Transplant.site'), summarize,
                   mean= mean(Bleaching.Score, na.rm=T), #mean pnet
                   N = sum(!is.na(Bleaching.Score)), # sample size
                   se = sd(Bleaching.Score, na.rm=T)/sqrt(N)) #SE
All.Means

pd <- position_dodge(0.1)
legend.title <- "Treatment"

blch.plot <- ggplot(All.Means, aes(x=Transplant.site, y=mean, color=Origin, group=Origin)) + #set up plot information
  geom_errorbar(aes(x=Transplant.site, ymax=mean+se, ymin=mean-se), colour="black", width=.1, position = pd) + #add standard error bars about the mean
  geom_line(position=pd, color="black") +
  xlab("Transplant Site") + #label x axis
  ylab(expression("Bleaching Score")) + #label y axis
  #  ylim(-0.9, -0.1)+
  geom_point(position=pd, aes(fill=Origin), color ="black", pch=21, size=4)+
  scale_fill_manual(legend.title, values=c("#999999", "#000000"))+ #colour modification
  theme_bw() + theme(panel.border = element_rect(linetype = "solid", color = "black"), panel.grid.major = element_blank(), #Makes background theme white
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(face = 'bold', 
                                  size = 12, 
                                  hjust = 0),
        legend.position = c(0.78,0.1),
        legend.box.background = element_rect(colour = "black"))#set title attributes

blch.plot #view plot

ggsave(file="output/2019/2019_Photo_Bch.pdf", blch.plot, width = 8, height = 8, units = c("in"))
