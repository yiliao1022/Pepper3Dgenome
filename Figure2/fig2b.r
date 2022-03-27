library(ggplot2)
library(plyr)
setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/FINAL/All")
fre<-read.table("All.domain.length.40k.bed.bed",header = T)
fre
mu <- ddply(fre, "Tissues", summarise, grp.mean=mean(Length))

ggplot(fre, aes(x=Length,color=Tissues)) + 
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Tissues),
                             linetype="dashed") +
  xlim(10000, 1000000)


library(ggplot2)
library(plyr)
setwd("/home/yiliao/Documents/Pepper_2021/ABcompartments/draw/2022_revision/4papper")
fre<-read.table("All.length.10k.bed.bed",header = T)
fre
mu <- ddply(fre, "Tissues", summarise, grp.mean=mean(Length))

ggplot(fre, aes(x=Length,color=Tissues)) + 
  geom_density() +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Tissues),
             linetype="dashed") +
  xlim(9000, 1000000)
