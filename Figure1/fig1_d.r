library("ggplot2")
library("dplyr")
setwd("/home/yiliao/Documents/Pepper_2021/Final/3_Hicheatmaps")
distance <- read.table("All.tsv.bed.bed", header=T)

distance
distance$Tissues <- as.factor(distance$Tissues)
mydata <-distance %>% filter(Tissues=='HL2' | Tissues=='GR1' | Tissues=='TZ1' | Tissues=='YP1') %>% filter (Distance=='Ratio')
ggplot(mydata, aes(x=Tissues, y=Frequency, fill=factor(Distance))) + geom_boxplot(outlier.shape = NA) +labs(fill = "Distance")

mydata1 <-distance %>% filter(Tissues=='HL1' | Tissues=='GR1') %>% filter (Distance=='Ratio')
mydata1
wilcox.test(Frequency~Tissues, data = mydata1, paired = T)

mydata2 <-distance %>% filter(Tissues=='HL1' | Tissues=='TZ1') %>% filter (Distance=='Ratio')
mydata2
wilcox.test(Frequency~Tissues, data = mydata2, paired = T)

mydata3 <-distance %>% filter(Tissues=='YP1' | Tissues=='GR1') %>% filter (Distance=='Ratio')
mydata3
wilcox.test(Frequency~Tissues, data = mydata3)

mydata4 <-distance %>% filter(Tissues=='YP1' | Tissues=='TZ1') %>% filter (Distance=='Ratio')
mydata4
wilcox.test(Frequency~Tissues, data = mydata4)

#### For juicer HiCMap

distance <- read.table("Juicer.all.bed.bed", header=T)

distance
distance$Tissues <- as.factor(distance$Tissues)
mydata <-distance %>% filter(Tissues=='HL1' | Tissues=='GR1' | Tissues=='TZ1' | Tissues=='YP1') %>% filter (Distance=='Ratio')
ggplot(mydata, aes(x=Tissues, y=Frequency, fill=factor(Distance))) + geom_boxplot(outlier.shape = NA) +labs(fill = "Distance")

mydata <-distance %>% filter(Tissues=='HL1' | Tissues=='GR1') %>% filter (Distance=='Ratio')
mydata
wilcox.test(Frequency~Tissues, data = mydata, paired = T)

#### For juicer plot

distance1 <- read.table("GR_HL_TZ_YP.bed.bed.out", header=T)
ggplot(distance1 %>%
         group_by(Tissue) %>%
         mutate(weight = 1 / n()),
       aes(x = Distance, fill = Tissue, colour=Tissue)) +
  geom_line(aes(weight = Count), stat = 'density', position='dodge')

