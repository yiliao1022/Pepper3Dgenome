library("ggplot2")
library("dplyr")

setwd("/home/yiliao/Documents/Pepper_2021/Final/6_Figure_loop/FINAL")



#####Loop size
mydata1 <-read.table("All.bed.bed.bed",head=T)
head(mydata1)
mydata1$Resolution<- as.factor(mydata1$Resolution)
mydata1$Tissues<- as.factor(mydata1$Tissues)
ggplot(mydata1, aes(x=Resolution, y=Size, fill=factor(Tissues))) + geom_boxplot(outlier.shape = NA) + stat_summary(fun=mean, geom="point", shape=20, size=4, color="red", fill="red") + ylim(0,8000000) +labs(fill = "Tissues")

#######################

### For Loop tau
setwd("/home/yiliao/Documents/Pepper_2021/Final/6_Figure_loop/FINAL/TADvsLoop")
mydata2 <-read.table("Figure6.Final.tau.loop.tissues.bed.bed",head=T)
head(mydata2)
mydata2$Tissues<- as.factor(mydata2$Tissues)
ggplot(mydata2, aes(x=Tissues, y=Tau, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(0,1) +labs(fill = "Types")

mydata2_filter <-mydata2 %>% filter(Tissues=='All') %>% filter(Types=='Specific' | Types=='Shared')
wilcox.test(Tau~Types, data = mydata2_filter,  alternative ="less", paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='GR') %>% filter(Types=='Specific' | Types=='Shared')
wilcox.test(Tau~Types, data = mydata2_filter,  alternative ="less", paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='HL') %>% filter(Types=='Specific' | Types=='Shared')
wilcox.test(Tau~Types, data = mydata2_filter,  alternative ="less", paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='TZ') %>% filter(Types=='Specific' | Types=='Shared')
wilcox.test(Tau~Types, data = mydata2_filter, alternative ="less", paired = FALSE)

mydata2_filter <-mydata2 %>% filter(Tissues=='YP') %>% filter(Types=='Specific' | Types=='Shared')
wilcox.test(Tau~Types, data = mydata2_filter, alternative ="less", paired = FALSE)


####################### For Loop Fold change 

mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/6_Figure_loop/20210902/Boundaries.loop.fold.out.out",head=T)

mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-1,15)+labs(fill = "Types") +
  labs(title="absolute log2FC for loops based HiCexplorer output", x= "Types", y= "abs(log2FC)")

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative ="less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative ="less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='TZvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative ="less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,  alternative ="less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsHL')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative ="less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,  alternative ="less", paired = FALSE) 






