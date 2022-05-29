library("ggplot2")
library("dplyr")

##### Arrowhead for Boundaries TAU
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_bound/Boundaries.three.tau.out.out",head=T)
head(mydata3)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-0.1,1.1)+labs(fill = "Types") +
  labs(title="TAU for boundaries based Arrowhead output", x= "Types", y= "Expression specifity index (Tau)")

mydata3_filter <-mydata3 %>% filter(Types=='conserve3' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Types=='conserve2' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Types=='conserve1' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

#mydata3_filter <-mydata3 %>% filter(Tissues=='GR')  %>% filter(Types=='conserved' | Types=='specific')
#wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 
#mydata3_filter <-mydata3 %>% filter(Tissues=='HL')  %>% filter(Types=='conserved' | Types=='specific')
#wilcox.test(Density~Types, data = mydata3_filter,alternative = "less", paired = FALSE) 
#mydata3_filter <-mydata3 %>% filter(Tissues=='TZ')  %>% filter(Types=='conserved' | Types=='specific')
#wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 
#mydata3_filter <-mydata3 %>% filter(Tissues=='YP')  %>% filter(Types=='conserved' | Types=='specific')
#wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 
##head(mydata3_filter)
##t.test(Density~Types, data = mydata3_filter, alternative = "less", var.equal = FALSE)


##### Arrowhead for Boundaries FC
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_bound/Boundaries.fold.out.out",head=T)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-1,15)+labs(fill = "Types") +
    labs(title="absolute log2FC for boundaries based Arrowhead output", x= "Types", y= "abs(log2FC)")

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,alternative = "less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='TZvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsHL')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 


################ TopDom for Boundaries FC
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_bound/Boundaries.fold.out.out",head=T)

mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_bound/40k/Boundaries.fold.out.out",head=T)
mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-1,7.5)+labs(fill = "Types") +
  labs(title="absolute log2FC for boundaries based TopDom output", x= "Types", y= "abs(log2FC)")

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,alternative = "less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='TZvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsHL')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less",paired = FALSE) 


################ TopDom for Boundaries TAU
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_bound/Boundaries.three.tau.out.out",head=T)

mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-0.1,1.1)+labs(fill = "Types") +
  labs(title="TAU for boundaries based TopDom output", x= "Types", y= "Expression specifity index (Tau)")

mydata3_filter <-mydata3 %>% filter(Types=='conserve3' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Types=='conserve2' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less", paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Types=='conserve1' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, alternative = "less", paired = FALSE) 


#####################################################
#####################################################
#####################################################

###############  TopDom for TAD  FC
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_TAD/Boundaries.fold.out1.active.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_TAD/Boundaries.fold.out1.inactive.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_TAD/Boundaries.fold.out.out",head=T)

mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-1,15)+labs(fill = "Types") +
  labs(title="absolute log2FC for all TADs based TopDom output", x= "Types", y= "abs(log2FC)")

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='TZvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsHL')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

###############  TopDom for TAD  TAU
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_TAD/Boundaries.three.tau.out1.active.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_TAD/Boundaries.three.tau.out1.inactive.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/TopDom_TAD/Boundaries.three.tau.out.out",head=T)

mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-0.1,1.1)+labs(fill = "Types") +
  labs(title="TAU for all TADs based TopDom output", x= "Types", y= "Expression specifity index (Tau)")

mydata3_filter <-mydata3 %>% filter(Types=='conserve3' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Types=='conserve2' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE)

mydata3_filter <-mydata3 %>% filter(Types=='conserve1' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)

############### Arrowhead for TAD  FC
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_TAD/Boundaries.fold.out1.inactive.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_TAD/Boundaries.fold.out1.active.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_TAD/Boundaries.fold.out.out",head=T)

mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-1,15)+labs(fill = "Types") +
  labs(title="absolute log2FC for active TADs based Arrowhead output", x= "Types", y= "abs(log2FC)")

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='HLvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='TZvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsGR')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsHL')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Tissues=='YPvsTZ')  %>% filter(Types=='shared' | Types=='speci')
wilcox.test(Density~Types, data = mydata3_filter,paired = FALSE) 

###############  Arrowhead for TAD  TAU
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_TAD/Boundaries.three.tau.out1.active.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_TAD/Boundaries.three.tau.out1.inactive.bed.bed",head=T)
mydata3 <-read.table("/home/yiliao/Documents/Pepper_2021/Final/5_Figure5_TAD/TAD_3/TAU/Arrowhead_TAD/Boundaries.three.tau.out.out",head=T)

mydata3$Tissues <- as.factor(mydata3$Tissues )
ggplot(mydata3, aes(x=Tissues, y=Density, fill=factor(Types))) + geom_boxplot(outlier.shape = NA) + ylim(-0.1,1.1)+labs(fill = "Types") + 
  labs(title="TAU for inactive TADs based TopDom output", x= "Types", y= "Expression specifity index (Tau)")

mydata3_filter <-mydata3 %>% filter(Types=='conserve3' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE) 

mydata3_filter <-mydata3 %>% filter(Types=='conserve2' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)

mydata3_filter <-mydata3 %>% filter(Types=='conserve1' | Types=='specific')
wilcox.test(Density~Types, data = mydata3_filter, paired = FALSE)
