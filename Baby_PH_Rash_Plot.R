##########################
####Truefc staill have problem
####################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
#filename <- args[1]
filename="Buttocks.txt"
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);

A$Rash_Buttocks_Num=as.numeric(A$Rash_Buttocks)/10
A$PH_Buttocks_Num= as.numeric(A$PH_Buttocks)/100
A$TEWL_Buttocks_Num= as.numeric(A$TEWL_Buttocks)/100
A$Rash_Perianal_Num=as.numeric(A$Rash_Perianal)/10
A$Rash_Perianal_Category = A$Rash_Perianal_Num
A$Rash_Perianal_Category[A$Rash_Perianal_Num>=0.05]="MildRash(0.5,1)"
A$Rash_Perianal_Category[A$Rash_Perianal_Num>1]="HighRash(>=1.5)"
A$Rash_Perianal_Category[A$Rash_Perianal_Num==0]="NoRash"

A$Rash_Buttocks_Category = A$Rash_Buttocks_Num
A$Rash_Buttocks_Category[A$Rash_Buttocks_Num>=0.05]="MildRash(0.5,1)"
A$Rash_Buttocks_Category[A$Rash_Buttocks_Num>1]="HighRash(>=1.5)"
A$Rash_Buttocks_Category[A$Rash_Buttocks_Num==0]="NoRash"

A$PH_Buttocks_Category = A$PH_Buttocks_Num
A$PH_Buttocks_Category[A$PH_Buttocks_Num>5.5]="HighPH >5.5"
A$PH_Buttocks_Category[A$PH_Buttocks_Num<=5.5]="LowPH <=5.5"

order_levels <- c("NoRash", "MildRash(0.5,1)", "HighRash(>=1.5)")
A$Rash_Buttocks_Category <- factor(A$Rash_Buttocks_Category, levels = order_levels, ordered = TRUE)
A$Rash_Perianal_Category <- factor(A$Rash_Perianal_Category, levels = order_levels, ordered = TRUE)

A$Rash2_Perianal_Category = A$Rash_Perianal_Num
A$Rash2_Perianal_Category[A$Rash_Perianal_Num==0.05]="NotClear"
A$Rash2_Perianal_Category[A$Rash_Perianal_Num>1]="Rash(>=1)"
A$Rash2_Perianal_Category[A$Rash_Perianal_Num==0]="NoRash"

A$Rash2_Buttocks_Category = A$Rash_Buttocks_Num
A$Rash2_Buttocks_Category[A$Rash_Buttocks_Num==0.05]="NotClear"
A$Rash2_Buttocks_Category[A$Rash_Buttocks_Num>=1]="Rash(>=1)"
A$Rash2_Buttocks_Category[A$Rash_Buttocks_Num==0]="NoRash"

B <- A %>%
  filter(!is.na(Rash_Perianal_Category))

C <- A %>%
  filter(Country =="GERMANY")

D <- B %>%
  filter(Country =="GERMANY")

E <- A %>%
  filter(Country =="US")

F <- B %>%
  filter(Country =="US")

G <- A %>%
  filter(Country =="CHINA")

H <- B %>%
  filter(Country =="CHINA")


png(filename=paste0("A.Buttock_TEWL.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
stat.test <- A %>%
  t_test(TEWL_Buttocks_Num ~ Rash_Buttocks_Category) %>%
  mutate(y.position = 5)
bxp <- ggboxplot(A, x = "Rash_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "Rash_Buttocks_Category")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none") + scale_x_discrete(limits = order_levels)
dev.off()


png(filename=paste0("A.Buttock_TEWL.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
stat.test <- A %>%
  t_test(TEWL_Buttocks_Num ~ PH_Buttocks_Category) %>%
  mutate(y.position = 5)
bxp <- ggboxplot(A, x = "PH_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "PH_Buttocks_Category")
bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
dev.off()


png(filename=paste0("A.Buttock_PH.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- A %>%
   t_test(PH_Buttocks_Num ~ Rash_Buttocks_Category) %>%
  mutate(y.position = 7.5)
 bxp <- ggboxplot(A, x = "Rash_Buttocks_Category", y = "PH_Buttocks_Num", color = "Rash_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()

 png(filename=paste0("A.Buttock_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- A %>%
   t_test(Rash_Buttocks_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 2.1)
 bxp <- ggboxplot(A, x = "PH_Buttocks_Category", y = "Rash_Buttocks_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()

 png(filename=paste0("A.Buttock_PH.Perianal_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- A %>%
   t_test(PH_Buttocks_Num ~ Rash_Perianal_Category) %>%
   mutate(y.position = 7.5)
 bxp <- ggboxplot(A, x = "Rash_Perianal_Category", y = "PH_Buttocks_Num", color = "Rash_Perianal_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()
 
 png(filename=paste0("A.Perianal_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- A %>%
   t_test(Rash_Perianal_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 2.1)
 bxp <- ggboxplot(A, x = "PH_Buttocks_Category", y = "Rash_Perianal_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()
 
 
 ####################################################################################
 
 png(filename=paste0("B.Buttock_TEWL.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- B %>%
   t_test(TEWL_Buttocks_Num ~ Rash_Buttocks_Category) %>%
   mutate(y.position = 5)
 bxp <- ggboxplot(B, x = "Rash_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "Rash_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none") + scale_x_discrete(limits = order_levels)
 dev.off()
 
 png(filename=paste0("B.Buttock_TEWL.Buttocks_Rash2cat.png"), width=1500, height=1600, res=300)
 stat.test <- B %>%
   t_test(TEWL_Buttocks_Num ~ Rash2_Buttocks_Category) %>%
   mutate(y.position = 5)
 bxp <- ggboxplot(B, x = "Rash2_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "Rash2_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none") + scale_x_discrete(limits = order_levels)
 dev.off()
 
 
 
 
 png(filename=paste0("B.Buttock_TEWL.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- B %>%
   t_test(TEWL_Buttocks_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 5)
 bxp <- ggboxplot(B, x = "PH_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()
 
 
 png(filename=paste0("B.Buttock_PH.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- B %>%
   t_test(PH_Buttocks_Num ~ Rash_Buttocks_Category) %>%
   mutate(y.position = 7.5)
 bxp <- ggboxplot(B, x = "Rash_Buttocks_Category", y = "PH_Buttocks_Num", color = "Rash_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()
 
 
 png(filename=paste0("B.Buttock_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- B %>%
   t_test(Rash_Buttocks_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 2.1)
 bxp <- ggboxplot(B, x = "PH_Buttocks_Category", y = "Rash_Buttocks_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()
 
 png(filename=paste0("B.Buttock_PH.Perianal_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- B %>%
   t_test(PH_Buttocks_Num ~ Rash_Perianal_Category) %>%
   mutate(y.position = 7.5)
 bxp <- ggboxplot(B, x = "Rash_Perianal_Category", y = "PH_Buttocks_Num", color = "Rash_Perianal_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()
 
 png(filename=paste0("B.Perianal_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- B %>%
   t_test(Rash_Perianal_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 2.1)
 bxp <- ggboxplot(B, x = "PH_Buttocks_Category", y = "Rash_Perianal_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()
 #########################################
 png(filename=paste0("C.Buttock_TEWL.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- C %>%
   t_test(TEWL_Buttocks_Num ~ Rash_Buttocks_Category) %>%
   mutate(y.position = 5)
 bxp <- ggboxplot(C, x = "Rash_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "Rash_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()
 
 
 png(filename=paste0("C.Buttock_TEWL.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- C %>%
   t_test(TEWL_Buttocks_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 5)
 bxp <- ggboxplot(C, x = "PH_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()
 

 png(filename=paste0("C.Buttock_PH.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- C %>%
   t_test(PH_Buttocks_Num ~ Rash_Buttocks_Category) %>%
   mutate(y.position = 7.5)
 bxp <- ggboxplot(C, x = "Rash_Buttocks_Category", y = "PH_Buttocks_Num", color = "Rash_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()
 
 png(filename=paste0("C.Buttock_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- C %>%
   t_test(Rash_Buttocks_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 2.1)
 bxp <- ggboxplot(C, x = "PH_Buttocks_Category", y = "Rash_Buttocks_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()
 
 png(filename=paste0("C.Buttock_PH.Perianal_Rashcat.png"), width=1500, height=1600, res=300)
 stat.test <- C %>%
   t_test(PH_Buttocks_Num ~ Rash_Perianal_Category) %>%
   mutate(y.position = 7.5)
 bxp <- ggboxplot(C, x = "Rash_Perianal_Category", y = "PH_Buttocks_Num", color = "Rash_Perianal_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
 dev.off()
 
 png(filename=paste0("C.Perianal_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
 stat.test <- C %>%
   t_test(Rash_Perianal_Num ~ PH_Buttocks_Category) %>%
   mutate(y.position = 2.1)
 bxp <- ggboxplot(C, x = "PH_Buttocks_Category", y = "Rash_Perianal_Num", color = "PH_Buttocks_Category")
 bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
   geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
 dev.off()

   ###################################
   
   png(filename=paste0("D.Buttock_TEWL.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
   stat.test <- D %>%
     t_test(TEWL_Buttocks_Num ~ Rash_Buttocks_Category) %>%
     mutate(y.position = 5)
   bxp <- ggboxplot(D, x = "Rash_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "Rash_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
   dev.off()
   
   png(filename=paste0("D.Buttock_TEWL.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
   stat.test <- D %>%
     t_test(TEWL_Buttocks_Num ~ PH_Buttocks_Category) %>%
     mutate(y.position = 5)
   bxp <- ggboxplot(D, x = "PH_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "PH_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
   dev.off()
   
   
   
   png(filename=paste0("D.Buttock_PH.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
   stat.test <- D %>%
     t_test(PH_Buttocks_Num ~ Rash_Buttocks_Category) %>%
     mutate(y.position = 7.5)
   bxp <- ggboxplot(D, x = "Rash_Buttocks_Category", y = "PH_Buttocks_Num", color = "Rash_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
   dev.off()
   
   
   png(filename=paste0("D.Buttock_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
   stat.test <- D %>%
     t_test(Rash_Buttocks_Num ~ PH_Buttocks_Category) %>%
     mutate(y.position = 2.1)
   bxp <- ggboxplot(D, x = "PH_Buttocks_Category", y = "Rash_Buttocks_Num", color = "PH_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
   dev.off() 
 
  png(filename=paste0("D.Buttock_PH.Perianal_Rashcat.png"), width=1500, height=1600, res=300)
   stat.test <- D %>%
       t_test(PH_Buttocks_Num ~ Rash_Perianal_Category) %>%
       mutate(y.position = 7.5)
   bxp <- ggboxplot(D, x = "Rash_Perianal_Category", y = "PH_Buttocks_Num", color = "Rash_Perianal_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
       geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
  dev.off()
 

     png(filename=paste0("D.Perianal_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
   stat.test <- D %>%
       t_test(Rash_Perianal_Num ~ PH_Buttocks_Category) %>%
       mutate(y.position = 2.1)
   bxp <- ggboxplot(D, x = "PH_Buttocks_Category", y = "Rash_Perianal_Num", color = "PH_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
       geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
   dev.off()
#######################################
   png(filename=paste0("E.Buttock_TEWL.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
   stat.test <- E %>%
     t_test(TEWL_Buttocks_Num ~ Rash_Buttocks_Category) %>%
     mutate(y.position = 5)
   bxp <- ggboxplot(E, x = "Rash_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "Rash_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
   dev.off()
   
   png(filename=paste0("E.Buttock_TEWL.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
   stat.test <- E %>%
     t_test(TEWL_Buttocks_Num ~ PH_Buttocks_Category) %>%
     mutate(y.position = 5)
   bxp <- ggboxplot(E, x = "PH_Buttocks_Category", y = "TEWL_Buttocks_Num", color = "PH_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
   dev.off()
   
   
   
   png(filename=paste0("E.Buttock_PH.Buttocks_Rashcat.png"), width=1500, height=1600, res=300)
   stat.test <- E %>%
     t_test(PH_Buttocks_Num ~ Rash_Buttocks_Category) %>%
     mutate(y.position = 7.5)
   bxp <- ggboxplot(E, x = "Rash_Buttocks_Category", y = "PH_Buttocks_Num", color = "Rash_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
   dev.off()
   
   
   png(filename=paste0("E.Buttock_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
   stat.test <- E %>%
     t_test(Rash_Buttocks_Num ~ PH_Buttocks_Category) %>%
     mutate(y.position = 2.1)
   bxp <- ggboxplot(E, x = "PH_Buttocks_Category", y = "Rash_Buttocks_Num", color = "PH_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
   dev.off() 
   
   png(filename=paste0("E.Buttock_PH.Perianal_Rashcat.png"), width=1500, height=1600, res=300)
   stat.test <- E %>%
     t_test(PH_Buttocks_Num ~ Rash_Perianal_Category) %>%
     mutate(y.position = 7.5)
   bxp <- ggboxplot(E, x = "Rash_Perianal_Category", y = "PH_Buttocks_Num", color = "Rash_Perianal_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")+ scale_x_discrete(limits = order_levels)
   dev.off()
   
   
   png(filename=paste0("E.Perianal_Rash.Buttocks_PHcat.png"), width=1000, height=1600, res=300)
   stat.test <- E %>%
     t_test(Rash_Perianal_Num ~ PH_Buttocks_Category) %>%
     mutate(y.position = 2.1)
   bxp <- ggboxplot(E, x = "PH_Buttocks_Category", y = "Rash_Perianal_Num", color = "PH_Buttocks_Category")
   bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
     geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none")
   dev.off()