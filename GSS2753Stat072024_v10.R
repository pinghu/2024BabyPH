##########################
####Truefc staill have problem
####################################
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
#filename="Buttock.7.RA"
#filename="Perianal.7.RA"
#filename="Buttock.short.7.RA"
#filename="Buttock_Perianal_Match_Samples_transposed.txt"
library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(tidyverse)

filename1="/Disk2/project/GSS2753LacticPH/2024July/Buttocks161_meta_data"
#filename1="Buttocks161_meta_data"
A1<-read.table(filename1, sep="\t", header=TRUE)
d <- dim(A1);

A1$Rash_Buttocks_Num=as.numeric(A1$Rash_Buttocks)/10
A1$PH_Num= as.numeric(A1$PH_Buttocks)/100
A1$TEWL_Num= as.numeric(A1$TEWL_Buttocks)/10
A1$Rash_Perianal_Num=as.numeric(A1$Rash_Perianal)/10

A1$PH_Category = A1$PH_Num
A1$PH_Category[A1$PH_Num>5.5]="HighPH >5.5"
A1$PH_Category[A1$PH_Num<=5.5]="LowPH <=5.5"
if (grepl("Buttock", filename)) {
  A1$Rash_Num <- A1$Rash_Buttocks_Num
}else if(grepl("Perianal", filename)){
  A1$Rash_Num <- A1$Rash_Perianal_Num
}

A1$Rash_Category = A1$Rash_Num
A1$Rash_Category[A1$Rash_Num>=0.05]="MildRash(0.5-1)"
A1$Rash_Category[A1$Rash_Num>1]="HighRash(>=1.5)"
A1$Rash_Category[A1$Rash_Num==0]="NoRash"
A1$RashPH_Category <-paste(A1$Rash_Category,A1$PH_Category )

A1$Rash2_Category = A1$Rash_Num
A1$Rash2_Category[A1$Rash_Num>=0.5]="Rash(>=0.5)"
A1$Rash2_Category[A1$Rash_Num==0]="NoRash(0)"
A1$Rash2PH_Category <-paste(A1$Rash2_Category,A1$PH_Category )

order_levels2 <- c("LowPH <=5.5", "HighPH >5.5")
A1$PH_Category <- factor(A1$PH_Category, levels = order_levels2, ordered = TRUE)

order_levels <- c("NoRash", "MildRash(0.5-1)", "HighRash(>=1.5)")
A1$Rash_Category <- factor(A1$Rash_Category, levels = order_levels, ordered = TRUE)

A1$RashPH_Category <-paste(A1$Rash_Category,A1$PH_Category )
order_levels3 <- c("NoRash LowPH <=5.5",
                   "NoRash HighPH >5.5",
                   "MildRash(0.5-1) LowPH <=5.5",
                   "MildRash(0.5-1) HighPH >5.5",
                   "HighRash(>=1.5) LowPH <=5.5",
                   "HighRash(>=1.5) HighPH >5.5")
A1$RashPH_Category <- factor(A1$RashPH_Category, levels = order_levels3, ordered = TRUE)                   

my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}

my.wilcox.p.value <- function(...) {
    obj<-try(wilcox.test(...), silent=TRUE)
     if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.kruskal.p.value <- function(...) {
  obj<-try(kruskal.test(...), silent=TRUE)
  if (inherits(obj, "try-error")) return(NA) else return(obj$p.value)
}



truefc<-function(VVV){
  #print(VVV)
  if (is.finite(VVV )){
	  XXX=VVV
	  if(VVV==0){
	      XXX=NA
   	  }else if(VVV<1){
	      XXX=-1/VVV
    	}
	  return(XXX)
  }else{
    return("NA")
  }
}

A<-read.table(filename, sep="\t", header=TRUE)
d <- dim(A);
B=A[1:d[1], 2:d[2]]
ZZ=as.numeric(min(B[B>0&!is.na(B)]))/100
C=B+ZZ

Cname=colnames(A)[2:d[2]]
Clen=length(Cname) ##there are 3 annotation columns
Site=rep("NA", Clen)
Country=rep("NA", Clen)
Sex=rep("NA", Clen)
SID=rep("NA", Clen)
Age=rep("NA",Clen)
ID=rep("NA", Clen)
OID=rep("NA", Clen)
SiteCountry=rep("NA", Clen)
splitname<-strsplit(Cname, "[.]")


for(mm in  1:Clen ){
  Site[mm]=splitname[[mm]][1]
  Age[mm]=splitname[[mm]][5]
  Country[mm]=splitname[[mm]][2]
  Sex[mm]=splitname[[mm]][4]
  SID[mm]=splitname[[mm]][3]
  OID[mm]=splitname[[mm]][9]
}
SID=as.numeric(SID)

matching_dataset <- A1 %>%
  filter(SubjctID_Buttocks %in% SID) %>%
  arrange(match(SubjctID_Buttocks, SID))

gene_means1 <- matching_dataset %>%
  group_by(RashPH_Category) %>%
  summarize(mean_gene = mean(Age, na.rm = TRUE)) 
tx4=paste(gene_means1$RashPH_Category, collapse =",")

result=paste("genename", "taxonName", "percentZeroNA", "meanALL","mLowPH(<=5.5)", "mHighPH(>5.5)",
              "mNoRash","mMildRash",  "mHighRash", tx4,
             "KP_PH_Category", "Pwilcox_HighPH_LowPH", "PCor_PH_spearman", "PCor_PH_pearson",
             "KP_Rash_Category","Pwilcox_HighRash_NoRash", "Pwilcox_HighRash_MildRash", "Pwilcox_MildRash_NoRash","Pcor_Rash_spearman", "PCor_Rash_pearson",
             "KP_RashPH_Category","Pwilcox_HighRashHighPH_NoRashLowPH","TEWL_PCor_spearman", "TEWL_PCor_pearson",
             "TFC_HighPH_LowPH","Cor_PH_spearman", "Cor_PH_pearson",
	     "TFC_HighRash_NoRash","TFC_HighRash_MildRash","TFC_MildRash_NoRash", 
	     "Cor_Rash_spearman", "Cor_Rash_pearson", 
             "TFC_HighRashHighPH_NoRashLowPH",
             "TEWL_COR_spearman", "TEWL_COR_pearson", 
             sep=",") 


print(result)
for (i in 1:d[1]){
  
  genename=A[i,1]
  gene<-as.numeric(C[i,])
  loggene<-log(gene)
  splitG<-strsplit(as.character(genename), "[.]")
  LLL=length(splitG[[1]])
  genus=splitG[[1]][LLL]
  cleaned_genus <- gsub("^\\w__", "", genus)  # Removes the prefix
  cleaned_genus <- gsub("_", " ", cleaned_genus)  # Replaces underscores with spaces
  percentZeroNA=(sum(B[i,]==0)+sum(is.na(B[i,])))/(d[2]-1)
  mydata0=data.frame(gene,loggene,  Site, Country,Sex, Age,SID)
  mydata <- inner_join(mydata0, matching_dataset, by = c("SID" = "SubjctID_Buttocks"))
  
  KP_PH_Category=kruskal.test(gene ~ PH_Category, data = mydata)$p.value
  KP_Rash2_Category=kruskal.test(gene ~ Rash2_Category, data = mydata)$p.value
  KP_Rash2PH_Category=kruskal.test(gene ~ Rash2PH_Category, data = mydata)$p.value
  KP_Rash_Category=kruskal.test(gene ~ Rash_Category, data = mydata)$p.value
  KP_RashPH_Category=kruskal.test(gene ~ RashPH_Category, data = mydata)$p.value
  
  # print(paste0("gene KP_PH_Category=", KP_PH_Category,
  #              "; KP_Rash2_Category=", KP_Rash2_Category,
  #              "; KP_Rash2PH_Category=", KP_Rash2PH_Category))
  minP=min(KP_PH_Category,
           KP_Rash_Category, KP_RashPH_Category
           )
  
    LowPH<-mydata$gene[mydata$PH_Category=="LowPH <=5.5"]; mLowPH=mean(LowPH);
    HighPH <- mydata$gene[mydata$PH_Category == "HighPH >5.5" ]; mHighPH=mean(HighPH);
    Pwilcox_HighPH_LowPH=my.wilcox.p.value(HighPH, LowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighPH_LowPH=truefc(mHighPH/mLowPH)
    
  
    NoRash<-mydata$gene[mydata$Rash2_Category=="NoRash(0)"]; mNoRash=mean(NoRash);
    Rash2<-mydata$gene[mydata$Rash2_Category =="Rash(>=0.5)"]; mRash2=mean(Rash2);
    Pwilcox_Rash2_NoRash=my.wilcox.p.value(Rash2,NoRash, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_Rash2_NoRash=truefc(mRash2/mNoRash)
    
    NoRash<-mydata$gene[mydata$Rash_Category=="NoRash"]; mNoRash=mean(NoRash);
    MildRash<-mydata$gene[mydata$Rash_Category=="MildRash(0.5-1)"]; mMildRash=mean(MildRash);
    HighRash<-mydata$gene[mydata$Rash_Category=="HighRash(>=1.5)"]; mHighRash=mean(HighRash);
    Pwilcox_HighRash_NoRash=my.wilcox.p.value(HighRash,NoRash, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighRash_NoRash=truefc(mHighRash/mNoRash)
    Pwilcox_HighRash_MildRash=my.wilcox.p.value(HighRash,MildRash, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighRash_MildRash=truefc(mHighRash/mMildRash)
    Pwilcox_MildRash_NoRash=my.wilcox.p.value(MildRash,NoRash, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_MildRash_NoRash=truefc(mMildRash/mNoRash)
   
    
     gene_means1 <- mydata %>%
      group_by(RashPH_Category) %>%
      summarize(mean_gene = mean(gene, na.rm = TRUE)) 
    tx4=paste(gene_means1$mean_gene, collapse =",")
    HighRashHighPH<-mydata$gene[mydata$RashPH_Category=="HighRash(>=1.5) HighPH >5.5"]; mHighRashHighPH=mean(HighRashHighPH);
    NoRashLowPH<-mydata$gene[mydata$RashPH_Category=="NoRash LowPH <=5.5"]; mNoRashLowPH=mean(NoRashLowPH);
    Pwilcox_HighRashHighPH_NoRashLowPH=my.wilcox.p.value(HighRashHighPH, NoRashLowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
    TFC_HighRashHighPH_NoRashLowPH = truefc(mHighRashHighPH/mNoRashLowPH)
    
   
    
    PH_COR_spearman<-cor(mydata$gene, mydata$PH_Num,method="spearman")
    PH_PCor_spearman<-cor.test( mydata$gene, mydata$PH_Num,method="spearman", use="pairwise.complete.obs")$p.value
    PH_COR_pearson<-cor(mydata$gene, mydata$PH_Num)
    PH_PCor_pearson<-cor.test( mydata$gene, mydata$PH_Num, use="pairwise.complete.obs")$p.value

    Rash_COR_spearman<-cor(mydata$gene, mydata$Rash_Num,method="spearman")
    Rash_PCor_spearman<-cor.test( mydata$gene, mydata$Rash_Num,method="spearman", use="pairwise.complete.obs")$p.value
    Rash_COR_pearson<-cor(mydata$gene, mydata$Rash_Num)
    Rash_PCor_pearson<-cor.test( mydata$gene, mydata$Rash_Num, use="pairwise.complete.obs")$p.value    
    

    TEWL_COR_spearman<-cor(mydata$gene, mydata$TEWL_Num,method="spearman")
    TEWL_PCor_spearman<-cor.test( mydata$gene, mydata$TEWL_Num,method="spearman", use="pairwise.complete.obs")$p.value
    TEWL_COR_pearson<-cor(mydata$gene, mydata$TEWL_Num)
    TEWL_PCor_pearson<-cor.test( mydata$gene, mydata$TEWL_Num, use="pairwise.complete.obs")$p.value     
    
    
      RashHighPH<-mydata$gene[mydata$Rash2PH_Category=="Rash(>=0.5) HighPH >5.5"]; mRashHighPH=mean(RashHighPH);
      NoRashLowPH<-mydata$gene[mydata$Rash2PH_Category=="NoRash(0) LowPH <=5.5"]; mNoRashLowPH=mean(NoRashLowPH);
      RashLowPH<-mydata$gene[mydata$Rash2PH_Category=="Rash(>=0.5) LowPH <=5.5"]; mRashLowPH=mean(RashLowPH);
      NoRashHighPH<-mydata$gene[mydata$Rash2PH_Category=="NoRash(0) HighPH >5.5"]; mNoRashHighPH=mean(NoRashHighPH);
      Pwilcox_RashHighPH_NoRashLowPH=my.wilcox.p.value(RashHighPH, NoRashLowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
      TFC_RashHighPH_NoRashLowPH = truefc(mRashHighPH/mNoRashLowPH)
      
      Pwilcox_RashHighPH_NoRashHighPH=my.wilcox.p.value(RashHighPH, NoRashHighPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
      TFC_RashHighPH_NoRashHighPH = truefc(mRashHighPH/mNoRashHighPH)
      
      Pwilcox_RashLowPH_NoRashLowPH=my.wilcox.p.value(RashLowPH, NoRashLowPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
      TFC_RashLowPH_NoRashLowPH = truefc(mRashLowPH/mNoRashLowPH)
      
      Pwilcox_RashLowPH_NoRashHighPH=my.wilcox.p.value(RashLowPH, NoRashHighPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
      TFC_RashLowPH_NoRashHighPH = truefc(mRashLowPH/mNoRashHighPH)
      
      Pwilcox_RashLowPH_RashHighPH=my.wilcox.p.value(RashLowPH, RashHighPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
      TFC_RashLowPH_RashHighPH = truefc(mRashLowPH/mRashHighPH) 
      
      Pwilcox_NoRashLowPH_NoRashHighPH=my.wilcox.p.value(NoRashLowPH, NoRashHighPH, na.rm=TRUE, paired = F, alternative = "two.sided", na.action=na.omit)
      TFC_NoRashLowPH_NoRashHighPH = truefc(mNoRashLowPH/mNoRashHighPH) 
      
      
      
      
      result=paste(genename, cleaned_genus,percentZeroNA, 
                   mean(gene),mLowPH, mHighPH,
                   mNoRash,mMildRash,mHighRash,  tx4,
                   KP_PH_Category, Pwilcox_HighPH_LowPH, PH_PCor_spearman, PH_PCor_pearson,
                   KP_Rash_Category,Pwilcox_HighRash_NoRash, Pwilcox_HighRash_MildRash, Pwilcox_MildRash_NoRash, Rash_PCor_spearman, Rash_PCor_pearson,
                   KP_RashPH_Category,Pwilcox_HighRashHighPH_NoRashLowPH,TEWL_PCor_spearman, TEWL_PCor_pearson,
                   TFC_HighPH_LowPH,PH_COR_spearman, PH_COR_pearson,
                   TFC_HighRash_NoRash,TFC_HighRash_MildRash,TFC_MildRash_NoRash,Rash_COR_spearman, Rash_COR_pearson, 
                   TFC_HighRashHighPH_NoRashLowPH,
                   TEWL_COR_spearman, TEWL_COR_pearson,
                   sep=",") 
     
      print(result)
##########################################################
      if(minP<=0.1){
        mydata$PH_Category <- as.character(mydata$PH_Category)
        mydata$Rash_Category <- as.character(mydata$Rash_Category)
        mydata$RashPH_Category <- as.character(mydata$RashPH_Category)
        mydata$Rash2_Category <- as.character(mydata$Rash2_Category)
        mydata$Rash3_Category <- as.character(mydata$Rash2PH_Category)
	
        long_data <- pivot_longer(mydata, cols = c("PH_Category", "Rash_Category", "RashPH_Category"
						     ),
                                  names_to = "Category", values_to = "CategoryValue")
        
        
        
        
        
      
      # Assuming long_data is already preprocessed
      summary_data <- long_data %>%
        group_by(Category, CategoryValue) %>%
        summarise(Mean = mean(gene), SE = sd(gene)/sqrt(n()), Count = n(), .groups = 'drop')
      
      order_levels_new <- c("LowPH <=5.5", "HighPH >5.5", "NoRash", "MildRash(0.5-1)", "HighRash(>=1.5)",
                            "NoRash LowPH <=5.5", "NoRash HighPH >5.5",
                            "MildRash(0.5-1) LowPH <=5.5", "MildRash(0.5-1) HighPH >5.5", 
                            "HighRash(>=1.5) LowPH <=5.5", "HighRash(>=1.5) HighPH >5.5"
      )
      
      summary_data$CategoryValue <- factor(summary_data$CategoryValue, levels = order_levels_new, ordered = TRUE)
      
     # png(filename = paste0(genename, ".", filename, ".error_bar.png"), width = 2500, height = 1200, res = 300)
      png(filename = paste0(genus, ".", filename, ".error_bar.png"), width = 1250, height = 600, res = 80)
      # Define the dodge width for consistent bar width across facets
      dodge_width <- 0.9
      
      # Create the bar plot
      plot <- ggplot(summary_data, aes(y = CategoryValue, x = Mean, fill = Category)) +
        geom_col(position = position_dodge(dodge_width), width = 0.7) +
        geom_errorbar(aes(xmin = Mean - SE, xmax = Mean + SE), width = 0.2, position = position_dodge(dodge_width)) +
       # geom_text(aes(label = Count, x = -0.02), vjust = 0.5, hjust = 1, position = position_dodge(dodge_width)) +       
        labs(y = "", x = genename, 
             title = paste(genename, "\n", filename,
                           "\npPH=", sprintf("%.2f", KP_PH_Category), 
                           "pRash=", sprintf("%.2f", KP_Rash_Category),
                           "pRashPH=", sprintf("%.2f", KP_RashPH_Category), 
                           "\n",cleaned_genus
             )) +
        theme_minimal() +
        geom_hline(yintercept = c(2.5, 5.5), linetype = "dotted", color = "grey50") +
        geom_vline(xintercept = 0, color = "black") +
        theme(legend.position = "none",
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white", colour = "black"), 
              plot.title = element_text(hjust = 0.5)  # Center-align the title
        )
      
      # Print and save the plot
      print(plot)
      dev.off()
      
      }

      # if( KP_PH_Category <=0.05){
      # 	    stat.test <- mydata %>%
      #          wilcox_test(gene ~ PH_Category) %>%
      #       mutate(y.position = max(gene))
      #     if(min(stat.test$p) <=0.05){
      #       png(filename=paste0(filename, ".", genus, ".PH.png"), width=600, height=800, res=80)
      #       bxp <- ggboxplot(mydata, x = "PH_Category", y = "gene", color = "PH_Category")
      #       bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
      #         ggtitle(paste0(filename, "\n", genus, " PH", "\nKpPH=", sprintf("%.2f", KP_PH_Category)))+
      #         geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
      #       dev.off()
      #     }
      # } 
      # if( KP_RashPH_Category <=0.05){
      #        mydata$RashPH_Category <- as.factor(mydata$RashPH_Category)
      #        mydata$RashPH_Category <- as.character(mydata$RashPH_Category)
      # 	      stat.test <- mydata %>%
      # 	        wilcox_test(gene ~ RashPH_Category) %>%
      # 	        mutate(y.position = max(gene))
      # 	      if(min(stat.test$p) <=0.05){
      # 	        png(filename=paste0(filename, ".", genus, ".RashPH.png"), width=600, height=800, res=80)
      # 	        bxp <- ggboxplot(mydata, x = "RashPH_Category", y = "gene", color = "RashPH_Category")
      # 	        bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
      # 	          ggtitle(paste0(filename, "\n", genus, " RashPH", "\nKpRashPH=", sprintf("%.2f", KP_RashPH_Category)))+
      # 	          geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
      # 	        dev.off()
      # 	        
      # 	      }
      # }  	      
      # if( KP_Rash_Category <=0.05){
      # 	        stat.test <- mydata %>%
      # 	          wilcox_test(gene ~ Rash_Category) %>%
      # 	          mutate(y.position = max(gene))
      # 	        if(min(stat.test$p) <=0.05){
      # 	          png(filename=paste0(filename, ".", genus, ".Rash.png"), width=600, height=800, res=80)
      # 	          bxp <- ggboxplot(mydata, x = "Rash_Category", y = "gene", color = "Rash_Category")
      # 	          bxp + stat_pvalue_manual(stat.test,   label = "p", tip.length = 0.01, step.increase = 0.1)+
      # 	            ggtitle(paste0(filename, "\n", genus, " RashPH", "\nKpRash=", sprintf("%.2f", KP_Rash_Category)))+
      # 	            geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.7) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
      # 	          dev.off()
      # 	          
      # 	        }
      # }
}
