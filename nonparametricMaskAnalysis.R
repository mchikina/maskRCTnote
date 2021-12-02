library(readr)
library(dplyr)
library(ggpubr)


##################################
## household selection analysis ## 
##################################

hhdata.raw=read_delim("baseHH_data.csv.gz",delim = ',', col_names = T)


#create surgical, cloth group variable
hhdata=hhdata.raw
hhdata$treatGroup=character(nrow(hhdata))
hhdata$treatGroup[]="surgical"
hhdata$treatGroup[hhdata$cloth==1|hhdata$cloth_pair==1]="cloth"
table(hhdata$treatment, hhdata$treatGroup)

#create a unified "both" group
hhdataUnified=hhdata
hhdataUnified$treatGroup="both"
#append the two together 
#this duplicates each entry but simplifies scripting
hhdata=rbind(hhdata, hhdataUnified)
table(hhdata$treatGroup, hhdata$treatment)


#define variables
myvarsHH=c("consentFrac", "consentNotNA",  "vSize") #internal
mynamesHH=c("fraction consented", "fraction reached", "mapped households") #full names for plotting

names(mynamesHH)=myvarsHH

treatTypes=c("surgical", "cloth", "both")

resMat=matrix(ncol=2, nrow=length(myvarsHH)*3)
rownames(resMat)=paste(rep(treatTypes, each=length(myvarsHH)), rep(myvarsHH,3))
rownames(resMat)
colnames(resMat)=c("p.val", "effect")
plhh=list() #this will hold the plots


for( tType in c("surgical", "cloth", "both")){
  
  #treatment and pair ID are defined per village 
  #using mean instead of first for treatment and pair_id so that we can catch any errors downstream
  groupedDataHH=hhdata %>% filter(treatGroup==tType) %>% #filter (above_60==1) %>% #uncomment to filter for above_60
    group_by(union) %>% #union is village
    summarize(consentNotNA=mean(!is.na(consent)), #fraction reached
              consentFrac = mean(consent=="Yes", na.rm=T),  #consented/reached
              #consent_sum=sum(consent=="Yes", na.rm=T),
              treat=mean(treatment), 
              vSize=n(), #village size in HH units, number of entries per union
              pair_id=mean(pairID))
#  groupedDataHH$consent_rateNA=groupedDataHH$consent_sum/groupedDataHH$vSize
  stopifnot(all(table(groupedDataHH$pair_id)==2)) #check we hav exactly two of each pair_id
  
  uIDs=unique(groupedDataHH$pair_id)
  
  iicontrol=which(groupedDataHH$treat==0)
  iicontrol=iicontrol[match(uIDs, groupedDataHH$pair_id[iicontrol])]
  
  iitreat=which(groupedDataHH$treat==1)
  iitreat=iitreat[match(uIDs, groupedDataHH$pair_id[iitreat])]
  
  
#make separate control ant treatment tables
  groupedDataC=groupedDataHH[iicontrol,]
  groupedDataT=groupedDataHH[iitreat,]

    #check that they are matched
  stopifnot(all(groupedDataC$pair_id==groupedDataT$pair_id))
  
  #put them back together in matched order for the plot
  plotData=rbind(groupedDataC, groupedDataT)
  plotData$treat=c("control", "treatment")[plotData$treat+1]
  
  for( i in 1:length(myvarsHH)){
    
    var=myvarsHH[i]
    var.name=mynamesHH[i]
    test.name=paste(tType,var)
    #record the p-value
    resMat[test.name, "p.val"]=(wilcox.test(groupedDataC[,var, drop=T], 
                                           groupedDataT[,var, drop=T], paired = T))$p.value
    #record the effect
    meanC=mean(groupedDataC[, var, drop=T])
    meanT=mean(groupedDataT[, var, drop=T])
    resMat[test.name, "effect"]=(meanT-meanC)/(meanC)
    
    plhh[[test.name]]=ggboxplot(plotData, x="treat", y=var, outlier.shape=NA, color="black", fill="treat")+ylim(c(min(plotData[, var]), max(plotData[, var])*1.03))+
      geom_dotplot(aes_string(x="treat", y=var),alpha=1, binaxis = "y", stackdir = "centerwhole", binwidth = diff(range(plotData[, var]))/80)+
      stat_compare_means(method="wilcox.test", paired=T)+
      #ggtitle(paste(tType,": ", var.name, sep=""))+
      ggtitle(var.name)+
      xlab("")+ylab(var.name)
  }
}


resMatHH=resMat

###############################
## individual level analysis ## 
###############################

bdata.raw=read_delim("endlineBlood_data.csv.gz",delim = ',', col_names = T)


#assign groups 
bdata=bdata.raw
bdata$treatGroup=character(nrow(bdata))
bdata$treatGroup[]="surgical"
bdata$treatGroup[bdata$treat_cloth==1|bdata$cloth_pair==1]="cloth"

#make unified dataset
bdataUnified=bdata
bdataUnified$treatGroup="both"
bdata=rbind(bdata, bdataUnified)




myvars=c("meanPos", "sumPos", "meanSymp", "sumSymp", "SDprop","vSize", "propMask")
mynames=c("symptomatic seropositive rate", "symptomatic seropositive count",
          "symptomatic rate", "symptomatic count", "social distancing proportion",
          "consenting population size", "proper mask fraction")

names(mynames)=myvars


resMat=matrix(ncol=2, nrow=length(myvars)*3)
rownames(resMat)=paste(rep(treatTypes, each=length(myvars)), rep(myvars,3))
rownames(resMat)
colnames(resMat)=c("p.val", "effect")

pl=list() #plot List


for( tType in c("surgical", "cloth", "both")){
  
  #treatment and pair ID are defined per village 
  #using mean instead of first for treatment and pair_id so that we can catch any errors downstream
  groupedData=bdata %>% filter(treatGroup==tType) %>%  
    group_by(union) %>%
    summarize(meanPos = mean(posXsymp), sumPos=sum(posXsymp),
              meanSymp=mean(symp), sumSymp=sum(symp), 
              SDprop=mean(soc_dist_prop), 
              propMask=mean(proper_mask_prop),
              treat=mean(treatment), vSize=n(), pair_id=mean(pairID))

 
  stopifnot(all(table(groupedData$pair_id)==2)) #check we hav exactly two of each pair_id
  
  uIDs=unique(groupedData$pair_id)
  
  iicontrol=which(groupedData$treat==0)
  iicontrol=iicontrol[match(uIDs, groupedData$pair_id[iicontrol])]
  
  iitreat=which(groupedData$treat==1)
  iitreat=iitreat[match(uIDs, groupedData$pair_id[iitreat])]
  
  
  #make separate control and treatment tables
  groupedDataC=groupedData[iicontrol,]
  groupedDataT=groupedData[iitreat,]
  
  #check that they are matched
  stopifnot(all(groupedDataC$pair_id==groupedDataT$pair_id))
  
  #put them back together in matched order for the plot
  plotData=rbind(groupedDataC, groupedDataT)
  plotData$treat=c("control", "treatment")[plotData$treat+1]
  
  for( i in 1:length(myvars)){
    
    var=myvars[i]
    var.name=mynames[i]
    test.name=paste(tType,var)
    #record the p-value
    resMat[test.name, "p.val"]=(wilcox.test(groupedDataC[,var, drop=T], 
                                            groupedDataT[,var, drop=T], paired = T))$p.value
    #record the effect
    meanC=mean(groupedDataC[, var, drop=T])
    meanT=mean(groupedDataT[, var, drop=T])
    resMat[test.name, "effect"]=(meanT-meanC)/(meanC)
    
    pl[[test.name]]=ggboxplot(plotData, x="treat", y=var, outlier.shape=NA, color="black", fill="treat")+ylim(c(min(plotData[, var]), max(plotData[, var])*1.03))+
      geom_dotplot(aes_string(x="treat", y=var),alpha=1, binaxis = "y", stackdir = "centerwhole", binwidth = diff(range(plotData[, var]))/80)+
      stat_compare_means(method="wilcox.test", paired=T)+
      #ggtitle(paste(tType,": ", var.name, sep=""))+
      ggtitle(var.name)+
      xlab("")+ylab(var.name)
  }
}

resMatInd=resMat



#resMatHH and resMatInd have the stats
#pl and plHH are the plots
names(pl)
pl[[1]]

ii=grep("surgical", names(pl))
ggarrange(plotlist = pl[ii[c(2,4,6,1,3,5)]], common.legend = T, nrow=2, ncol=3)


ii=grep("cloth", names(pl))
ggarrange(plotlist = pl[ii[c(2,4,6,1,3,5)]], common.legend = T, nrow=2, ncol=3)

ii=grep("both", names(pl))
ggarrange(plotlist = pl[ii[c(2,4,6,1,3,5)]], common.legend = T, nrow=2, ncol=3)


ggarrange(plotlist = plhh[c(3,2,1)+6], common.legend = T, nrow=1, ncol=3)

