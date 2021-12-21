library(readr)
library(dplyr)
library(ggpubr)

library('inline')
library('Rcpp')

# R version of the C++ function
pair_permute=function(A,B, nperm){
  np=length(A)
  effectR=double(np)
  for (iter in 1:nperm){
    totalA=0
    totalB=0
    for(i in 1:np){
      ru=runif(1)
      if(ru>0.5){
        totalA=totalA+A[i]
        totalB=totalB+B[i]
      }
      else{
        totalA=totalA+B[i]
        totalB=totalB+A[i]
      }
     
    }
    effectR[iter]=(totalB-totalA)/totalA
    
  }
  effectR
}


rcppPlugin <- getPlugin( "Rcpp" )

rcppPlugin$env$PKG_CXXFLAGS="-O3"

sourceCpp("PairPermute.cpp")


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




myvars=c("meanPos", "sumPos", "meanSymp", "sumSymp", "SDprop","vSize", "propMask", "posOverSymp")
mynames=c("symptomatic seropositive rate", "symptomatic seropositive count",
          "symptomatic rate", "symptomatic count", "social distancing proportion",
          "consenting population size", "proper mask fraction", "postive conditioned on symptoms")


#myvars=c("meanSymp", "sumSymp", "vSize")
#mynames=c("symptomatic rate", "symptomatic count", 
#          "consenting population size")


names(mynames)=myvars



matList=list()
permRes=list()
nperm=10000
npermCpp=1e6
for( tType in c("surgical", "cloth", "both")){
  resMat=matrix(ncol=nperm+1, nrow=length(myvars))

  rownames(resMat)=myvars
  colnames(resMat)=c("effect", paste0("effectR", 1:nperm))
  #treatment and pair ID are defined per village 
  #using mean instead of first for treatment and pair_id so that we can catch any errors downstream
  groupedData=bdata %>% filter(treatGroup==tType) %>% # filter(above_60==1)%>%
    group_by(union) %>%
    summarize(meanPos = mean(posXsymp), sumPos=sum(posXsymp),
              meanSymp=mean(symp), sumSymp=sum(symp), 
              SDprop=mean(soc_dist_prop),   
              propMask=mean(proper_mask_prop),
              treat=mean(treatment), vSize=n(), pair_id=mean(pairID))

 groupedData$posOverSymp=(groupedData$sumPos)/(groupedData$sumSymp)
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
  

  
  for( i in 1:length(myvars)){

    var=myvars[i]
    var.name=mynames[i]
    test.name=paste(tType,var)
    
    #record the effect
    meanC=mean(groupedDataC[, var, drop=T])
    meanT=mean(groupedDataT[, var, drop=T])
    resMat[var, "effect"]=(meanT-meanC)/(meanC)
    

    
  
    resMat[var,-1]=pair_permute(groupedDataC[, var, drop=T], groupedDataT[, var, drop=T], nperm)
#        tmp=pair_permute(groupedDataC[, var, drop=T], groupedDataT[, var, drop=T], nperm)
    matList[[tType]]=resMat
    #do some real permutations
    x=resMat[var, -1]
    #pp=pretty(x*1.5, n=45)
   pp=seq(-0.25, 0.25, length.out = 51)
   pp  
   ppmids=pp[-length(pp)]+(step<-(pp[2]-pp[1])/2)
    tmp=permute(groupedDataC[, var,T] , groupedDataT[,var,T], pp[1],pp[length(pp)],length(pp)-1, npermCpp)
    permRes[[tType]][[var]]$pval=tmp$second
    permRes[[tType]][[var]]$counts=tmp$first
    permRes[[tType]][[var]]$mids=ppmids
        permRes[[tType]][[var]]$step=step
  
    if(F){
      #permuted values from R
      df=data.frame(x=matList[[tType]][var, -1])

      #binned values from cpp
      df2=data.frame(x=permRes[[tType]][[var]]$mids, y=permRes[[tType]][[var]]$counts)
  
      df2=df2[df2$y>0,]
      
      #histogram data
      tmp=cbind((table(cut(df$x, df2$x-step))))
    
      #subsample the R daata
      iis=sample(1:length(df$x), length(df$x)/2)  
      tmp1=cbind((table(cut(df$x[iis], df2$x-step))))
      tmp2=cbind((table(cut(df$x[-iis], df2$x-step))))
      
     if(F){ #check that we are getting same results from R and cpp up to random error
      par(mfrow=c(1,2))
      plot(tmp, df2[-nrow(df2),"y"], col=((1:length(tmp))>length(tmp)/2)+1, xlab="R", ylab="CPP", main="Histogram counts");abline(a=0, b=1, lwd=2)
      plot(tmp1, tmp2, col=((1:length(tmp))>length(tmp)/2)+1,xlab="sample1", ylab="sample2", main="Histogram Counts");abline(a=0, b=1, lwd=2)
       #print(ggarrange(p1+ggtitle(var),p2, ncol=1))
     }
      }
    
    }
  
}
if(F){#R data
pl=list()
tType="surgical"
for (var in myvars){
  df=data.frame(x=matList[[tType]][var, -1])
  #compute p-value
  thisE=matList[[tType]][var, 1]
    randomE=matList[[tType]][var, -1]
  if(thisE<0){
    p=(sum(randomE<thisE)+1)/(sum(randomE<0)+1)
  }
    else{
      p=(sum(randomE>thisE)+1)/(sum(randomE>0)+1)
    }
  
      pl[[var]]=          gghistogram(df, x="x", fill="grey")+geom_vline(xintercept = matList[[tType]][var,1], 
                                                        color="red", size=3)+ggtitle(paste(mynames[var], "p=", format.pval(p)))+xlab("")
   
      }

}

pl=list()
tType="both"
for (var in myvars){
  df=data.frame(x=matList[[tType]][var, -1])
  #compute p-value
  thisE=(matList[[tType]][var, 1])*100
  show(thisE)
  p=permRes[[tType]][[var]]$pval
  mids=permRes[[tType]][[var]]$mids*100
    counts=permRes[[tType]][[var]]$counts
  dfplot=data.frame(effect=mids, counts=counts)
  dfplot[dfplot$counts>0,]

  x=matList[[tType]][var, -1]*100
  qq=quantile(x, c(0.05, 0.95))
  qq=round(qq,1)
  show(qq)
dfqq=data.frame(q1=qq[1], q2=qq[2])
    pl[[var]]= ggplot(dfplot)+
    geom_rect(data = dfqq,aes(xmin=q1, xmax=q2, ymin=0, ymax=Inf), fill="#0000FF52")+
    geom_bar(aes(x=effect, y=counts), 
             stat="identity", fill="lightgrey", color="black")+
  theme_bw(base_size = 12) +
    geom_vline(xintercept =thisE,  color="red", size=2)+
    ggtitle(paste(mynames[var], "p=", format.pval(p, eps = 1/npermCpp)))+
    xlab("")+xlim(c(-0.2, 0.2)*100)
  
}

#plot all the variables
#ggarrange(plotlist = pl, nrow=4, ncol=2)
#plot manuscript variables
ggarrange(pl[[1]],pl[[3]], pl[[6]]+xlab("effect size (%)"), nrow=3, ncol=1)
#change to serif
