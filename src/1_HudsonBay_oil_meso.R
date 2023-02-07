## Title:      Hudson Bay Data
## Created by: Paul Julian (pauljulianphd@gmail.com)
## Created on: 01/18/2023

## BAD ## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
#Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

# Libraries
# devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)
library(openxlsx)

library(flextable)
library(magrittr)
#Paths
wd="C:/Julian_LaCie/_GitHub/CNRC_Hudson"
paths=paste0(wd,c("/Plots/HudsonBay/","/Export/","/Data/","/GIS"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]
GIS.path=paths[4]

# -------------------------------------------------------------------------
# chem data ---------------------------------------------------------------
sample.xwalk=data.frame(Sample.Name=c("T=0, Treatment 2a", "T=0, Treatment 2b", "T=0, Treatment 2c", 
                                      "Treatment 4a", "Treatment 4b", "Treatment 4c", 
                                      "Treatment 3a", "Treatment 3b", "Treatment 3c",
                                      "Treatment 2a", "Treatment 2b", "Treatment 2c", 
                                      "Treatment 1a", "Treatment 1b", "Treatment 1c"),
                        Treatment=c(rep("Control",3),rep("T4",3),rep("T3",3),rep("T2",3),rep("T1",3)),
                        rep=rep(c("a","b","c"),5))

## Recovery Values
rec.frac=read.xlsx(paste0(data.path,"HudsonBay_20230118/Results_Hudson Bay Summarized Chemistry Data from Alexa-NFtagged.xlsx"),sheet=3)
colnames(rec.frac)=c("Sample.Name","PAH.rec.frac","Alkane.rec.frac")

##Alkanes
alk.dat=read.xlsx(paste0(data.path,"HudsonBay_20230118/Results_Hudson Bay Summarized Chemistry Data from Alexa-NFtagged.xlsx"),sheet=2)
alk.dat.melt=melt(alk.dat,id.vars = 'Sample.Name')
alk.dat.melt[alk.dat.melt$value=="<DL","value"]=0;# zero for now till replaced with MDLs or other method
alk.dat.melt$value=as.numeric(alk.dat.melt$value)

alk.dat.melt=merge(alk.dat.melt,sample.xwalk,"Sample.Name")
alk.dat.melt=merge(alk.dat.melt,rec.frac[,c("Sample.Name","Alkane.rec.frac")],"Sample.Name")

alk.dat.melt$PreRec=with(alk.dat.melt,value/Alkane.rec.frac)

alk.var.list=unique(as.character(alk.dat.melt$variable))
## PAH
PAH.dat=read.xlsx(paste0(data.path,"HudsonBay_20230118/Results_Hudson Bay Summarized Chemistry Data from Alexa-NFtagged.xlsx"),sheet=1)

PAH.dat.melt=melt(PAH.dat,id.vars = 'Sample.Name')
PAH.dat.melt$cen=with(PAH.dat.melt,ifelse(value=="<DL",T,F))
PAH.dat.melt[PAH.dat.melt$value=="<DL","value"]=0;# zero for now till replaced with MDLs or other method
PAH.dat.melt$value=as.numeric(PAH.dat.melt$value)

PAH.dat.melt$value.original=as.numeric(PAH.dat.melt$value)
## censored data
library(NADA)
cen.vars=unique(subset(PAH.dat.melt,cen==T)$variable)

ddply(PAH.dat.melt,"variable",summarise,N.cen=sum(cen,na.rm=T))

PAH.dat.melt2=data.frame()
for(i in 1:length(cen.vars)){
tmp=subset(PAH.dat.melt,variable==cen.vars[i])
myros=ros(tmp$value, tmp$cen)
myros=as.data.frame(myros)
myros.max=max(subset(myros,censored==T)$modeled)

tmp$value=with(tmp,ifelse(cen==T,myros.max,value));# replaced censored data with max ros value
PAH.dat.melt2=rbind(PAH.dat.melt2,tmp)
}



####

PAH.var.list=unique(as.character(PAH.dat.melt$variable))

TPAH.var.list=PAH.var.list[grepl("methy",PAH.var.list)==F&PAH.var.list!="17.alpha.(H),21.beta.(H)-Hopane(30)"]
methPAH.var.list=PAH.var.list[grepl("methy",PAH.var.list)]


## Outliers identifed by Lab QA/QC and/or project staff
treatment2a.outlier.vars=c("dibenzothiophene","C1-dibenzothiophene","C-2.dibenzothiophene","C-3.dibenzothiophene","C-1.Fluoranthene")
subset(PAH.dat.melt,Sample.Name=="Treatment 2a"&variable%in%treatment2a.outlier.vars)
PAH.dat.melt$value=with(PAH.dat.melt,ifelse(Sample.Name=="Treatment 2a"&
                                              variable%in%treatment2a.outlier.vars,
                                            NA,value))

treatment1c.outlier.vars=c("C-1.Fluoranthene","C-2.Fluoranthene")
subset(PAH.dat.melt,Sample.Name=="Treatment 1c"&variable%in%treatment1c.outlier.vars)
PAH.dat.melt$value=with(PAH.dat.melt,ifelse(Sample.Name=="Treatment 1c"&
                                              variable%in%treatment1c.outlier.vars,
                                            NA,value))

subset(PAH.dat.melt,Sample.Name=="Treatment 1b"&variable%in%"C-3.dibenzothiophene")
PAH.dat.melt$value=with(PAH.dat.melt,ifelse(Sample.Name=="Treatment 1b"&
                                              variable%in%"C-3.dibenzothiophene",
                                            NA,value))

subset(PAH.dat.melt,Sample.Name=="Treatment 2b"&variable%in%"C-2.Fluoranthene")
PAH.dat.melt$value=with(PAH.dat.melt,ifelse(Sample.Name=="Treatment 2b"&
                                              variable%in%"C-2.Fluoranthene",
                                            NA,value))

PAH.dat.melt=merge(PAH.dat.melt,sample.xwalk,"Sample.Name")
PAH.dat.melt=merge(PAH.dat.melt,rec.frac[,c("Sample.Name","PAH.rec.frac")],"Sample.Name")

PAH.dat.melt$PreRec=with(PAH.dat.melt,value/PAH.rec.frac)

vars=c("Sample.Name", "variable", "value", "Treatment", "rep","PreRec")
dat.melt.all=rbind(alk.dat.melt[,vars],
                   PAH.dat.melt[,vars])


hop.dat=subset(dat.melt.all,variable=="17.alpha.(H),21.beta.(H)-Hopane(30)")
hop.dat=hop.dat[,c("Sample.Name","PreRec")]
colnames(hop.dat)=c("Sample.Name","Hopane")
dat.melt.all=merge(dat.melt.all,hop.dat)

dat.melt.all$Norm=with(dat.melt.all,PreRec/Hopane)


total.PAH=ddply(subset(dat.melt.all,variable%in%TPAH.var.list),
                c("Sample.Name","Treatment", "rep"),summarise,param="T.PAH",
                TValue=sum(Norm,na.rm=T))
meth.PAH=ddply(subset(dat.melt.all,variable%in%methPAH.var.list),
                c("Sample.Name","Treatment", "rep"),summarise,param="methly.PAH",
                TValue=sum(Norm,na.rm=T))          
alk.PAH=ddply(subset(dat.melt.all,variable%in%alk.var.list),
               c("Sample.Name","Treatment", "rep"),summarise,param="alkanes",
               TValue=sum(Norm,na.rm=T))          

total.vals=rbind(total.PAH,meth.PAH,alk.PAH)

sum.stats=ddply(total.vals,c('Treatment',"param"),summarise,
                mean.val=mean(TValue,na.rm=T),
                sd.val=sd(TValue,na.rm=T))
sum.stats$Treatment=factor(sum.stats$Treatment,levels=c("Control","T4","T3","T2","T1"))
sum.stats=sum.stats[order(sum.stats$Treatment),]

cols=viridisLite::viridis(5,alpha=0.4,option="E")
# cols=grey.colors(5)
# png(filename=paste0(plot.path,"PAH_treat.png"),width=5.5,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.5,1.5),oma=c(3.5,4,1.5,0.25),lwd=0.75);
layout(matrix(1:3,3,1,byrow = T))

x.labs=c("ANSo, T=0","Untreated","Fertilizer","T2 ANSo","T1 ANSo\nFertilizer")
ylim.val=c(0,90);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(subset(sum.stats,param=="T.PAH")$mean.val,col=cols,ylim=ylim.val,axes=F)
with(subset(sum.stats,param=="T.PAH"),errorbars(x,mean.val,sd.val,"black",length=0.05))
axis_fun(2,ymaj,ymin,ymaj,lwd=1)
axis_fun(1,x,x,NA,padj=1,line=-0.5,lwd=1);box(lwd=1)
mtext(side=2,line=2.5,"Total PAHs (\u03BCg g\u207B\u00B9)")

ylim.val=c(0,30);by.y=10;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(subset(sum.stats,param=="methly.PAH")$mean.val,col=cols,ylim=ylim.val,axes=F)
with(subset(sum.stats,param=="methly.PAH"),errorbars(x,mean.val,sd.val,"black",length=0.05))
axis_fun(2,ymaj,ymin,ymaj,lwd=1)
axis_fun(1,x,x,NA,padj=1,line=-0.5,lwd=1);box(lwd=1)
mtext(side=2,line=2.5,"Total Methylated\nPAHs (\u03BCg g\u207B\u00B9)")

ylim.val=c(0,140);by.y=40;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
x=barplot(subset(sum.stats,param=="alkanes")$mean.val,col=cols,ylim=ylim.val,axes=F)
with(subset(sum.stats,param=="alkanes"),errorbars(x,mean.val,sd.val,"black",length=0.05))
axis_fun(2,ymaj,ymin,ymaj,lwd=1)
axis_fun(1,x,x,x.labs,padj=1,line=-0.75,lwd=1);box(lwd=1)
mtext(side=2,line=2.5,"Total Alkanes (\u03BCg g\u207B\u00B9)")
mtext(side=1,line=3,"Treatments")
dev.off()


# genomics ----------------------------------------------------------------

dat.meta=read.table(
  paste0(data.path,"Genomics/ISCHudBay_MG.feature_table_normalized_bacteriaArchaea_L7.tsv"),
  fill=T,header=T,sep="\t")
names(dat.meta)

treat.xwalk=data.frame(meta.colnm=names(dat.meta)[2:ncol(dat.meta)])
val1=strsplit(treat.xwalk$meta.colnm,"\\_")
treat.xwalk$rep.val=sapply(val1,"[",7)

treat.xwalk$treat=with(treat.xwalk,ifelse(grepl("no_no",meta.colnm),paste0("Treatment 4",rep.val),
                                          ifelse(grepl("no_NH4H2PO4.NaNO3",meta.colnm),paste0("Treatment 3",rep.val),
                                                 ifelse(grepl("ANSo_no",meta.colnm),paste0("Treatment 2",rep.val),
                                                        ifelse(grepl("ANSo_NH4H2PO4.NaNO3",meta.colnm),paste0("Treatment 1",rep.val),NA)))))

dat.16S=read.xlsx(
  paste0(data.path,"Genomics/ISCHudBay_16S.feature_table_final_rarefied_L6.xlsx"))
names(dat.16S)

## Table variables for crosswalk
tables.val=data.frame(Chemistry=PAH.dat$Sample.Name)
val=sapply(strsplit(tables.val$Chemistry," " ),"[",2)
val[val=="Treatment"]=NA
tables.val$Genetic.16s=ifelse(is.na(val),NA,paste("WilliamKennedy",val,sep="_"))
tables.val=base::merge(tables.val,
                       treat.xwalk[,c("meta.colnm","treat")],
                       by.x="Chemistry",by.y="treat",sort=F,all.x=T)
colnames(tables.val)=c(colnames(tables.val)[1:2],"Genetic.MG")

tables.val%>%
  flextable()%>%
  padding(padding=1.5,part="all")%>%
  autofit()%>%
  font(fontname="Times New Roman",part="all")# %>%print("docx")



## 16s analysis -----------------------------------------------------------

str.val=strsplit(dat.16S$Taxon,";")
taxon.16s=data.frame(king=sapply(str.val,"[",1),
                     phyl=sapply(str.val,"[",2),
                     class=sapply(str.val,"[",3),
                     ord=sapply(str.val,"[",4),
                     fam=sapply(str.val,"[",5),
                     genera=sapply(str.val,"[",6))
for(i in 1:ncol(taxon.16s)){
  taxon.16s[,i]=ifelse(taxon.16s[,i]=="Other","Other",sapply(strsplit(taxon.16s[,i],"__"),"[",2))
}
head(taxon.16s)
taxon.16s$Taxon=dat.16S$Taxon
## artifical ID number
taxon.16s$sppID=paste0("SPP",1:nrow(taxon.16s))

dat.16S=merge(dat.16S,taxon.16s,"Taxon")
vars=c(paste("WilliamKennedy",paste0(sort(rep(1:4,3)),letters[1:3]),sep="_"),"sppID")
dat.16S.melt=melt(dat.16S[,vars])

dat.16S.2=dcast(dat.16S.melt,variable~sppID,value.var="value",sum)
dat.16S.2$Sample_TotalSeqs <- rowSums(dat.16S.2[,2:ncol(dat.16S.2)])
## data is reads per 1000 reads (rarefied)
dat.16S.2.gen=dat.16S.2[,2:(ncol(dat.16S.2)-1)]
## Diversity indicies
## proportional abundance
x.val=dat.16S.2$Sample_TotalSeqs# apply(dat.16S.2.gen[,2:ncol(dat.16S.2.gen)],1,sum)
dat.16S.2.gen[,2:ncol(dat.16S.2.gen)]=sweep(dat.16S.2.gen[,2:ncol(dat.16S.2.gen)],1,x.val,"/")

## Diversity (alpha) Indices
# Shannon Index
shannon=-dat.16S.2.gen[,2:ncol(dat.16S.2.gen)]*log(dat.16S.2.gen[,2:ncol(dat.16S.2.gen)])
shannon_H=apply(shannon,1,sum,na.rm=T)

barplot(shannon_H)

# Simpson
simpson_D=1-apply(dat.16S.2.gen[,2:ncol(dat.16S.2.gen)]^2,1,sum,na.rm=T)
simpson_invD=1/apply(dat.16S.2.gen[,2:ncol(dat.16S.2.gen)]^2,1,sum,na.rm=T)

plot(simpson_D)
plot(simpson_invD)
#species richness
richness=vegan::specnumber(dat.16S.2.gen[,2:ncol(dat.16S.2.gen)])

# Pielou's evenness J' #https://www.rpubs.com/roalle/mres_2019
pielou_even=shannon_H/log(richness)

cbind(treat=dat.16S.2$variable,
      shannon_H,
      simpson_D,
      simpson_invD,
      richness,
      pielou_even)

### Top taxa

dat.16S.2.melt=melt(dat.16S.2[,1:(ncol(dat.16S.2)-1)],id.vars = c("variable"))
colnames(dat.16S.2.melt)=c("treat","variable","value")
dat.16S.2.melt.tax=merge(dat.16S.2.melt,taxon.16s,by.x="variable",by.y="sppID",all.x=T)

tmp=dcast(dat.16S.2.melt.tax,treat~fam,value.var = "value",sum)
top.names.fam=names(sort(-rank(colSums(tmp[,2:ncol(tmp)])))[1:4])
top.names.fam2=names(sort(-rank(colSums(tmp[,2:ncol(tmp)])))[1:20])

tmp=dcast(dat.16S.2.melt.tax,treat~genera,value.var = "value",sum)
top.names=names(sort(-rank(colSums(tmp[,2:ncol(tmp)])))[1:5])
top.names.gen2=names(sort(-rank(colSums(tmp[,2:ncol(tmp)])))[1:20])
top.names.gen3=names(sort(-rank(colSums(tmp[,2:ncol(tmp)])))[1:10])
top.names=top.names[top.names!="Other"]                
# subset(taxon.16s2,genera=="SM1A02")

taxon.16s2=taxon.16s
taxon.16s2$gen.group=with(taxon.16s,ifelse(genera%in%top.names,genera,"Other"))                
taxon.16s2$gen.group2=with(taxon.16s,ifelse(genera%in%top.names.gen2,genera,"Other"))                
taxon.16s2$gen.group3=with(taxon.16s,ifelse(genera%in%top.names.gen3,genera,"Other"))                
taxon.16s2$fam.group=with(taxon.16s,ifelse(fam%in%top.names.fam,fam,"Other")) 
taxon.16s2$fam.group2=with(taxon.16s,ifelse(fam%in%top.names.fam2,fam,"Other")) 


dat.16S.2.melt.tax2=merge(dat.16S.2.melt,taxon.16s2,by.x="variable",by.y="sppID",all.x=T)

top4_abund.gen=dcast(dat.16S.2.melt.tax2,treat~gen.group,value.var = "value",sum)
top4_abund.gen=top4_abund.gen[,c("treat","Other",unique(subset(taxon.16s2,gen.group!="Other")$genera))]
top4_abund.gen[,2:ncol(top4_abund.gen)]=sweep(top4_abund.gen[,2:ncol(top4_abund.gen)],1,rep(1000,nrow(top4_abund.gen)),"/")

top20_abund.gen=dcast(dat.16S.2.melt.tax2,treat~gen.group2,value.var = "value",sum)
top20_abund.gen=top20_abund.gen[,c("treat","Other",unique(subset(taxon.16s2,gen.group2!="Other")$genera))]
top20_abund.gen[,2:ncol(top20_abund.gen)]=sweep(top20_abund.gen[,2:ncol(top20_abund.gen)],1,rep(1000,nrow(top20_abund.gen)),"/")

top10_abund.gen=dcast(dat.16S.2.melt.tax2,treat~gen.group3,value.var = "value",sum)
top10_abund.gen=top10_abund.gen[,c("treat","Other",unique(subset(taxon.16s2,gen.group3!="Other")$genera))]
top10_abund.gen[,2:ncol(top10_abund.gen)]=sweep(top10_abund.gen[,2:ncol(top10_abund.gen)],1,rep(1000,nrow(top10_abund.gen)),"/")

top4_abund.fam=dcast(dat.16S.2.melt.tax2,treat~fam.group,value.var = "value",sum)
top4_abund.fam=top4_abund.fam[,c("treat","Other",unique(subset(taxon.16s2,fam.group!="Other")$fam))]
top4_abund.fam[,2:ncol(top4_abund.fam)]=sweep(top4_abund.fam[,2:ncol(top4_abund.fam)],1,rep(1000,nrow(top4_abund.fam)),"/")

top20_abund.fam=dcast(dat.16S.2.melt.tax2,treat~fam.group2,value.var = "value",sum)
top20_abund.fam=top20_abund.fam[,c("treat","Other",unique(subset(taxon.16s2,fam.group2!="Other")$fam))]
top20_abund.fam[,2:ncol(top20_abund.fam)]=sweep(top20_abund.fam[,2:ncol(top20_abund.fam)],1,rep(1000,nrow(top20_abund.fam)),"/")


# png(filename=paste0(plot.path,"PAH_treat_16s_family.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(4,1.5,0.5,3.5),oma=c(1,2,0.75,0.5),lwd=0.5);
layout(matrix(1:2,1,2,byrow = T),widths=c(1,0.25))
ylim.val=c(0,1);by.y=0.1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
cols=c("grey",viridis::plasma(4,alpha=0.5))
xlabs=rev(c("Untreated","Fertilizer","T2 ANSo","T1 ANSo\nFertilizer"))
leg.vals=names(top4_abund.fam[2:ncol(top4_abund.fam)])
x=barplot(t(top4_abund.fam[2:ncol(top4_abund.fam)]),space=c(0,0,0,0.5,0,0,0.5,0,0,0.5),col=cols,axes=F,yaxs="i",width=0.75,ylim=ylim.val,xpd=F)
axis_fun(1,x[c(2,5,8,11)],x[c(2,5,8,11)],xlabs,padj=1,line=-0.75,lwd=1);box(lwd=1)
axis_fun(2,ymaj,ymin,format(ymaj))
mtext(side=1,line=3,"Treatment")
mtext(side=2,line=2.5,"Proportion")

par(mar=c(4,0.5,0.5,0.5))
plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.1,0.5,legend=rev(leg.vals),
       pch=22,pt.bg=rev(cols),pt.cex = 1.5,
       lty=c(NA),lwd=c(0.01),col=rev(cols),
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Family")
dev.off()

# png(filename=paste0(plot.path,"PAH_treat_16s_genera.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(4,1.5,0.5,3.5),oma=c(1,2,0.75,0.5),lwd=0.5);
layout(matrix(1:2,1,2,byrow = T),widths=c(1,0.25))

cols=c("grey",viridis::plasma(4,alpha=0.5))
xlabs=rev(c("Untreated","Fertilizer","T2 ANSo","T1 ANSo\nFertilizer"))
leg.vals=names(top4_abund.gen[2:ncol(top4_abund.gen)])
x=barplot(t(top4_abund.gen[2:ncol(top4_abund.gen)]),space=c(0,0,0,0.5,0,0,0.5,0,0,0.5),col=cols,axes=F,yaxs="i",width=0.75,ylim=ylim.val,xpd=F)
axis_fun(1,x[c(2,5,8,11)],x[c(2,5,8,11)],xlabs,padj=1,line=-0.75,lwd=1);box(lwd=1)
axis_fun(2,ymaj,ymin,format(ymaj))
mtext(side=1,line=3,"Treatment")
mtext(side=2,line=2.5,"Proportion")

par(mar=c(4,0.5,0.5,0.5))
plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.1,0.5,legend=rev(leg.vals),
       pch=22,pt.bg=rev(cols),pt.cex = 1.5,
       lty=c(NA),lwd=c(0.01),col=rev(cols),
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Genera")
dev.off()


# RDA ---------------------------------------------------------------------
library(vegan)

## Top 4 Fam --------------------------------------------------------------
top4_abund.fam
vals=strsplit(as.character(top4_abund.fam$treat),"_")
top4_abund.fam$treat.rep=paste0("T",sapply(vals,"[",2))

total.vals
total.vals$treat.rep=with(total.vals,paste0(Treatment,rep))
total.vals2=subset(total.vals,Treatment!="Control")
total.vals2=dcast(total.vals2,treat.rep+Treatment+rep~param,value.var = "TValue",mean)

gen.vars=c("treat.rep",names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)])
chem.gen.fam=merge(top4_abund.fam[,gen.vars],total.vals2,"treat.rep")

env_st<-chem.gen.fam[,c("alkanes","methly.PAH","T.PAH")]
# abotu<-top4_abund.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]]
# abotu<-decostand(top4_abund.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]],"standardize")
 abotu<-decostand(chem.gen.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]],"hell")
                       
spe.rda <- rda(abotu~., data=env_st)
ord.step.vals=ordiR2step(rda(abotu~1, data=env_st), scope= formula(spe.rda), direction= "both", R2scope=TRUE, pstep=100)

env_st2=env_st[,c("T.PAH","alkanes")]

spe.rda2 <- rda(abotu~., data=env_st2)
summary(spe.rda2, display=NULL)
vif.cca(spe.rda2)

anova.cca(spe.rda2, step=1000)
anova.cca(spe.rda2, by='axis', step=1000)
rslt.terms=anova.cca(spe.rda2, by='terms', step=1000)
rslt.terms
(R2adj <- RsquareAdj(spe.rda2)$adj.r.squared)

rslt.terms=data.frame(rslt.terms)
rslt.terms$variables=rownames(rslt.terms)
rownames(rslt.terms)=1:nrow(rslt.terms)
rslt.terms=rslt.terms[,c(5,1:4)]
colnames(rslt.terms)=c("variable","df","var","fvalue","pvalue")

rslt.terms%>%
  flextable%>%
  colformat_double(j=3,digits=3,na_str = " ")%>%
  colformat_double(j=4,digits=2,na_str = " ")%>%
  compose(j="pvalue",i=~pvalue<0.05,value=as_paragraph('< 0.05'))%>%
  compose(j="pvalue",i=~pvalue<0.01,value=as_paragraph('< 0.01'))%>%
  italic(j="pvalue",i=~pvalue<0.05)%>%
  compose(j="variable",i=~variable=="ATemp.DegC",value=as_paragraph('Temperature'))%>%
  compose(j="variable",i=~variable=="SPC.uScm",value=as_paragraph('Specific Conductance'))%>%
  compose(j="variable",i=~variable=="SDD.cm",value=as_paragraph('Secchi'))%>%
  compose(j="variable",i=~variable=="TN.mgL",value=as_paragraph('Total Nitrogen'))%>%
  compose(j="variable",i=~variable=="DP.ugL",value=as_paragraph('Dissolved Phosphorus'))%>%
  compose(j="variable",i=~variable=="TN_TP",value=as_paragraph('TN:TP (molar rato)'))%>%
  set_header_labels("variable"="Parameter",
                    "df"="DF",
                    "var" = "Variance",
                    "fvalue" = "F",
                    "pvalue"="\u03C1-value")%>%
  width(width=c(1.5,0.5,0.75,0.75,0.75))%>%
  align(align="center",part="header")%>%
  align(j=2:5,align="center",part="body")%>%
  padding(padding=1.5,part="all")%>%
  footnote(j=1,ref_symbols = " ",part="header",
           value=as_paragraph("R\u00B2 Adj: ",format(round(R2adj,2),nsmall=2)))%>%
  font(fontname="Times New Roman",part="all")%>%
  fontsize(size=11,part="all")%>%
  fontsize(size=12,part="header")# %>%print("docx")

eig <- spe.rda2$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
eig.pca

# png(filename=paste0(plot.path,"PAHtreat_16s_RDA_fam_scree.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,1,0.25,0.5));

ylim.val=c(0,120);by.y=25;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2);#set y limit and delineates the major and minor ticks
x=barplot(eig.pca$variance,ylim=ylim.val,col="white",border=0,yaxt="n")# inital plot to get the measurements
abline(h=ymaj,lty=3,col="grey")#makes vertical lines from y axis
x=barplot(eig.pca$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)# the real plot that matters
lines(x,eig.pca$cumvariance,col="indianred1",lwd=2)# adds the cumulative variance for each factor
points(x,eig.pca$cumvariance,pch=21,bg="indianred1",cex=1.25)
abline(h=80,lty=2,col="red",lwd=2)
axis_fun(1,line=-0.5,x,x,seq(1,length(x),1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Components")
mtext(side=2,line=1.75,"Percentage of Variances")
legend.text=c("Absolute","Cumulative");#helper vaiable for legend
pt.col=c("grey","indianred1")#helper vaiable for legend
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

plot(spe.rda2)
plot(spe.rda2,scaling="sites")
plot(spe.rda2,scaling="species")
plot(spe.rda2,scaling="symmetric")

total.vals2
scrs2<-scores(spe.rda2,display=c("sites","species"),choices=c(1,2,3),scaling="sites");
scrs.arrows<-scores(spe.rda2,choices=c(1,2,3),"bp",scaling="sites");
points(scrs2$sites[,1:2],pch=21,bg=cols,cex=2)
labs=rownames(scrs$species)
labs.arrows=rownames(scrs.arrows)

x.labs=c("T1 ANSo Fertilizer", "T2 ANSo", "Fertilizer", "Untreated")
cols=viridis::plasma(4,alpha=0.5)
cols2=c(rep(cols[1],3),rep(cols[2],3),rep(cols[3],3),rep(cols[4],3))

# png(filename=paste0(plot.path,"PAHtreat_16s_RDA_fam.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,2,0.25,0.5));
layout(matrix(1:2,1,2,byrow = T),widths=c(1,0.4))

xlim.val=c(-0.8,0.6);by.x=0.4;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-0.5,0.5);by.y=0.25;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs2$sites[,c(1,2)],pch=21,bg=cols2,col="black",cex=2,lwd=0.5)
arrows(0,0,scrs.arrows[,1],scrs.arrows[,2],length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
text(scrs.arrows[,1],scrs.arrows[,2],labels=labs.arrows,cex=0.75,font=1,col="red")
with(scrs2,text(species[,1],species[,2],labels=labs,cex=0.75,font=3,col="blue"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("RDA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.25,paste0("RDA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance

plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.25,0.5,legend=x.labs,
       pch=21,pt.bg=cols,pt.cex = 1.5,
       lty=c(NA),lwd=c(0.01),col=rev(cols),
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Treatment")
dev.off()



## Top 10 genus level -----------------------------------------------------
top10_abund.gen
vals=strsplit(as.character(top10_abund.gen$treat),"_")
top10_abund.gen$treat.rep=paste0("T",sapply(vals,"[",2))

total.vals
total.vals$treat.rep=with(total.vals,paste0(Treatment,rep))
total.vals2=subset(total.vals,Treatment!="Control")
total.vals2=dcast(total.vals2,treat.rep+Treatment+rep~param,value.var = "TValue",mean)

gen.vars=c("treat.rep",names(top10_abund.gen)[2:(ncol(top10_abund.gen)-1)])
chem.gen=merge(top10_abund.gen[,gen.vars],total.vals2,"treat.rep")

env_st<-chem.gen[,c("alkanes","methly.PAH","T.PAH")]
# abotu<-top4_abund.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]]
# abotu<-decostand(top4_abund.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]],"standardize")
abotu<-decostand(chem.gen[,names(top10_abund.gen)[2:(ncol(top10_abund.gen)-1)]],"hell")

spe.rda <- rda(abotu~., data=env_st)
ord.step.vals=ordiR2step(rda(abotu~1, data=env_st), scope= formula(spe.rda), direction= "both", R2scope=TRUE, pstep=100)

env_st2=env_st[,c("T.PAH","alkanes")]

spe.rda2_genus <- rda(abotu~., data=env_st2)
summary(spe.rda2_genus, display=NULL)
vif.cca(spe.rda2_genus)

anova.cca(spe.rda2_genus, step=1000)
anova.cca(spe.rda2_genus, by='axis', step=1000)
rslt.terms=anova.cca(spe.rda2_genus, by='terms', step=1000)
rslt.terms
(R2adj <- RsquareAdj(spe.rda2_genus)$adj.r.squared)

rslt.terms=data.frame(rslt.terms)
rslt.terms$variables=rownames(rslt.terms)
rownames(rslt.terms)=1:nrow(rslt.terms)
rslt.terms=rslt.terms[,c(5,1:4)]
colnames(rslt.terms)=c("variable","df","var","fvalue","pvalue")

rslt.terms%>%
  flextable%>%
  colformat_double(j=3,digits=3,na_str = " ")%>%
  colformat_double(j=4,digits=2,na_str = " ")%>%
  compose(j="pvalue",i=~pvalue<0.05,value=as_paragraph('< 0.05'))%>%
  compose(j="pvalue",i=~pvalue<0.01,value=as_paragraph('< 0.01'))%>%
  italic(j="pvalue",i=~pvalue<0.05)%>%
  set_header_labels("variable"="Parameter",
                    "df"="DF",
                    "var" = "Variance",
                    "fvalue" = "F",
                    "pvalue"="\u03C1-value")%>%
  width(width=c(1.5,0.5,0.75,0.75,0.75))%>%
  align(align="center",part="header")%>%
  align(j=2:5,align="center",part="body")%>%
  padding(padding=1.5,part="all")%>%
  footnote(j=1,ref_symbols = " ",part="header",
           value=as_paragraph("R\u00B2 Adj: ",format(round(R2adj,2),nsmall=2)))%>%
  font(fontname="Times New Roman",part="all")%>%
  fontsize(size=11,part="all")%>%
  fontsize(size=12,part="header")# %>%print("docx")

eig <- spe.rda2_genus$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca_genus <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
eig.pca_genus

# png(filename=paste0(plot.path,"PAHtreat_16s_RDA_fam_scree.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,1,0.25,0.5));

ylim.val=c(0,120);by.y=25;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2);#set y limit and delineates the major and minor ticks
x=barplot(eig.pca_genus$variance,ylim=ylim.val,col="white",border=0,yaxt="n")# inital plot to get the measurements
abline(h=ymaj,lty=3,col="grey")#makes vertical lines from y axis
x=barplot(eig.pca_genus$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)# the real plot that matters
lines(x,eig.pca_genus$cumvariance,col="indianred1",lwd=2)# adds the cumulative variance for each factor
points(x,eig.pca_genus$cumvariance,pch=21,bg="indianred1",cex=1.25)
abline(h=80,lty=2,col="red",lwd=2)
axis_fun(1,line=-0.5,x,x,seq(1,length(x),1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Components")
mtext(side=2,line=1.75,"Percentage of Variances")
legend.text=c("Absolute","Cumulative");#helper vaiable for legend
pt.col=c("grey","indianred1")#helper vaiable for legend
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

plot(spe.rda2_genus)
plot(spe.rda2_genus,scaling="sites")
plot(spe.rda2_genus,scaling="species")
plot(spe.rda2_genus,scaling="symmetric")


total.vals2
scrs2<-scores(spe.rda2_genus,display=c("sites","species"),choices=c(1,2,3),scaling="sites");
scrs.arrows<-scores(spe.rda2_genus,choices=c(1,2,3),"bp",scaling="sites");
points(scrs2$sites[,1:2],pch=21,bg=cols,cex=2)
labs=rownames(scrs2$species)
labs.arrows=rownames(scrs.arrows)

x.labs=c("T1 ANSo Fertilizer", "T2 ANSo", "Fertilizer", "Untreated")
cols=viridis::plasma(4,alpha=0.5)
cols2=c(rep(cols[1],3),rep(cols[2],3),rep(cols[3],3),rep(cols[4],3))

# png(filename=paste0(plot.path,"PAHtreat_16s_RDA_Top10genus.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,2,0.25,0.5));
layout(matrix(1:2,1,2,byrow = T),widths=c(1,0.4))

xlim.val=c(-1,0.4);by.x=0.2;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-0.5,0.25);by.y=0.25;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs2$sites[,c(1,2)],pch=21,bg=cols2,col="black",cex=2,lwd=0.5)
arrows(0,0,scrs.arrows[,1],scrs.arrows[,2],length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
text(scrs.arrows[,1],scrs.arrows[,2],labels=labs.arrows,cex=0.75,font=1,col="red")
with(scrs2,text(species[,1],species[,2],labels=labs,cex=0.75,font=3,col="blue"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("RDA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.25,paste0("RDA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance


plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.25,0.5,legend=x.labs,
       pch=21,pt.bg=cols,pt.cex = 1.5,
       lty=c(NA),lwd=c(0.01),col=rev(cols),
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Treatment")
dev.off()






## Top 4 genus ------------------------------------------------------------
top4_abund.gen
vals=strsplit(as.character(top4_abund.gen$treat),"_")
top4_abund.gen$treat.rep=paste0("T",sapply(vals,"[",2))

total.vals
total.vals$treat.rep=with(total.vals,paste0(Treatment,rep))
total.vals2=subset(total.vals,Treatment!="Control")
total.vals2=dcast(total.vals2,treat.rep+Treatment+rep~param,value.var = "TValue",mean)

gen.vars=c("treat.rep",names(top4_abund.gen)[2:(ncol(top4_abund.gen)-1)])
chem.gen=merge(top4_abund.gen[,gen.vars],total.vals2,"treat.rep")

env_st<-chem.gen[,c("alkanes","methly.PAH","T.PAH")]
# abotu<-top4_abund.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]]
# abotu<-decostand(top4_abund.fam[,names(top4_abund.fam)[2:(ncol(top4_abund.fam)-1)]],"standardize")
abotu<-decostand(chem.gen[,names(top4_abund.gen)[2:(ncol(top4_abund.gen)-1)]],"hell")

env_st2=env_st[,c("T.PAH","alkanes")]

spe.rda2_genus <- rda(abotu~., data=env_st2)
summary(spe.rda2_genus, display=NULL)
vif.cca(spe.rda2_genus)

anova.cca(spe.rda2_genus, step=1000)
anova.cca(spe.rda2_genus, by='axis', step=1000)
rslt.terms=anova.cca(spe.rda2_genus, by='terms', step=1000)
rslt.terms
(R2adj <- RsquareAdj(spe.rda2_genus)$adj.r.squared)

rslt.terms=data.frame(rslt.terms)
rslt.terms$variables=rownames(rslt.terms)
rownames(rslt.terms)=1:nrow(rslt.terms)
rslt.terms=rslt.terms[,c(5,1:4)]
colnames(rslt.terms)=c("variable","df","var","fvalue","pvalue")

rslt.terms%>%
  flextable%>%
  colformat_double(j=3,digits=3,na_str = " ")%>%
  colformat_double(j=4,digits=2,na_str = " ")%>%
  compose(j="pvalue",i=~pvalue<0.05,value=as_paragraph('< 0.05'))%>%
  compose(j="pvalue",i=~pvalue<0.01,value=as_paragraph('< 0.01'))%>%
  italic(j="pvalue",i=~pvalue<0.05)%>%
  set_header_labels("variable"="Parameter",
                    "df"="DF",
                    "var" = "Variance",
                    "fvalue" = "F",
                    "pvalue"="\u03C1-value")%>%
  width(width=c(1.5,0.5,0.75,0.75,0.75))%>%
  align(align="center",part="header")%>%
  align(j=2:5,align="center",part="body")%>%
  padding(padding=1.5,part="all")%>%
  footnote(j=1,ref_symbols = " ",part="header",
           value=as_paragraph("R\u00B2 Adj: ",format(round(R2adj,2),nsmall=2)))%>%
  font(fontname="Times New Roman",part="all")%>%
  fontsize(size=11,part="all")%>%
  fontsize(size=12,part="header")# %>%print("docx")

eig <- spe.rda2_genus$CA$eig
variance <- eig*100/sum(eig)
cumvar <- cumsum(variance)
eig.pca_genus <- data.frame(eig = eig, variance = variance,cumvariance = cumvar)
eig.pca_genus

# png(filename=paste0(plot.path,"PAHtreat_16s_RDA_fam_scree.png"),width=5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,1,0.25,0.5));

ylim.val=c(0,120);by.y=25;ymaj=seq(ylim.val[1],100,by.y);ymin=seq(ylim.val[1],100,by.y/2);#set y limit and delineates the major and minor ticks
x=barplot(eig.pca_genus$variance,ylim=ylim.val,col="white",border=0,yaxt="n")# inital plot to get the measurements
abline(h=ymaj,lty=3,col="grey")#makes vertical lines from y axis
x=barplot(eig.pca_genus$variance,ylim=ylim.val,col="grey",yaxt="n",add=T)# the real plot that matters
lines(x,eig.pca_genus$cumvariance,col="indianred1",lwd=2)# adds the cumulative variance for each factor
points(x,eig.pca_genus$cumvariance,pch=21,bg="indianred1",cex=1.25)
abline(h=80,lty=2,col="red",lwd=2)
axis_fun(1,line=-0.5,x,x,seq(1,length(x),1),0.7)
axis_fun(2,ymaj,ymin,ymaj,0.75);box(lwd=1)
mtext(side=1,line=1.5,"Components")
mtext(side=2,line=1.75,"Percentage of Variances")
legend.text=c("Absolute","Cumulative");#helper vaiable for legend
pt.col=c("grey","indianred1")#helper vaiable for legend
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col=c("black",pt.col[2]),lty=c(0,1),lwd=1.5,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,text.col="white")
legend("topleft",legend=legend.text,pch=c(22,21),pt.bg=pt.col,col="black",lty=0,lwd=0.5,pt.cex=1,ncol=2,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5)
dev.off()

plot(spe.rda2_genus)
plot(spe.rda2_genus,scaling="sites")
plot(spe.rda2_genus,scaling="species")
plot(spe.rda2_genus,scaling="symmetric")

total.vals2
scrs2<-scores(spe.rda2_genus,display=c("sites","species"),choices=c(1,2,3),scaling="sites");
scrs.arrows<-scores(spe.rda2_genus,choices=c(1,2,3),"bp",scaling="sites");
points(scrs2$sites[,1:2],pch=21,bg=cols,cex=2)
labs=rownames(scrs2$species)
labs.arrows=rownames(scrs.arrows)

x.labs=c("T1 ANSo Fertilizer", "T2 ANSo", "Fertilizer", "Untreated")
cols=viridis::plasma(4,alpha=0.5)
cols2=c(rep(cols[1],3),rep(cols[2],3),rep(cols[3],3),rep(cols[4],3))

# png(filename=paste0(plot.path,"PAHtreat_16s_RDA_Top4genus.png"),width=6.5,height=4.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,2,0.25,0.5));
layout(matrix(1:2,1,2,byrow = T),widths=c(1,0.4))

xlim.val=c(-1,0.4);by.x=0.2;xmaj=c(0,seq(xlim.val[1],xlim.val[2],by.x));xmin=seq(xlim.val[1],xlim.val[2],by.x/2);
ylim.val=c(-0.4,0.5);by.y=0.2;ymaj=c(0,seq(ylim.val[1],ylim.val[2],by.y));ymin=seq(ylim.val[1],ylim.val[2],by.y/2);
plot(xlim.val,ylim.val,type="n",yaxt="n",xaxt="n",ylab=NA,xlab=NA);
abline(h=0,v=0,lty=3,col="grey");
points(scrs2$sites[,c(1,2)],pch=21,bg=cols2,col="black",cex=2,lwd=0.5)
arrows(0,0,scrs.arrows[,1],scrs.arrows[,2],length = 0.05, angle = 15, code = 2,col="indianred1",lwd=1.5);# makes the arrows
text(scrs.arrows[,1],scrs.arrows[,2],labels=labs.arrows,cex=0.75,font=1,col="red")
with(scrs2,text(species[,1],species[,2],labels=labs,cex=0.75,font=3,col="blue"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj),1); #adds x axis ticks
axis_fun(2,ymaj,ymin,format(ymaj),1); #adds y axis ticks
mtext(side=1,line=1.8,paste0("RDA 1 (",round(eig.pca$variance[1],1),"%)"));#adds x axis label with percent variance
mtext(side=2,line=2.25,paste0("RDA 2 (",round(eig.pca$variance[2],1),"%)"));#adds y axis label with percent variance


plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.25,0.5,legend=x.labs,
       pch=21,pt.bg=cols,pt.cex = 1.5,
       lty=c(NA),lwd=c(0.01),col=rev(cols),
       ncol=1,cex=0.75,bty="n",y.intersp=1,x.intersp=1,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj=0,title="Treatment")
dev.off()
