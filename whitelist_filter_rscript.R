library(data.table, quietly=T)
gList<-fread("NEJM_2017_genes_01262020.txt")
whitelist.mis<-fread("CHIP_missense_vars_agb_01262020.txt")
whitelist.splice<-fread("CHIP_splice_vars_agb_01262020.txt")
whitelist.LoF<-fread("CHIP_nonsense_FS_vars_agb_01262020.txt")

args <- commandArgs(trailingOnly=TRUE)
sample_id<-args[2]
vars1<-fread("annovar_out.hg38_multianno.txt")
vars<-transform(vars1, Sample=sample_id)

colnames=c("Chr","Start","End","Ref","Alt","Func.refGene","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
names(vars)[1:10]<-colnames
rm(colnames)

varsOI<-vars[Gene.refGene%in%gList$Gene,]


varsOI.func<-varsOI[Gene.refGene%in%gList$Gene&
                      (grepl("exonic",varsOI$Func.refGene,fixed=T)|
                         grepl("splicing",varsOI$Func.refGene,fixed=T)),]
varsOI.func<-merge(varsOI.func,gList[,c("Gene", "Accession")],by.x="Gene.refGene", by.y="Gene")

#Func.refGene # (exonic, splicing or some wierd merged combination)
#GeneDetail.refGene # NM_017940:exon16:c.1380-2A>G
#ExonicFunc.refGene #"nonsynonymous SNV"
#AAChange.refGene #"ASXL1:NM_015338:exon11:c.G1718A:p.R573Q"
#Accession          "NM_015338"
extractTranscript<-function(GeneDetail,AAChange,Accession){
  splice<-grep(Accession,strsplit(GeneDetail,",")[[1]], value=T, fixed=T)
  nonsyn<-grep(Accession,strsplit(AAChange,",")[[1]], value=T, fixed=T)
  if(length(splice)>0){return(splice)} else {
    if(length(nonsyn)>0){return(nonsyn)} else{
      return("nan")}
  }
}

#Apply function
aa<-apply(varsOI.func[,c('GeneDetail.refGene', 'AAChange.refGene','Accession')],
          1, function(x) {extractTranscript(x[1],x[2],x[3]) })
aa2<-lapply(aa, `[[`, 1)
varsOI.func$transcriptOI<-unlist(aa2)
varsOI.func<-varsOI.func[varsOI.func$transcriptOI!="nan",]
rm(aa,aa2)

extractNonsyn<-function(transcriptOI){
  protChange<-grep("p.",strsplit(transcriptOI,":")[[1]], value=T,fixed=T)
  if(length(protChange)>0){return(gsub("p\\.","",protChange))} else {return("nan")}
}

bb<-apply(varsOI.func[,c('transcriptOI')],1, function(x) {extractNonsyn(x[1]) })
bb2<-lapply(bb, `[[`, 1)
varsOI.func$NonsynOI<-unlist(bb2)
rm(bb,bb2)

#Annotate and filter
varsOI.func<-transform(varsOI.func,whitelist=F,
                       wl.mis=F,wl.lof=F,wl.splice=F,wl.exception=F,
                       manualreview=F)
#1) Handle missense vars
vmis<-varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"
vmis_wl<-paste(varsOI.func$Gene.refGene,varsOI.func$NonsynOI,sep="_")%in% paste(whitelist.mis$Gene,whitelist.mis$AAChange,sep="_")
varsOI.func[vmis&vmis_wl,"whitelist"]=T
varsOI.func[vmis&vmis_wl,"wl.mis"]=T

#2) Handle LoF and frame shift vars
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vLOFgene<-varsOI.func$Gene.refGene%in%whitelist.LoF$Gene #genes are Lof Genes
varsOI.func[(vlof|vFS)&vLOFgene,"whitelist"]=T
varsOI.func[(vlof|vFS)&vLOFgene,"wl.lof"]=T

#3) Handle Splice vars
#possible splicing variant
vSplice<-grepl("splicing",varsOI.func$Func.refGene, fixed=T)
#genes are splice Genes
vSplicegene<-varsOI.func$Gene.refGene%in%whitelist.splice$Gene
#confirm that the splicing refers to the correct transcript
vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
                                1, function(x) {grepl(x[2],x[1],fixed=T)})
varsOI.func[vSplice&vSplicegene&vSpliceCorrectTranscript,"whitelist"]=T
varsOI.func[vSplicegene&vSplicegene&vSpliceCorrectTranscript,"wl.splice"]=T
#If not the correct transcript, flag for manual review
varsOI.func[(vSplice&vSplicegene) & (!vSpliceCorrectTranscript),"manualreview"]=T

#4) Handle the following exceptions
#ASXL1	Frameshift/nonsense/splice-site in exon 11-12
#LoF
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon11<-grepl("exon11",varsOI.func$transcriptOI, fixed=T)
vexon12<-grepl("exon12",varsOI.func$transcriptOI, fixed=T)
asxl1Exception<-(varsOI.func$Gene.refGene=="ASXL1")&(vlof|vFS)&(vexon11|vexon12)
varsOI.func[asxl1Exception,"whitelist"]=T
varsOI.func[asxl1Exception,"wl.lof"]=T
varsOI.func[asxl1Exception,"wl.exception"]=T
#Splicing
vSplice<-grepl("splicing",varsOI.func$Func.refGene, fixed=T)
vSpliceCorrectTranscript<-apply(varsOI.func[,c('GeneDetail.refGene','Accession')],
                                1, function(x) {grepl(x[2],x[1],fixed=T)})
vexon11splice<-grepl("exon11",varsOI.func$GeneDetail.refGene, fixed=T)
vexon12splice<-grepl("exon12",varsOI.func$GeneDetail.refGene, fixed=T)
asxl1ExceptionSplice<-(varsOI.func$Gene.refGene=="ASXL1")&
  vSplice&vSpliceCorrectTranscript&
  (vexon11splice|vexon12splice)
varsOI.func[asxl1ExceptionSplice,"whitelist"]=T
varsOI.func[asxl1ExceptionSplice,"wl.splice"]=T
varsOI.func[asxl1ExceptionSplice,"wl.exception"]=T
#ASXL2	Frameshift/nonsense/splice-site in exon 11-12
#Lof
asxl2Exception<-(varsOI.func$Gene.refGene=="ASXL2")&(vlof|vFS)&(vexon11|vexon12)
varsOI.func[asxl2Exception,"whitelist"]=T
varsOI.func[asxl2Exception,"wl.lof"]=T
varsOI.func[asxl2Exception,"wl.exception"]=T
#Splice
asxl2ExceptionSplice<-(varsOI.func$Gene.refGene=="ASXL2")&
  vSplice&vSpliceCorrectTranscript&
  (vexon11splice|vexon12splice)
varsOI.func[asxl2ExceptionSplice,"whitelist"]=T
varsOI.func[asxl2ExceptionSplice,"wl.splice"]=T
varsOI.func[asxl2ExceptionSplice,"wl.exception"]=T

#PPM1D	Frameshift/nonsense in exon 5 or 6
vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vexon5<-grepl("exon5",varsOI.func$transcriptOI, fixed=T)
vexon6<-grepl("exon6",varsOI.func$transcriptOI, fixed=T)
ppm1dException<-(varsOI.func$Gene.refGene=="PPM1D")&(vlof|vFS)&(vexon5|vexon6)
varsOI.func[ppm1dException,"whitelist"]=T
varsOI.func[ppm1dException,"wl.lof"]=T
varsOI.func[ppm1dException,"wl.exception"]=T

#TET2	missense mutations in catalytic domains (p.1104-1481 and 1843-2002)
TETidx<-which(varsOI.func$Gene.refGene=="TET2"&
                varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                nchar(varsOI.func$NonsynOI)==6)
for(i in TETidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,5))
  if((AApos>=1104&AApos<=1481)|(AApos>=1843&AApos<=2002))
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}

#CBL	RING finger missense p.381-421
CBLidx<-which(varsOI.func$Gene.refGene=="CBL"&
                varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                nchar(varsOI.func$NonsynOI)==5)
for(i in CBLidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
  if(AApos>=381&AApos<=421)
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}
#CBLB	RING finger missense p.372-412
CBLBidx<-which(varsOI.func$Gene.refGene=="CBLB"&
                 varsOI.func$ExonicFunc.refGene=="nonsynonymous SNV"&
                 nchar(varsOI.func$NonsynOI)==5)
for(i in CBLBidx){
  AApos<-as.numeric(substr(varsOI.func$NonsynOI[i],2,4))
  if(AApos>=372&AApos<=412)
  {
    varsOI.func[i,"whitelist"]=T
    varsOI.func[i,"wl.mis"]=T
    varsOI.func[i,"wl.exception"]=T
  }
}

#5) flag remaining exceptions for manual review

vlof<-grepl("X",varsOI.func$NonsynOI, fixed=T) #stop gain or stop loss
vFS<-grepl("fs",varsOI.func$NonsynOI, fixed=T) #frameshift
vSplice<-grepl("splicing",varsOI.func$Func.refGene)

#GATA3	Frameshift/nonsense/splice-site ZNF domain
varsOI.func[(vlof|vFS|vSplice)&(varsOI.func$Gene.refGene=="GATA3"),"manualreview"]=T
#CREBBP	S1680del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift deletion" &
              (varsOI.func$Gene.refGene=="CREBBP"),"manualreview"]=T
#CSF3R	truncating c.741-791
varsOI.func[(vlof|vFS|vSplice)&(varsOI.func$Gene.refGene=="CSF3R"),"manualreview"]=T
#DNMT3A	F732del,	F752del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift deletion" &
              (varsOI.func$Gene.refGene=="DNMT3A"),"manualreview"]=T
#EP300	VF1148_1149del
varsOI.func[varsOI.func$ExonicFunc.refGene=="nonframeshift deletion" &
              (varsOI.func$Gene.refGene=="EP300"),"manualreview"]=T
#FLT3	FY590-591GD,	del835
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
               varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="FLT3"),"manualreview"]=T
#JAK2	del/ins537-539L, del/ins538-539L, del/ins540-543MK, del/ins540-544MK, del542-543, del543-544,	ins11546-547
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
               varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="JAK2"),"manualreview"]=T
#KDM6A	del419
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion")
            &(varsOI.func$Gene.refGene=="KDM6A"),"manualreview"]=T
#KIT	ins503,	del560,	del579,	del551-559
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
               varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="KIT"),"manualreview"]=T
#MPL	del513 W515-518KT
varsOI.func[(varsOI.func$ExonicFunc.refGene=="nonframeshift deletion"|
               varsOI.func$ExonicFunc.refGene=="nonframeshift insertion")
            &(varsOI.func$Gene.refGene=="MPL"),"manualreview"]=T
#NPM1	Frameshift p.W288fs (insertion at c.859_860, 860_861, 862_863, 863_864)
varsOI.func[(varsOI.func$ExonicFunc.refGene=="frameshift deletion"|
               varsOI.func$ExonicFunc.refGene=="frameshift insertion")
            &(varsOI.func$Gene.refGene=="NPM1"),"manualreview"]=T


#########
#New code May 11 2020
#Add LoF filtering
#Function to extract number
library(stringr)
numextract <- function(string){
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

#Function to extract amino acid length of protein
library("biomaRt")
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
protein = getSequence(id=gList$Accession,
                      type="refseq_mrna",
                      seqType="peptide",
                      mart=mart)
protein=transform(protein, aalen=nchar(protein$peptide))
varsOI.func2<-merge(varsOI.func,protein[!duplicated(protein$refseq_mrna),-1],
                    by.x="Accession",by.y="refseq_mrna",all.x=T)

varsOI.func2$aapos=-9
varsOI.func2$aapos=numextract(varsOI.func2$NonsynOI)
varsOI.func2$aafirst10pctPeptide=varsOI.func2$aalen *0.1
varsOI.func2$aalast10pctPeptide=varsOI.func2$aalen *0.9

#Flag first and last 10% of protein
varsOI.func2$LOFfirst10pct=varsOI.func2$wl.lof&   (as.numeric(varsOI.func2$aapos )<varsOI.func2$aafirst10pctPeptide)
varsOI.func2$LOFlast10pct=varsOI.func2$wl.lof&   (as.numeric(varsOI.func2$aapos ) > varsOI.func2$aalast10pctPeptide)

check_vars <- data.frame(Sample=sample_id,
                         total_num_variants=length(varsOI.func2$Sample),
                         total_num_whitelist=length(varsOI.func2[whitelist==T,]$Sample),
                         total_num_manualreview=length(varsOI.func2[manualreview==T,]$Sample))

#write out files
write.csv(check_vars,paste(sample_id, ".annovar.varsOI.check.csv", sep=""))
write.csv(varsOI.func2,paste(sample_id, ".annovar.varsOI.allvariants.csv", sep=""))
write.csv(varsOI.func2[whitelist==T,],paste(sample_id, ".annovar.varsOI.wl.csv", sep=""))
write.csv(varsOI.func[manualreview==T,],paste(sample_id, ".annovar.varsOI.manualreview.csv", sep=""))
