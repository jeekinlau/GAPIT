#'
#' GAPIT.SUPER.GS
#'
#' @description 
#' Perform GPS with SUPER and Compress method.
#'
#' @param Y Phenotype data.frame,
#' @param GD = NULL,
#' @param GM = NULL,
#' @param KI = NULL,
#' @param Z = NULL,
#' @param CV = NULL,
#' @param GK = NULL,
#' @param kinship.algorithm = NULL,
#' @param bin.from = 10000,
#' @param bin.to = 10000,
#' @param bin.by = 1000,
#' @param inclosure.from = 10,
#' @param inclosure.to = 10,
#' @param inclosure.by = 10,
#' @param group.from = 1000000,
#' @param group.to = 1000000,
#' @param group.by = 10,
#' @param kinship.cluster = "average", 
#' @param kinship.group = 'Mean',
#' @param PCA.total = 0,
#' @param GT = NULL,
#' @param PC = NULL,
#' @param GI = NULL,
#' @param Timmer  =  NULL, 
#' @param Memory  =  NULL,
#' @param model = "",
#' @param sangwich.top = NULL,
#' @param sangwich.bottom = NULL,
#' @param QC = TRUE,
#' @param GTindex = NULL,
#' @param LD = 0.05,
#' @param file.output = TRUE,
#' @param cutOff = 0.01
#'
#'
#' @author Zhiwu Zhang and Jiabo Wang
#'
#'
#' @export
`GAPIT.SUPER.GS`<-
function(Y,
         GD = NULL,
         allGD=NULL,
         GM = NULL,
         KI = NULL,
         Z = NULL,
         CV = NULL,
         allCV=NULL,
         GK = NULL,
         kinship.algorithm = NULL,
         bin.from = 10000,
         bin.to = 10000,
         bin.by = 1000,
         inclosure.from = 10,
         inclosure.to = 10,
         inclosure.by = 10,
				 group.from = 1000000,
				 group.to = 1000000,
				 group.by = 10,
				 kinship.cluster = "average", 
				 kinship.group = 'Mean',
				 PCA.total = 0,
         GT = NULL,
				 PC = NULL,
				 GI = NULL,
				 Timmer = NULL, 
				 Memory = NULL,
				 model = "",
				 sangwich.top = NULL,
				 sangwich.bottom = NULL,
				 QC = TRUE,
         QTN.gs=NULL,
				 GTindex = NULL,
				 LD = 0.05,
				 file.output = TRUE,
         GAPIT3.output=TRUE,
         CV.Extragenetic=0,
				 cutOff = 0.01
                        ){
 
#Object: To perform GPS with SUPER and Compress method
#Designed by Zhiwu Zhang
#Writen by Jiabo Wang
#Last update: Novber 6, 2015 		
######################################################
print("--------------------- Welcome to GAPIT SUPER GS----------------------------")
Timmer=GAPIT.Timmer(Infor="GAPIT.SUPER.GS")
Memory=GAPIT.Memory(Infor="GAPIT.SUPER.GS")
#  if(!require(EMMREML)) install.packages("EMMREML")
#  library(EMMREML)

shortcut=FALSE
LL.save=1e10
# print(head(Y))
#In case of null Y and null GP, return genotype only  
thisY=Y[,2]
thisY=thisY[!is.na(thisY)]
name.of.trait=colnames(Y)[2]
if(length(thisY) <3){
 shortcut=TRUE
 }else{
  if(stats::var(thisY) ==0) shortcut=TRUE
}
if(shortcut){
print(paste("Y is empty. No GWAS/GS performed for ",name.of.trait,sep=""))
return (list(compression=NULL,kinship.optimum=NULL, kinship=KI,PC=PC,GWAS=NULL, GPS=NULL,Pred=NULL, REMLs=NULL,Timmer=Timmer,Memory=Memory))
}
print("------------Examining data (QC)------------------------------------------")
# if(is.null(Y)) stop ("GAPIT says: Phenotypes must exist.")
if(is.null(KI)&missing(GD) & kinship.algorithm!="SUPER") stop ("GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created.")
if(is.null(GD) & is.null(GT)) {
	GT=as.matrix(Y[,1])
	GD=matrix(1,nrow(Y),1)
  rownames(GD)=as.character(GT)
  GI=as.data.frame(matrix(0,1,3) )
  colnames(GI)=c("SNP","Chromosome","Position")
}
# print(cbind(CV,PC))
# if(PCA.total>0&!is.null(CV))CV=GAPIT.CVMergePC(CV,PC)
# if(PCA.total>0&is.null(CV))CV=PC

if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & is.null(Z)){
taxa=as.character(Y[,1])
Z=as.data.frame(diag(1,nrow(Y)))
Z=rbind(taxa,Z)
taxa=c('Taxa',as.character(taxa))
Z=cbind(taxa,Z)
}
if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & !is.null(Z))
{
  if(nrow(Z)-1<nrow(Y)) Z=GAPIT.ZmatrixFormation(Z=Z,Y=Y)
}
# noCV=FALSE
# if(is.null(CV)){
#   noCV=TRUE
#   if(!is.null(GD))
#   {
#     taxa.gd=as.character(rownames(GD))
#     CV=cbind(as.character(taxa.gd),as.data.frame(GD[,1]))
#     }else{
#       if(!is.null(KI))
#       {
#         CV=KI[,1:2]
#       }else{
#         CV=Y[,1:2]  
#       }
#     }
  
#   CV[,2]=1
#   colnames(CV)=c("taxa","overall")
# }
# print(CV)
#Remove duplicat and integragation of data
print("QC is in process...")

my_allCV=allCV

# if(QC)
# {
# # print(head(Y))
# # print(Z[,1:3])

#   qc <- GAPIT.QC2(Y=Y,KI=KI, GT=GT,CV=CV,Z=Z,GK=GK)
#   # GTindex=qc$GTindex
#   Y=qc$Y
#   KI=qc$KI
#   CV=qc$CV
#   Z=qc$Z
#   GK=qc$GK
# }

my_taxa=as.character(KI[,1])
my_allKI=KI


print("The value of QC is")
# print(QC)
# rm(qc)
gc()
print("------------Examining data (QC) done-------------------------------------")
super_pass=FALSE
SUPER_myKI=NULL
SUPER_optimum_GD=NULL
if (!is.null(sangwich.top)) super_pass=TRUE
if(super_pass)
{
  print("-------------------start SUPER BREAD-----------------------------------")
#Create GK if not provided
#print(memory.size())
  if(is.null(GK))
  {
    nY=floor(nrow(Y)*.9)
    nG=ncol(GD)
    # snpsam=0
    if(nG>nY){snpsam=sample(1:nG,nY)}else{snpsam=1:nG}
    GK=GD[,snpsam]
    SNPVar=apply(as.matrix(GK), 2, stats::var)
	#print(snpsam)
    # if(snpsam==1)stop ("GAPIT says: SUPER_GS must putin GD and GM.")
    GK=GK[,SNPVar>0]
    GK=cbind(as.data.frame(GT),as.data.frame(GK)) #add taxa 
  }
  #print(head(CV))
  #myGD=cbind(as.data.frame(GT),as.data.frame(GD)) 

  # file.output.temp=file.output
  # file.output=FALSE
#  print(memory.size())
  GP=GAPIT.Bread(Y=Y,CV=CV,Z=Z,KI=KI,GK=GK,GD=cbind(as.data.frame(GT),as.data.frame(GD)),GM=GI,method=sangwich.top,LD=LD,file.output=FALSE,CV.Extragenetic=CV.Extragenetic)$GWAS
  # file.output=file.output.temp
#  print(memory.size())
  GK=NULL

  if(inclosure.to>nrow(Y))   ##########removed by Jiabo Wang ,unlimited number of inclosures
  {
    inclosure.to=nrow(Y)-1
    print("the number of choosed inclosure is more than number of individuals")
    print("Set the number of choosed incolosure max equal to individuals")
  }
  if(inclosure.from>inclosure.to)   ##########removed by Jiabo Wang ,unlimited number of inclosures
  {
    inclosure.from=inclosure.to
  }
  bin.level=seq(bin.from,bin.to,by=bin.by)
  inclosure=seq(inclosure.from,inclosure.to,by=inclosure.by)
#print(inclosure)
  e=1 #################################number of bins and inclosure
  count=0
  num_selection=length(bin.level)*length(inclosure)
  SUPER_selection=matrix(,num_selection,6)
  colnames(SUPER_selection)=c("bin","pseudo_QTNs","Max_pQTNs","REML","VA","VE")
#for (bin in bin.level){bin=bin.level[e]}
#for (inc in inclosure){inc=inclosure[e]}
  for (bin in bin.level)
  {
    for (inc in inclosure)
    {
      count=count+1
      mySpecify=GAPIT.Specify(GI=GI,GP=GP,bin.size=bin,inclosure.size=inc)
      SNP.QTN=mySpecify$index
      num_pseudo_QTN=length(mySpecify$CB)
      num_bins=mySpecify$num_bins
#print(paste("bin---",bin,"---inc---",inc,sep=""))      
      GK=GD[,SNP.QTN]
      SUPER_GD=GD[,SNP.QTN]
      SNPVar=apply(as.matrix(GK), 2, stats::var)
      GK=GK[,SNPVar>0]
      SUPER_GD=SUPER_GD[,SNPVar>0]
      GK=cbind(as.data.frame(GT),as.data.frame(GK)) #add taxa
      SUPER_GD=cbind(as.data.frame(GT),as.data.frame(SUPER_GD)) #add taxa
      myBurger=GAPIT.Burger(Y=Y,CV=CV,GK=GK)  #modifed by Jiabo Wang
      myREML=myBurger$REMLs
      myVG=myBurger$vg
      myVE=myBurger$ve
      SUPER_selection[count,1]=bin
      SUPER_selection[count,2]=num_pseudo_QTN
      SUPER_selection[count,3]=num_bins
      SUPER_selection[count,4]=myREML
      SUPER_selection[count,5]=myVG
      SUPER_selection[count,6]=myVE
  #print(SUPER_selection[count,])
      if(count==1)
      {
        GK.save=GK
        LL.save=myREML
        SUPER_optimum_GD=SUPER_GD     ########### get SUPER GD
      }else{
        if(myREML<LL.save)
        {
          GK.save=GK
          LL.save=myREML
          SUPER_optimum_GD=SUPER_GD     ########### get SUPER GD
        }
      }
      if (num_bins==num_pseudo_QTN) break
    }# inc end
  }# bin end
  SUPER_selection<-SUPER_selection[!is.na(SUPER_selection[,1]),]
  print(SUPER_selection)
  print("-----select optimum pseudo QTNs from all the bins-------")
  if(is.null(dim(SUPER_selection)))
  {
    optimum_SUPER=SUPER_selection
  }else{
    optimum_SUPER=SUPER_selection[which(as.numeric(SUPER_selection[,4])==min(as.numeric(SUPER_selection[,4]))),]
  }
  print(optimum_SUPER)
  ########################BUILD SUPER KINSHIP
  ##########################################################
  colnames(SUPER_optimum_GD)=c("taxa",colnames(SUPER_optimum_GD)[-1])
  SUPER_taxa=as.character(SUPER_optimum_GD[,1])
  SUPER_X=SUPER_optimum_GD[,-1]
  if(kinship.algorithm=="Loiselle")SUPER_myKI_test= GAPIT.kinship.loiselle(snps=t(as.matrix(.5*(SUPER_optimum_GD[,-1]))), method="additive", use="all")
 # if(kinship.algorithm=="VanRaden")SUPER_myKI_test= GAPIT.kinship.VanRaden(snps=as.matrix(SUPER_optimum_GD[,-1])) 
  if(kinship.algorithm=="Zhang")SUPER_myKI_test= GAPIT.kinship.Zhang(snps=as.matrix(SUPER_optimum_GD[,-1])) 
  if(kinship.algorithm!="Loiselle"&kinship.algorithm!="Zhang")SUPER_myKI_test= GAPIT.kinship.VanRaden(snps=as.matrix(SUPER_optimum_GD[,-1])) 

# SUPER_myKI_test=GAPIT.kinship.VanRaden(snps=as.matrix(SUPER_optimum_GD[,-1]))     #  build kinship
  colnames(SUPER_myKI_test)=SUPER_taxa
  SUPER_myKI=cbind(SUPER_taxa,as.data.frame(SUPER_myKI_test))
  print("select optimum number of marker effect in GD")
  Z=NULL
  if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & is.null(Z))
  {
  taxa=as.character(Y[,1])
  Z=as.data.frame(diag(1,nrow(Y)))
  Z=rbind(taxa,Z)
  taxa=c('Taxa',as.character(taxa))
  Z=cbind(taxa,Z)
  }
if(kinship.algorithm!="None" & kinship.algorithm!="SUPER" & !is.null(Z))
{
  if(nrow(Z)-1<nrow(Y)) Z=GAPIT.ZmatrixFormation(Z=Z,Y=Y)
}
print("QC2 is in process...")
GK=NULL
# CVI <- CV
# print("@@@")
# print(dim(Y))
# print(dim(SUPER_myKI))
# print(length(GT))
# if(QC)
# {
#   qc <- GAPIT.QC2(Y=Y,KI=SUPER_myKI, GT=GT,CV=CV,Z=Z,GK=GK)
#   # GTindex=qc$GTindex
#   Y=qc$Y
#   KI=qc$KI
#   CV=qc$CV
#   Z=qc$Z
#   GK=qc$GK
# }

# rm(qc)
gc()
}# super_pass end

nk=1000000000
if(!is.null(KI)) nk=min(nk,nrow(KI))
if(!is.null(GK)) nk=min(nk,nrow(GK))
if(!is.null(KI))
{
  if(group.to>nk) {
    #group.to=min(nrow(KI),length(GTindex)) #maximum of group is number of rows in KI
    group.to=nk #maximum of group is number of rows in KI
    #warning("The upper bound of groups is too high. It was set to the size of kinship!") 
    print("The upper bound of groups is too high. It was set to the size of kinship!") 
  }
	if(group.from>nk){ 
    group.from=nk
    #warning("The lower bound of groups is too high. It was set to the size of kinship!") 
    print("The lower bound of groups is too high. It was set to the size of kinship!") 
  } 
}

if(!is.null(CV)){
 	if(group.to<=ncol(CV)+1) {
	#The minimum of group is number of columns in CV
	  group.from=ncol(CV)+2
	  group.to=ncol(CV)+2
	  #warning("The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!")
	  print("The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!")
	}
}
# print(group.to)
# print(group.from)
  GROUP=seq(group.to,group.from,by=-group.by)#The reverse order is to make sure to include full model
if(missing("kinship.cluster")) kinship.cluster=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid")
if(missing("kinship.group")) kinship.group=c("Mean", "Max", "Min", "Median")
numSetting=length(GROUP)*length(kinship.cluster)*length(kinship.group)
ys=as.matrix(Y[2])
# print(dim(CV))
# print(dim(allCV))
X0=as.matrix(CV[,-1,drop=FALSE])
if(min(X0[,1])!=max(X0[,1])) X0 <- cbind(1, X0) #do not add overall mean if X0 has it already at first column
hold_Z=Z

 # library("EMMREML")
order_count=0
storage_reml=NULL
Compression=matrix(,numSetting,6)
colnames(Compression)=c("Type","Cluster","Group","REML","VA","VE")

for (group in GROUP)
{
  for (ca in kinship.cluster)
  {
  for (kt in kinship.group)
  {
  #if(group=1) group=2
#if(!optOnly) {print("Compressing and Genome screening..." )}
order_count=order_count+1
if(order_count==1)print("-------Mixed model with Kinship-----------------------------")
# if(group<ncol(GD)+1) group=2 # the emma function (emma.delta.REML.dLL.w.Z) does not allow K has dim less then CV. turn to GLM (group=1)
# print(Sys.time())
if(nrow(KI)==group)
{
  # print(head(cbind(as.data.frame(KI[,1]),1:nrow(KI))))
  GA=cbind(as.data.frame(KI[,1]),1:nrow(KI))
  colnames(GA)=c("X1","X2")
  KG=as.matrix(KI[,-1])
  bk <- GAPIT.Block(Z=hold_Z,GA=GA,KG=KG)

  }else{
cp <- GAPIT.Compress(KI=KI,kinship.cluster=ca,kinship.group=kt,GN=group,Timmer=Timmer,Memory=Memory)
# print(cp$GA)

bk <- GAPIT.Block(Z=hold_Z,GA=cp$GA,KG=cp$KG)
}
# print(Sys.time())

zc <- GAPIT.ZmatrixCompress(Z=hold_Z,GAU =bk$GA)
# print(Sys.time())
zrow=nrow(zc$Z)
zcol=ncol(zc$Z)-1
K = as.matrix(bk$KW)
# print(colnames(K))
#if (nrow(as.matrix(bk$KW))==1)
Z=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol)
if(is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
if(is.null(X0)) X0 <- matrix(1, ncol(ys), 1)
#handler of special Z and K
if(!is.null(Z)){ if(ncol(Z) == nrow(Z)) Z = NULL }
if(!is.null(K)) {if(length(K)<= 1) K = NULL}
X <-  X0 #covariate variables such as population structure
j=1
  if (is.null(Z)) Z=diag(x=1,nrow(K),ncol(K))
  if (group==1)   K=1
  XX=GAPIT.Licols(X0)
  X=XX$Xsub
  X.idx=XX$idx
  # print(dim(X))
  # print(dim(K))
  # print(dim(Z))
  # print(dim(Z))
   emma_test <- EMMREML::emmreml(as.numeric(ys), X=as.matrix(X), K=as.matrix(K), Z=Z,varbetahat=TRUE,varuhat=TRUE, PEVuhat=TRUE, test=TRUE)  
   # emma_test <- emmreml(as.numeric(ys), X=as.matrix(X), K=as.matrix(K), Z=Z,varbetahat=TRUE,varuhat=TRUE, PEVuhat=TRUE, test=TRUE)  

   print(paste(order_count, "of",numSetting,"--","Vg=",round(emma_test$Vu,4), "VE=",round(emma_test$Ve,4),"-2LL=",round(-2*emma_test$loglik,2), "  Clustering=",ca,"  Group number=", group ,"  Group kinship=",kt,sep = " "))
  emma_test_reml=-2*emma_test$loglik
  storage_reml=append(storage_reml,-2*emma_test$loglik)
Compression[order_count,1]=kt
Compression[order_count,2]=ca
Compression[order_count,3]=group
Compression[order_count,4]=emma_test_reml
Compression[order_count,5]=emma_test$Vu
Compression[order_count,6]=emma_test$Ve
  if(order_count==1){
   save_remle=emma_test_reml
   optimum_group=group
   optimum_Clustering=ca
   optimum_groupK=kt
   optimum_h2=emma_test$Vu/(emma_test$Vu+emma_test$Ve)
}else{
  if(emma_test_reml<save_remle){
   save_remle=emma_test_reml
   optimum_group=group
   optimum_Clustering=ca
   optimum_groupK=kt
   optimum_h2=emma_test$Vu/(emma_test$Vu+emma_test$Ve)

  }
}
}   # kt end

  } # ka end
  } # group end
  Compression=Compression[order(as.numeric(Compression[,4]),decreasing = FALSE),]
  Compression=matrix(Compression,ncol=6,byrow=F)
  colnames(Compression)=c("Type","Cluster","Group","REML","VA","VE")
  print(Compression)

  # write.csv(Compression,paste("GAPIT.",Compression,".csv",sep=""), row.names = FALSE,col.names = TRUE)

 if(optimum_group==1)  
{
optimum_group=2
}
#print(colnames(KI)[53:62])
if(nrow(Compression)>1)
{
cp <- GAPIT.Compress(KI=KI,kinship.cluster=optimum_Clustering,kinship.group=optimum_groupK,GN=optimum_group,Timmer=Timmer,Memory=Memory)
bk <- GAPIT.Block(Z=hold_Z,GA=cp$GA,KG=cp$KG)

zc <- GAPIT.ZmatrixCompress(Z=hold_Z,GAU =bk$GA)
zrow=nrow(zc$Z)
zcol=ncol(zc$Z)-1

K = as.matrix(bk$KW)
Z=matrix(as.numeric(as.matrix(zc$Z[,-1])),nrow=zrow,ncol=zcol)

if(is.null(dim(ys)) || ncol(ys) == 1)  ys <- matrix(ys, 1, length(ys))
if(is.null(X0)) X0 <- matrix(1, ncol(ys), 1)
  # X <-  X0 #covariate variables such as population structure
  XX=GAPIT.Licols(X0)
  X=XX$Xsub
  X.idx=XX$idx
  if (is.null(Z)) Z=diag(x=1,nrow(K),ncol(K))
   emma_REMLE <- EMMREML::emmreml(y=as.numeric(ys), X=as.matrix(X), K=as.matrix(K), Z=Z,varbetahat=TRUE,varuhat=TRUE, PEVuhat=TRUE, test=TRUE)  
  }else{
   emma_REMLE=emma_test
   print("gBLUP with only one time emma")
  } 
  my_allX=cbind(1,as.matrix(my_allCV[,-1]))
  XCV=my_allX[,X.idx,drop=FALSE]
  QTN.gs=QTN.gs-(ncol(my_allX)-length(X.idx))
  #print("!!!!")
  # print(table(X.idx))
  print(dim(XCV))
  print(QTN.gs)

#CV.Extragenetic specified
    if(ncol(XCV)>1&(ncol(XCV)-QTN.gs)!=1) XCVI=XCV[,c((2+CV.Extragenetic):(ncol(XCV)-QTN.gs)),drop=FALSE]
    XCVN=XCV[,c(1:(1+CV.Extragenetic)),drop=FALSE]
    if(QTN.gs!=0)XCVqtn=XCV[,c((ncol(XCV)-QTN.gs):ncol(XCV)),drop=FALSE]
    if(ncol(XCV)>1)beta.I=emma_REMLE$betahat[c((2+CV.Extragenetic):(ncol(XCV)-QTN.gs))]
    beta.N=emma_REMLE$betahat[c(1:(1+CV.Extragenetic))]
    if(QTN.gs!=0)beta.QTN=emma_REMLE$betahat[c((ncol(XCV)-QTN.gs):ncol(XCV))]
    # print(dim(XCVI))
    # print(length(beta.I))
    BLUE.N=XCVN%*%beta.N
    BLUE.QTN=rep(0,length(BLUE.N))    
    if(QTN.gs!=0)BLUE.QTN=XCVqtn%*%beta.QTN
    BLUE.I=rep(0,length(BLUE.N))
    if(ncol(XCV)>1&(ncol(XCV)-QTN.gs)!=1)BLUE.I=XCVI%*%beta.I
    #Interception only
   # print(dim(my_allX))
   # print(length(emma_REMLE$betahat))
   BLUE=cbind(BLUE.N,BLUE.I,BLUE.QTN)
   BLUE=data.frame(cbind(data.frame(my_allCV[,1]),data.frame(BLUE)))
   colnames(BLUE)=c("Taxa","BLUE.N","BLUE.I","QTNs")
   # print(dim(BLUE))
   # emma_BLUE=as.matrix(my_allX)%*%as.matrix(emma_REMLE$betahat)
   # emma_BLUE=as.data.frame(cbind(as.character(my_allCV[,1]),emma_BLUE))
   # colnames(emma_BLUE)=c("Taxa","emma_BLUE")
   gs <- GAPIT.GS(KW=bk$KW,KO=bk$KO,KWO=bk$KWO,GAU=bk$GAU,UW=cbind(emma_REMLE$uhat,emma_REMLE$PEVuhat))
   BB= merge(BLUE,gs$BLUP, by.x = "Taxa", by.y = "Taxa",all.x=T,sort=F)
   # print(dim(gs$BLUP))
   #print(dim(BB)) 
   BB[is.na(BB)]=0
   gBreedingValue=BB[,3]+BB[,4]+BB[,8]
   Prediction=BB[,2]+BB[,3]+BB[,4]+BB[,8]


   all_gs=cbind(BB,gBreedingValue,Prediction)
   colnames(all_gs)=c("Taxa","BLUE.N","BLUE.I","QTNs","Group","RefInf","ID","BLUP","PEV","gBreedingValue","Prediction")
   # colnames(all_gs)=c("Taxa","Group","RefInf","ID","BLUP","PEV","BLUE","Prediction","Pred_Heritable")
   # print(head(all_gs))
   if(GAPIT3.output)utils::write.csv(all_gs,paste("GAPIT.Association.Prediction_results.",model,".",name.of.trait,".csv",sep=""), row.names = FALSE)
  
  print("GAPIT SUPER GS completed successfully for multiple traits. Results are saved")
  return (list(GPS=BB,Pred=all_gs,Compression=Compression,kinship=my_allKI,SUPER_kinship=SUPER_myKI,SUPER_GD=SUPER_optimum_GD ,PC=my_allCV,Timmer=Timmer,Memory=Memory,GWAS=NULL,h2=optimum_h2 ))

}
