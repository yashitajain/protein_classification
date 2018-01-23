#Reading training and testing protein sequences in fasta format

library(mlr)
library(ggplot2)
library(Peptides)
library(Biostrings)
library(e1071)
set.seed(40)

fasta_data=function(fas){
  s = readDNAStringSet(fas)
  seq_name = names(s)
  sequence = paste(s)
  ss=data.frame(seq_name, sequence)
  return(ss)
}
traindata=fasta_data('target.faa')
testdata=fasta_data('sg_putative_10000.faa')

#Building feature variables by using different physiochemical properties like aaComp, hydrophobicity, molecular weight etc. using Peptide module

physiochem=function(ss,x){
  aaCom=NULL
  CP=NULL
  fasgai=NULL
  ai=NULL
  boman=NULL
  kid=NULL
  hydro=NULL
  insta=NULL
  mow=NULL
  pi=NULL
  blosum=NULL
  len=NULL
  ms=NULL
  prot=NULL
  st=NULL
  tscale=NULL
  zscale=NULL
  for (i in 1:x){
    aaCom=rbind(aaCom,t(data.frame(aaComp(ss[i,2])))[2,])
    CP=rbind(CP,t(data.frame(crucianiProperties(ss[i,2]))))
    fasgai=rbind(fasgai,t(data.frame(fasgaiVectors(ss[i,2]))))
    ai=rbind(ai,t(data.frame(aIndex(ss[i,2]))))
    boman=rbind(boman,t(data.frame(boman(ss[i,2]))))
    kid=rbind(kid,t(data.frame(kideraFactors(ss[i,2]))))
    hydro=rbind(hydro,t(data.frame(hydrophobicity(ss[i,2],scale = "BlackMould"))))
    insta=rbind(insta,t(data.frame(instaIndex(ss[i,2]))))
    mow=rbind(mow,t(data.frame(mw(ss[i,2],monoisotopic = F))))
    pi=rbind(pi,t(data.frame(pI(ss[i,2],pKscale = "EMBOSS"))))
    blosum=rbind(blosum,t(data.frame(blosumIndices(ss[i,2]))))
    len=rbind(len,t(data.frame(lengthpep(ss[i,2]))))
    ms=rbind(ms,t(data.frame(mswhimScores(ss[i,2]))))
    prot=rbind(prot,t(data.frame(protFP(ss[i,2]))))
    st=rbind(st,t(data.frame(stScales(ss[i,2]))))
    tscale=rbind(tscale,t(data.frame(tScales(ss[i,2]))))
    zscale=rbind(zscale,t(data.frame(zScales(ss[i,2]))))
  }
  result=cbind(ss,aaCom,fasgai,ai,boman,kid,hydro,insta,mow,pi,blosum,len,ms,prot,st,tscale,zscale)
  return(result)
  
}
train=physiochem(traindata,264)
test=physiochem(testdata,10000)


label=factor(c(rep('1',44),rep('0',220)))
names(label)='label'
train=cbind(train,label)

#Making training and test data set for the model.

x=subset(train,select=c(-sequence,-seq_name,-label))
test_x=subset(test,select=c(-sequence,-seq_name))

data.train=cbind(x,label)
data.test=test_x

#Applying Decision Tree model using MLR module

set.seed(26)

task = makeClassifTask(data = data.train, target = "label")

table(getTaskTargets(task))

lrn = makeLearner("classif.rpart", predict.type = "prob")
mod = train(lrn, task)


p1=predict(mod, newdata = data.test)

testdata=cbind(testdata,p1)
test1=testdata[which(testdata$response=='1'),]
w=test1[,1]
write.csv(w,file="ct_model.csv")

#Validating the Decision Tree model using Cross-Validation. 

#validation- CV 10-fold training set
rdesc = makeResampleDesc("CV", iters = 10)
ms = list(auc, mmce)
mr = benchmark(lrn, tasks = task, resampling = rdesc, measures = ms, show.info = FALSE)
df = generateThreshVsPerfData(mr, measures = list(fpr, tpr, mmce),aggregate = F)
#performance
plotThreshVsPerf(df) +
  theme(strip.text.x = element_text(size = 12))+geom_tile(color='darkblue')
#roc_curve
plotROCCurves(df)+theme(strip.text.x = element_text(size = 12))+geom_path(color='darkblue')


