## pava is the pool adjacent violator algorithm to perform isotonic transformation
pava <- function (x, wt = rep(1, length(x))) 
{
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}


# Simulation of clinical trial results
teqrOC<-function(sim,firstdose=2,probt,cohortSize=3,MaxNoCohorts=30,MTDss=12,pTarget,eq1,eq2,tootoxic){
results.old<-c(NA,NA,NA,NA,NA,NA)
for (j in 1:sim){
#print(j)
tox.old<-0
cumtox.old<-0
newdose<-firstdose
stopdose<-length(probt)
#print(newdose)
doselevel.old<-0
for (i in 1:MaxNoCohorts){
tox=rbinom(1,cohortSize, prob=probt[newdose])
tox.old<-rbind(tox.old,tox)
doselevel.old<-rbind(doselevel.old,newdose)
output<-data.frame(cbind(doselevel.old,tox.old), row.names=NULL)
colnames(output)<-c('doselevel','tox')
output<-output[c(-1),]
#print(paste('i=', i))
#print(output)
currentlevel<-output[output$doselevel==newdose,]
cumtox<-sum(currentlevel$tox)/(cohortSize*length(currentlevel$tox))
cumtox.old<-rbind(cumtox.old,cumtox)
#print(paste('newdose', newdose))
#print(paste('cumtox',cumtox))
upperlimit<-pTarget+eq2
lowerlimit<-pTarget-eq1
stopdosem1<-stopdose-1
#print(paste('stopdosem1', stopdosem1))
if (newdose<stopdose & newdose>1){
if (cumtox>lowerlimit) doselevel<-newdose
if (cumtox>upperlimit) doselevel<-newdose-1
if (cumtox>=tootoxic)  stopdose<-newdose-1
if (cumtox<lowerlimit)  doselevel<-newdose+1 
}
if (newdose==1){
if (cumtox>upperlimit) break
if (cumtox<upperlimit)  doselevel<-newdose 
if (cumtox<lowerlimit)  doselevel<-newdose+1 
}
if (newdose==stopdose){
if (cumtox>upperlimit & newdose>1) doselevel<-newdose-1
if (cumtox>=tootoxic)  stopdose<-newdose-1
if (cumtox<upperlimit)   doselevel<-newdose 
}
newdose<-doselevel
if (cohortSize*length(currentlevel$tox)>=MTDss & cumtox<tootoxic) break}
cumtox.m1<-cumtox.old[-c(1)]
simNo<-rep(j, length(output$tox))
results<-cbind(simNo,output, cumtox.m1)
results.old<-rbind(results.old, results)
}
simresults<-data.frame(results.old[-c(1),])
simData<-(list(simresults=simresults, cohortSize=cohortSize,probt=probt, MTDss=MTDss, pTarget=pTarget, lowerlimit=lowerlimit,upperlimit=upperlimit,tootoxic=tootoxic, sim=sim))
DLdata<-teqr.DLdata(simData=simData)
MTDdata<-teqr.MTDdata(simData=simData, DLdata=DLdata)
OperChar<-teqr.OperChar(simData=simData, DLdata=DLdata,MTDdata=MTDdata)
OperChar
}

#orders the clinical trial simulation results by doselevel
teqr.DLdata<-function(simData=simData){
sim<-simData$sim
SDR<-simData$simresults
re.old<-rep(NA,7)
for (i in 1:sim){
#print(i)
output<-SDR[SDR$simNo==i, ]
mindose<-min(output$doselevel)
maxdose<-max(output$doselevel)
for (j in mindose:maxdose){
outputd<-output[output$doselevel==j, ]
simNo<-outputd$simNo
doselevel<-outputd$doselevel
stox<-sum(outputd$tox)
dllength<-simData$cohortSize*length(outputd$tox)
if (dllength==0) next
toxl<-binom.test(stox,dllength)$conf[1]
toxu<-binom.test(stox,dllength)$conf[2]
toxest<-binom.test(stox,dllength)$estimate
re1<-cbind(simNo,doselevel,stox,dllength,toxl,toxu,toxest)
re<-re1[1,]
re.old<-rbind(re.old,re)
}
}
results<-data.frame(re.old[-c(1),], row.names=NULL)
names(results)<-c('simNo','doselevel','stox','dllength','toxl','toxu','toxest')
return(results=results)
}

#MTDdataset
teqr.MTDdata<-function(simData=simData,DLdata=DLdata){
doselevel=rep(NA,simData$sim)
for (i in 1:simData$sim){
output<-DLdata[DLdata$simNo==i, ]
iv<-rep(i,length(output$doselevel))
dat<-data.frame(cbind(c(iv),output$doselevel,round(c(pava(output$stox/output$dllength)),2)))
names(dat)<-c('i','dl','tox')
dat$tt<- ifelse(dat$tox<simData$tootoxic,1,0)
#print(dat)
if (sum(dat$tt)>0) {
newdat<-dat[dat$tt>0,] 
newdat$diff<-abs(newdat$tox-simData$pTarget)
newdat$mindiff<-min(newdat$diff)
mindiffdat<-newdat[newdat$diff==newdat$mindiff,]
doselevel[i]<-max(mindiffdat$dl)
}
if (sum(dat$tt)==0) doselevel[i]<-100
}
dose0<-DLdata$simNo[DLdata$doselevel==1 & DLdata$toxest>simData$upperlimit]
simNo<-seq(1:simData$sim)
maxpava<-data.frame(cbind(cbind(simNo),cbind(doselevel)))
for(i in dose0) maxpava$doselevel[i]<-0
#maxpava
Remax<-merge(maxpava,DLdata)
o<-order(Remax$simNo)
Remax.o<-Remax[o,]
return(Remax.o=Remax.o)
}

teqr.OperChar<-function(simData=simData,DLdata=DLdata,MTDdata=MTDdata){
totalN<-rep(NA,simData$sim)
for (i in 1:simData$sim){ 
totalN[i]<-sum(DLdata$dllength[DLdata$simNo==i])
}
MedN<-median(totalN)out<-table(simData$simresults$sim, simData$simresults$doselevel)
NoPatients<-(simData$cohortSize*margin.table(out,2))/simData$sim
sstox<-rep(NA,simData$sim)
for (i in 1:simData$sim) sstox[i]<-sum(DLdata$stox[DLdata$simNo==i])
sdllength<-rep(NA,simData$sim)
for (i in 1:simData$sim) sdllength[i]<-sum(DLdata$dllength[DLdata$simNo==i])
DLTrate.dat<-data.frame(cbind(sstox,sdllength))
DLTrate<-sstox/sdllength
MeanDLTrate<-mean(DLTrate)
MeanToxRate<-mean(MTDdata$toxest)
MeanCIlength<-mean(MTDdata$toxu-MTDdata$toxl)
#percent of simulation achieving a defined samplesize
stopNo.m1=simData$MTDss-1
stopNo<-ifelse(MTDdata$dllength>stopNo.m1,1,0)
PropObd<-prop.table(table(stopNo))[2]
NoMTD<-simData$sim-length(MTDdata$toxest)
Retab<-table(MTDdata$doselevel)
NoTrialsMTD<-prop.table(Retab)
NoTrialsMTD1<-prop.table(cbind(t(Retab),NoMTD))
oc<-list(sim=simData$sim,MedN=MedN,NoPatients=NoPatients,MeanDLTrate=MeanDLTrate,NoTrialsMTD=NoTrialsMTD,NoTrialsMTD1=NoTrialsMTD1,MeanToxRate=MeanToxRate,MeanCIlength=MeanCIlength,PropObd=PropObd, NoMTD=NoMTD,simData=simData, DLdata=DLdata, MTDdata=MTDdata)class(oc) <- "teqrOC"
    oc
}

#dose escalation guidelines
teqrDG<-function(TotalN,pTarget,eq1,eq2,tootoxic){
Nplus1<-TotalN+1
upperlimit<-pTarget+eq2
lowerlimit<-pTarget-eq1
prob<-matrix(NA, ncol=TotalN, nrow=Nplus1, byrow=TRUE)
for (i in 1:Nplus1){
prob[i,]<-((i-1)/seq(1:TotalN))
}
row<-seq(1:Nplus1)-1
col<-seq(1:TotalN)
row.names(prob)<-row
colnames(prob)<-col
prob
prob2<-matrix(NA, ncol=TotalN, nrow=Nplus1, byrow=TRUE)
for (i in 1:Nplus1){
for (j in 1:TotalN) {
if (prob[i,j]<lowerlimit) prob2[i,j]<-'E'
if (prob[i,j]>= lowerlimit & prob[i,j] <= upperlimit) prob2[i,j]<-'S'
if (prob[i,j]>upperlimit) prob2[i,j]<-'D'
if (prob[i,j]>=tootoxic) prob2[i,j]<-'DU'
if (prob[i,j]>1.00)  prob2[i,j]<-'  '
}
} 
row<-seq(1:Nplus1)-1
col<-seq(1:TotalN)
row.names(prob2)<-row
colnames(prob2)<-col
prob2
dg<-(list(probTable=round(prob,3), DoseGuideTable=prob2))
class(dg) <- "teqrDG"
   dg
}





#print function
print.teqrOC<-function (x,...) 
{
    cat("is an object of class teqrOC.   \n")
    cat("                              \n")  
    cat("Ave. no of Patients studied at each dose level\n")
    print(round(x$NoPatients,2)) 
    cat("                              \n")  
    cat("Rate dose level is chosed as the MTD \n")
    print(round(x$NoTrialsMTD1,2))
    cat("                              \n") 
    cat("Median Study Sample Size:", x$MedN,"\n")
    cat("Mean Study DLT Rate:", round(x$MeanDLTrate,2),"\n")
    cat("Mean Toxicity Rate at the MTD:", round(x$MeanToxRate,2),"\n")
    cat("Mean 95% binomial confidence interval length at the MTD:", round(x$MeanCIlength,2),"\n") 
    cat("Proportion of trials with MTD dose level sample size at or above the desired number:", round(x$PropObd,2),"\n")
    cat("No of simulated trials that do not determine an MTD:", x$NoMTD,"\n")  
    cat("No of simulated trials:", x$sim,"\n")  
    cat("The following simulation data sets are also contained within the TEQR object.\n")
    cat("The simulation level data: simData \n")
    cat("The dose level data: DLdata \n")
    cat("The MTD level data: MTDdata \n")
    invisible(NULL)
}

#print function
print.teqrDG<-function (x,...) 
{
    cat("is an object of class teqrDG which contains the dose escalation/de-escalation\n") 
    cat("guidelines table.Note the rows represent number of subjects that have experienced\n")
    cat("a DLT and the columns represent the number of subjects on the current dose level.\n")
    cat("The letter codes are the guidelines and the letters are defined as follows.\n")
    cat("E-escalate, S-Stay, D- De-escalate, and DU- De-escalate\n") 
    cat("and do not return to this dose. \n")  
    print(x$DoseGuideTable) 
    cat("                              \n")      
    cat("To print the dosing escalation/de-escalation guidelines\n")
    cat("or the underlying toxicity probabilities as objects\n") 
    cat("type x$DoseGuideTable and x$probTable, respectively.\n")
    cat("In this example the object name is output,\n") 
    cat("so the user types output$DoseGuideTable or output$probTable\n")  
    invisible(NULL)
}


