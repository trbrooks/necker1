## Sep. 10 2015
## Un-nests psychopy data so it can be analyzed in R

getwd()
setwd("/Users/trbrooks/Documents/necker/data directory/necker_round1/data")
dir()->file.names
grep("csv", file.names)->temp
file.names[temp]->pydata1

# ** ALL filenames in the pydata vector should start with a subject number!!**
pydata1

#remove un-numbered cases
pydata1<-pydata1[grep("^[0-9].*", pydata1)]

#remove deleted cases
exclude<-c("001_","002_","004_")
exclude<-grep(paste(exclude,collapse="|"),pydata1,value=TRUE)
pydata<-pydata1[! pydata1 %in% exclude]

pydata

doubler <- function(x) {
  strsplit(as.character((ts1[,x])[1]),",")[[1]]->temp
  as.numeric(gsub("\\[|\\]", temp, replacement=""))->out
  options(scipen=999)
  return (out)
}
facer <- function(x) {
  strsplit(as.character((ts1[,x])[1]), ",")[[1]]->temp
  (gsub("[[:punct:]]", temp, replacement="") == " True")->temp
  temp[1:setlength]->is_up
  is_up[2]->is_up[1]
  is_up<-c(is_up[1], is_up)
  return (is_up)
}
charact.Er <- function(x) {
  rep(as.character(ts1[[x]]),setlength)->out
  return (out)
}
lengthset.Er <- function(x) {
  if (length(x)>length(cstate)){
  out <- x[1:length(cstate)]
  }
  else if (length(x)<length(cstate)){
    length(cstate)-length(x)->temp
    out <- c(x,rep(x[length(x)],temp))
  }
  else {
    out <- x
  }
  return (out)
}


#create the dataframe from the python data
#i<-6
for (i in 1:length(pydata)){

  options(stringsAsFactors=FALSE)
  read.csv(pydata[i], header=TRUE)->ts1
  

  doubler(1)->cstate
  doubler(2)->ctime
  doubler(3)->tstamp
  length(cstate)->setlength
  facer(4)->is_up
  lengthset.Er(charact.Er(5))->iscong
  lengthset.Er(charact.Er(6))->istest
  
  #time (as in time of day) and date
  rep(strsplit(ts1[,7], "_")[[1]][4],length(cstate))->timeofday
  rep(as.numeric(as.numeric(ts1[[8]])),setlength)->frame.rate
  
  ##date
  c(strsplit(ts1[,7], "_")[[1]][1:3])->temp
  rep(paste(temp, collapse = '_'), length(cstate))->date
  
  ##
  charact.Er(9)->trial
  lengthset.Er(charact.Er(10))->session
  lengthset.Er(charact.Er(11))->participant
  
  ifelse(iscong[1]==1,!is_up->temp,is_up->temp)
  lengthset.Er(temp)->temp
  (as.numeric(temp) == cstate)->temp
  as.numeric(temp)->mismatch
  #as.numeric(c(temp, rep(NA, length(cstate) - length(temp))))->mismatch
  
  change<-0
  for (j in 2:length(cstate)) {
    ifelse(cstate[j]==cstate[j-1], change[j]<-0, change[j]<-1)
  }
  change[1]<-change[2]
  for (j in length(change):length(cstate)) {
    change[j]<-0
  }
  
  lengthset.Er(i)->trial
  
  x <- data.frame(ctime, cstate, date, timeofday, 
                  session, participant, mismatch, 
                  iscong, istest, change, trial,
                  deparse.level=1, check.rows=TRUE)
  
  
  ifelse(i>1,rbind(df1,x)->df1,x->df1)
  print(i)
}

#Did the loop run to the end?
(length(pydata)==i)
