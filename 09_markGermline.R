# File: 09_markGermline.R
# Auth: umar.niazi@kcl.ac.uk
# Desc: mark germline and somatic mutations
# Date: 17/12/2022

## set variables and source libraries
source('header.R')

## connect to mysql database to get sample information
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'File')
dbListFields(db, 'Sample')
# get the query
g_did
q = paste0('select Sample.id as sid, Sample.group1, Sample.group2, Sample.group3, Sample.title, Sample.description, File.id as fid, File.name as fname from Sample, File
           where (Sample.idData = 54) AND (File.idSample = Sample.id AND File.type like "%bam%")')
dfSample = dbGetQuery(db, q)
nrow(dfSample)
dfSample
# close connection after getting data
dbDisconnect(db)

library(VariantAnnotation)

## function to load csv results file from shearwater output
oVRloadCSV = function(csFile, p.cut=0.01){
  df = read.csv(csFile, header=T, row.names = 1, stringsAsFactors = F)
  df = df[df$p.adj <= p.cut,]
  v = VRanges(seqnames = df$chr, ranges = IRanges(start=df$pos, width = 1), ref = df$ref, alt = df$mut)
  mcols(v) = df
  n = paste(df$chr, df$pos, df$ref, df$mut, sep='.')
  names(v) = n
  return(v)
}

#i = which(dfSample$group1 == 'A005')
n = paste0(dfSample$title, '.csv')
setwd('dataExternal/remote/results/')
## check if files exist
l = list.files('.', '*.csv$')
table(n %in% l)
## load all the data
l = lapply(n, oVRloadCSV)
oVRLSample = VRangesList(l)
names(oVRLSample) = dfSample$title
rm(l)

cvSamples = dfSample$group1

## identify common variants across samples, to be used for Germline filtering
cvCommonSNPs = unique(do.call(c, lapply(oVRLSample, names)))
mCommonSNPs = matrix(NA, nrow=length(cvCommonSNPs), ncol=length(oVRLSample))
for (i in 1:ncol(mCommonSNPs)){
  mCommonSNPs[,i] = cvCommonSNPs %in% names(oVRLSample[[i]])
}
rownames(mCommonSNPs) = cvCommonSNPs
colnames(mCommonSNPs) = names(oVRLSample)


## process each sample group based on sample id
oVRoverlaps = function(oVR, iCounts){
  m = match(names(oVR), names(iCounts))
  stopifnot(identical(names(oVR), names(iCounts)[m]))
  oVR$overlaps = iCounts[m]
  return(oVR)
}

cvSamples.u = unique(cvSamples)

for (s in seq_along(cvSamples.u)){
  w = which(cvSamples == cvSamples.u[s])
  if (length(w) <= 1) next;
  i = rowSums(mCommonSNPs[,w])
  l = lapply(oVRLSample[w], oVRoverlaps, i)
  oVRLSample[w] = VRangesList(l)
}

lResults = list(vr=oVRLSample, meta=dfSample)
save(lResults, file='lResults_with_VRanges_overlaps.rds')

## write csv files
n = paste0(names(lResults$vr), '.csv')
setwd('results/csv/')
sapply(seq_along(n), function(i){
  d = mcols(lResults$vr[[i]])
  write.csv(d, file=n[i])
})

## figures for overlaps by subject
oVRall = unlist(oVRLSample)
dfOverlaps = data.frame(Sample=oVRall$sampleID, Overlaps=oVRall$overlaps) 
dfOverlaps$Subject = factor(gsub('_\\d+$', '', as.character(dfOverlaps$Sample)))
xtabs(~ Subject + Overlaps, data=dfOverlaps)

library(lattice)

histogram(~ Overlaps | Subject, data=dfOverlaps[dfOverlaps$Subject != 'A004',], scales=list(relation='free'))

# 
# m = match(names(v), names(i))
# # sanity check
# identical(names(v), names(i)[m])
# v$overlaps = i[m]

