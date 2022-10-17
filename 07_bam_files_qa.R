# File: 07_bam_files_qa.R
# Auth: umar.niazi@kcl.ac.uk
# Desc: Quality checks on the bam files
# Date: 10/10/2022

## set variables and source libraries
source('header.R')
library(Rsamtools)

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

dfMeta = dfSample

#### set working directory to appropriate location with bam files
setwd('results/bams/')
csFiles = list.files('.', pattern = '*sortc.bam$')
## match the file names to sample ids in metadata table
f = gsub('_q30_.+', '.bam', csFiles)
table(dfMeta$fname %in% f)
i = match(dfMeta$fname, f)
identical(dfMeta$fname, f[i])
dfMeta$fname_sorted = csFiles[i]

## load the list of files with duplicates removed
csFiles = list.files('.', pattern = '*sortc_rd.bam$')
f = gsub('_q30_.+', '.bam', csFiles)
table(dfMeta$fname %in% f)
i = match(dfMeta$fname, f)
identical(dfMeta$fname, f[i])
dfMeta$fname_sorted_rd = csFiles[i]

dim(dfMeta)
csFiles = dfMeta$fname_sorted

## perform the analysis one sample at a time
lAllBams = lapply(csFiles, function(x){
  c = countBam(x)
  cat(paste(x, 'done', '\n'))
  return(c)
})
dfBam = do.call(rbind, lAllBams)

## samples with duplicates removed
csFiles = dfMeta$fname_sorted_rd
lAllBams = lapply(csFiles, function(x){
  c = countBam(x)
  cat(paste(x, 'done', '\n'))
  return(c)
})
df = do.call(rbind, lAllBams)
dfBam = rbind(dfBam, df)

setwd(gcswd)
n2 = paste0('results/', 'bamFileCounts.csv')
write.csv(dfBam, file=n2)

### load the results if this part was run on cluster
dfBam = read.csv(n2, header=T, stringsAsFactors = F, row.names=1)

i = grep('sortc.bam', dfBam$file)
iReadCount.pre = round(dfBam$records[i] / 1e6, 2)
names(iReadCount.pre) = gsub('_q30_.+', '', dfBam$file[i])
iReadCount.post = round(dfBam$records[-i] / 1e6, 2)
names(iReadCount.post) = gsub('_q30_.+', '', dfBam$file[-i])
identical(names(iReadCount.pre), names(iReadCount.post))
## check if or otherwise sort names in same order as meta data
cvMeta = gsub('.bam', '', dfMeta$fname)
identical(cvMeta, names(iReadCount.pre))
identical(cvMeta, names(iReadCount.post))
names(iReadCount.pre) = dfMeta$title
names(iReadCount.post) = dfMeta$title

### create the plots of interest
setwd(gcswd)
pdf(file='results/bam.qa.pdf')

mReadCount = rbind(iReadCount.post, iReadCount.pre)
rownames(mReadCount) = c('Post', 'Pre')

barplot(mReadCount, beside=T, las=2, main='No. of reads aligned', ylab = 'No. of Reads in Millions', cex.names =0.3,
        ylim=c(0, ceiling(max(mReadCount)+2)))
legend('topright', legend = c('Dup Rem', 'Orig'), fill=c('black', 'grey'))

dev.off(dev.cur())
