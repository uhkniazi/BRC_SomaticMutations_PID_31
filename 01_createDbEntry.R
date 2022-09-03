# File: 01_createDbEntry.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: list the samples and create appropriate db entries
# Date: 3/9/2022


## set variables and source libraries
source('header.R')

## load the metadata file
dfMeta = read.csv('dataExternal/metadata.csv', header=T, stringsAsFactors = F)
str(dfMeta)

## some acrobatics to get data in right format
s = strsplit(dfMeta$Patient.and.location, split = ' ')
p = sapply(s, function(x) return(x[1]))
l = sapply(s, function(x) return(x[2]))
table(p); table(l)
table(p, l)

dfMeta = data.frame(f=dfMeta$FileName, p, l, s=dfMeta$Size, stringsAsFactors = F)

## connect to mysql database 
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# sample and file table
dbGetQuery(db, paste('describe Sample;'))
cSampleCol = dbGetQuery(db, paste('describe Sample;'))$Field[-1]

dbGetQuery(db, paste('describe File;'))
cFileCol = dbGetQuery(db, paste('describe File;'))$Field[-1]

## some acrobatics to create unique sample title
# combining patient id and number of replicates for the patient
# its not necessary but if you wish to make things complicated
t = split.data.frame(dfMeta, factor(p))
t = lapply(t, function(x) { 
  x$t = paste0(unique(x$p), '_', 1:length(x$p))
  return(x)})
names(t) = NULL
t = do.call(rbind, t)
dim(t)
i = t[as.character(1:73), ]
identical(dfMeta, i[,-5])
dfMeta = i
## create the entry for samples
cSampleCol

dfSamples = data.frame(idProject=g_pid, idData=g_did, title=dfMeta$t,
                       location='/scratch/prj/lynch_ls_pscc_seq/',
                       description= paste('group1 is Patient',
                                          'group2 is skin sample location',
                                          sep=';'),
                       group1 = dfMeta$p, group2= dfMeta$l)
# write this data to the database
rownames(dfSamples) = NULL

### NOTE: Do not execute this anymore as entry created
# write this table to database
# dbWriteTable(db, name='Sample', value=dfSamples, append=T, row.names=F)
# get this table again from database with ids added
g_did
dfSamples = dbGetQuery(db, paste('select * from Sample where Sample.idData = 54;'))

# create entries for these files in the database
dbListTables(db)
cn = dbListFields(db, 'File')[-1]
cn

identical(dfSamples$title, dfMeta$t)
dfSamples$f = dfMeta$f

# get the file names
dfFiles = data.frame(name=dfSamples$f, type='bam', idSample=dfSamples$id, group1='Aligned by sequencing centre', stringsAsFactors = F)
identical(dfFiles$name, dfMeta$f)
rownames(dfFiles) = NULL

# write this table to database
## note: do not execute as it is already done
# dbWriteTable(db, name='File', value=dfFiles, append=T, row.names=F)

dbDisconnect(db)
