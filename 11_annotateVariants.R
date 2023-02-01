# File: 11_annotateVariants.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: Annotate the identified somatic mutations using shearwater
# Date: 9/1/2023

source('header.R')
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(VariantAnnotation)

lResults = f_LoadObject(file.choose())

oVRall = unlist(lResults$vr)
## some filters to drop variants and samples
table(oVRall$sampleID, oVRall$overlaps)
oVRall = oVRall[oVRall$sampleID != 'A004_1']
oVRall = oVRall[oVRall$overlaps < 9]
sampleNames(oVRall) = as.character(oVRall$sampleID)
## find metadata matches for sample id
dfMeta = lResults$meta
i = match(oVRall$sampleID, dfMeta$title)
identical(dfMeta$title[i], oVRall$sampleID)
oVRall$patient = dfMeta$group1[i]
oVRall$skin_location = dfMeta$group2[i]

## use the snplocs database
oDBSnp = SNPlocs.Hsapiens.dbSNP144.GRCh37
## change the style to match snplocs package
seqlevelsStyle(oVRall); seqlevelsStyle(oDBSnp)
## do some acrobatics with seqlevels to perform search
## IF this will give error, otherwise ignore next few lines
oGPsnps = snpsByOverlaps(oDBSnp, oVRall)
  
# seqlevels(oGRsnps) = seqlevels(oDBSnp)
# seqinfo(oGRsnps) = seqinfo(oDBSnp) 
# oGPsnps = snpsByOverlaps(oDBSnp, oGRsnps)
# 
# manually check a few rs numbers if any hits at
# https://www.ncbi.nlm.nih.gov/snp/?term=rs3745946
#IUPAC_CODE_MAP
oGPsnps

## update the SNPs VRanges object
f = findOverlaps(oVRall, oGPsnps)
oVRall$RefSNP_id = NA_character_
oVRall$alleles_as_ambig = NA_character_
oVRall$RefSNP_id[queryHits(f)] = oGPsnps$RefSNP_id[subjectHits(f)]
oVRall$alleles_as_ambig[queryHits(f)] = oGPsnps$alleles_as_ambig[subjectHits(f)]

## some stats
d = data.frame(DPSnps=!is.na(oVRall$RefSNP_id), Patient=oVRall$patient)
xtabs( ~ DPSnps + Patient, data=d)

## location of the variants 
oTxdb = TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevelsStyle(oTxdb)
seqlevelsStyle(oVRall)
seqlevelsStyle(oTxdb) = seqlevelsStyle(oVRall)
genome(oVRall)
genome(oTxdb)
genome(oVRall) = genome(oTxdb)
oGRlocations = locateVariants(oVRall, oTxdb, AllVariants())
length(oGRlocations)
## as comparison is unstranded, remove those matches that do not
## have a gene id. As the variants are called relative to + strand,
## and they may match a gene which is on + or - strand in this case.
## If there is no gene id associated with the match, then remove those matches.
oGRlocations = oGRlocations[oGRlocations$GENEID %in% oVRall$ENTREZID]
length(oGRlocations)

library(org.Hs.eg.db)
## summarize each snp by location
dfSymbol = AnnotationDbi::select(org.Hs.eg.db, keys = oGRlocations$GENEID, columns = 'SYMBOL',
                                   keytype = 'ENTREZID')
table(dfSymbol$SYMBOL %in% oVRall$SYMBOL)
identical(oGRlocations$GENEID, dfSymbol$ENTREZID)
oGRlocations$SYMBOL = dfSymbol$SYMBOL  
table(oGRlocations$SYMBOL, oGRlocations$LOCATION)

## create a new VRanges SNPs object and copy
## the annotations from locations object
length(unique(oGRlocations$QUERYID)); length(oVRall);

## use the query id mapping variable from locations
oVRall.locations = oVRall[oGRlocations$QUERYID]
length(oVRall.locations)
# sanity check
table(overlapsAny(oVRall, oVRall.locations))
length(oVRall.locations); length(oGRlocations);
oVRall.locations$locationStrand = strand(oGRlocations)
## merge the metadata
colnames(mcols(oGRlocations))
mcols(oVRall.locations) = cbind(mcols(oVRall.locations), mcols(oGRlocations)[,-10])

lResults$gr_locations = oGRlocations
lResults$vr_with_locations = oVRall.locations
save(lResults, file='results/lResults_with_VRanges_overlaps_locations.rds')

## get the unique positions i.e. drop the extra transcripts
length(oVRall.locations)
i = unique(oVRall.locations$QUERYID)
i = match(i, oVRall.locations$QUERYID) 
oVRall.locations.slim = oVRall.locations[i]
# sanity check
table(names(oVRall.locations.slim) %in% names(oVRall.locations))

d = data.frame(Location=oVRall.locations.slim$LOCATION, 
               Patient=oVRall.locations.slim$patient, oVRall.locations.slim$SYMBOL)
xtabs(~ Location + Patient, data=d)

write.csv(mcols(oVRall.locations.slim), file='results/SNPsWithLocationsAllSamples.csv')
### extract the coding locations
oGRcoding = oGRlocations[oGRlocations$LOCATION == 'coding']
length(oGRcoding)

oVRall.coding = subsetByOverlaps(oVRall, oGRcoding)
v = DNAStringSet(oVRall.coding$mut)
library(BSgenome.Hsapiens.UCSC.hg19)
oBSGenome = BSgenome.Hsapiens.UCSC.hg19
seqlevelsStyle(oVRall.coding)
seqlevelsStyle(oBSGenome) = 'NCBI'

oVRall.coding.prot = predictCoding(oVRall.coding, oTxdb, oBSGenome, v, ignore.strand=T)
table(oVRall.coding.prot$CONSEQUENCE)

lResults$vr_coding_proteins = oVRall.coding.prot
save(lResults, file='results/lResults_with_VRanges_overlaps_locations_coding_consequences.rds')

f = oVRall.coding.prot$CONSEQUENCE == 'nonsynonymous'
table(f)
d = data.frame(Patient=oVRall.coding.prot$patient[f], 
               Gene=oVRall.coding.prot$SYMBOL[f])
table(d$Patient)
df = as.data.frame(table(d))
write.csv(df, file='results/nonsynonymous_frequency.csv')

write.csv(mcols(oVRall.coding.prot[f]), file='results/nonSynonymous_mutations.csv')
## get the unique positions i.e. drop the extra transcripts
length(oVRall.coding)
i = unique(oVRall.coding.prot$QUERYID)
i = match(i, oVRall.coding.prot$QUERYID) 
oVRall.coding.prot.slim = oVRall.coding.prot[i]
# sanity check
identical(names(oVRall.coding.prot.slim), names(oVRall.coding))
table(oVRall.coding.prot.slim$CONSEQUENCE)

################################################
### polyphen databsae
library('PolyPhen.Hsapiens.dbSNP131')
oDBpoly = PolyPhen.Hsapiens.dbSNP131
columns(oDBpoly)

iSel = which(oVRall.coding.prot.slim$CONSEQUENCE == 'nonsynonymous' & !is.na(oVRall.coding.prot.slim$RefSNP_id))
length(iSel)

dfPolyphen = select(oDBpoly, keys = oVRall.coding.prot.slim$RefSNP_id[iSel],
                    columns=c('PREDICTION', 'PPH2PROB'), keytype = 'RSID')

dfPolyphen = na.omit(dfPolyphen)

iSel = which(oVRall.coding.prot.slim$RefSNP_id %in% dfPolyphen$RSID)
length(iSel)
oVRall.polyphen = oVRall.coding.prot.slim[iSel] 
cvRsid = dfPolyphen$RSID[dfPolyphen$PREDICTION == 'probably damaging']
oVRall.polyphen[oVRall.polyphen$RefSNP_id %in% cvRsid]
