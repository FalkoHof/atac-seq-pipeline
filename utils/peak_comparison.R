library("GenomicRanges")
library("GenometriCorr")
library('ggbio')
library('VennDiagram')
library('limma')
library('seqLogo')

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

#path of track 1
track1.path = args[1]

#path of track 2
track2.path = args[2]

tss.path = args[3]


fseq.colnames <-c('chromosome','start','end','peak_id','score')

tair10.golden.path.length <- c(Ath_chr1 = 30427671, 
                               Ath_chr2 = 19698289, 
                               Ath_chr3 = 23459830, 
                               Ath_chr4 = 18585056, 
                               Ath_chr5 = 26975502)

track1.df <- read.table(track1.path, 
                       header = FALSE, 
                       sep = "\t", 
                       quote = "", 
                       col.names=fseq.colnames, 
                       stringsAsFactors=F)

track2.df <- read.table(track2.path, 
                        header = FALSE, 
                        sep = "\t", 
                        quote = "", 
                        col.names=fseq.colnames, 
                        stringsAsFactors=F)

tss.df <- read.table(tss.path, 
                        header = FALSE, 
                        sep = "\t", 
                        quote = "", 
                        col.names=fseq.colnames, 
                        stringsAsFactors=F)



track1.granges <- GRanges(Rle(track1.df$chromosome),
                          IRanges(start = track1.df$start, end = track1.df$end),
                          peakId = track1.df$peak_id, 
                          score = track1.df$score, 
                          seqlengths=tair.golden.path.length)

track2.granges <- GRanges(Rle(track2.df$chromosome),
                          IRanges(start = track2.df$start, end = track2.df$end), 
                          peakId = track2.df$peak_id, 
                          score = track2.df$score, 
                          seqlengths=tair.golden.path.length)

tss.granges <- GRanges(Rle(tss.df$Ath_chromosome),
                       IRanges(start = tss.df$Start, end = tss.df$End),
                       ModeLocation = tss.df$ModeLocation, 
                       Shape = tss.df$Shape, 
                       TranscriptID = tss.df$TranscriptID,
                       GeneName=tss.df$GeneName )


sum(countOverlaps(track1.granges,track2.granges))
sum(countOverlaps(track1.granges,tss.granges))
sum(countOverlaps(track2.granges,tss.granges))


plot(draw.pairwise.venn(dim(track1.df)[1], dim(track2.df)[1], sum(countOverlaps(track1.granges,track2.granges))))
plot(draw.pairwise.venn(dim(track1.df)[1], dim(tss.df)[1], sum(countOverlaps(track1.granges,tss.granges))))
plot(draw.pairwise.venn(dim(track2.df)[1], dim(tss.df)[1], sum(countOverlaps(track2.granges,tss.granges))))

draw.pairwise.venn(area1 = dim(test1.bt1)[1], area2 = dim(john.stam)[1], cross.area = sum(countOverlaps(test1.bt1.granges,john.stam.granges,minoverlap=10)),
                   category = c('ATAC-seq','DNase-seq'), fill = c("light blue", "pink"), cat.dist = c(0.02,0.02), cat.pos = c(0.025, 2),
                   alpha = rep(0.5, 2))
                                                                         
pn.area <- 10000
pn.dist <- 10000
pn.jacc <- 10000

track1.track2.corr <- GenometriCorrelation(track1.granges, track2.granges, 
                                           chromosomes.to.proceed=c('Ath_chr1','Ath_chr2','Ath_chr3','Ath_chr4','Ath_chr5'),
                                           ecdf.area.permut.number = pn.area,
                                           mean.distance.permut.number = pn.dist,
                                           jaccard.measure.permut.number = pn.jacc,
                                           keep.distributions = TRUE)

track1.tss.corr <- GenometriCorrelation(track1.granges, tss.granges, 
                                           chromosomes.to.proceed=c('Ath_chr1','Ath_chr2','Ath_chr3','Ath_chr4','Ath_chr5'),
                                           ecdf.area.permut.number = pn.area,
                                           mean.distance.permut.number = pn.dist,
                                           jaccard.measure.permut.number = pn.jacc,
                                           keep.distributions = TRUE)

track2.tss.corr <- GenometriCorrelation(track2.granges, tss.granges, 
                                        chromosomes.to.proceed=c('Ath_chr1','Ath_chr2','Ath_chr3','Ath_chr4','Ath_chr5'),
                                        ecdf.area.permut.number = pn.area,
                                        mean.distance.permut.number = pn.dist,
                                        jaccard.measure.permut.number = pn.jacc,
                                        keep.distributions = TRUE)
