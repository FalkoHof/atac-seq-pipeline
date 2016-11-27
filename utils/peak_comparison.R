require("GenomicRanges")
require("GenometriCorr")
require('ggbio')
require('VennDiagram')
require('limma')
require('seqLogo')

options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)

#path of track 1
#track1.path = args[1]

#path of track 2
#track2.path = args[2]

track1.path = '/home/falko/projects/atac-seq-pipeline/data/sub_nucl_peaks/test_1.unique_subnucl_fseq.npf'
track2.path = '/home/falko/projects/atac-seq-pipeline/data/sub_nucl_peaks/bc_50_proto_160416.unique_subnucl_fseq.npf'

track1.name <- unlist(strsplit(basename(track1.path), "[.]"))[1]
track2.name <- unlist(strsplit(basename(track2.path), "[.]"))[1]


#not sure if labeling of col 7-10 is correct
fseq.colnames <-c('chromosome','start','end','peak_id','score','strand', 'density', 'thickStart','thickEnd', 'signal')


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

track1.granges <- GRanges(Rle(track1.df$chromosome),
                          IRanges(start = track1.df$start, end = track1.df$end),
                          peakId = track1.df$peak_id, 
                          score = track1.df$score, 
                          seqlengths=tair10.golden.path.length)

track2.granges <- GRanges(Rle(track2.df$chromosome),
                          IRanges(start = track2.df$start, end = track2.df$end), 
                          peakId = track2.df$peak_id, 
                          score = track2.df$score, 
                          seqlengths=tair10.golden.path.length)


plot(draw.pairwise.venn(dim(track1.df)[1], dim(track2.df)[1], sum(countOverlaps(track1.granges,track2.granges))))
plot(draw.pairwise.venn(dim(track1.df)[1], dim(tss.df)[1], sum(countOverlaps(track1.granges,tss.granges))))
plot(draw.pairwise.venn(dim(track2.df)[1], dim(tss.df)[1], sum(countOverlaps(track2.granges,tss.granges))))




draw.pairwise.venn(area1 = dim(track1.df)[1], area2 = dim(track2.df)[1], cross.area = sum(countOverlaps(track1.granges,track2.granges,minoverlap=1L)),
                   category = c(track1.name,track2.name), fill = c("light blue", "pink"), cat.dist = c(0.02,0.02), cat.pos = c(0.025, 2),
                   alpha = rep(0.5, 2))


draw.pairwise.venn(area1 = dim(track1.df)[1], area2 = dim(track2.df)[1], cross.area = sum(countOverlaps(track1.granges,track2.granges,)),
                   category = c('ATAC-seq','DNase-seq'), fill = c("light blue", "pink"), cat.dist = c(0.02,0.02), cat.pos = c(0.025, 2),
                   alpha = rep(0.5, 2))





draw.pairwise.venn(area1 = dim(test1.bt1)[1], area2 = dim(john.stam)[1], cross.area = sum(countOverlaps(test1.bt1.granges,john.stam.granges,minoverlap=10)),
                   category = c('ATAC-seq','DNase-seq'), fill = c("light blue", "pink"), cat.dist = c(0.02,0.02), cat.pos = c(0.025, 2),
                   alpha = rep(0.5, 2))
                                                                         
pn.area <- 1000
pn.dist <- 1000
pn.jacc <- 1000

track1.track2.corr <- GenometriCorrelation(track1.granges, track2.granges, 
                                           chromosomes.to.proceed=c('Ath_chr1','Ath_chr2','Ath_chr3','Ath_chr4','Ath_chr5'),
                                           ecdf.area.permut.number = pn.area,
                                           mean.distance.permut.number = pn.dist,
                                           jaccard.measure.permut.number = pn.jacc,
                                           keep.distributions = TRUE)

track1.track2.corr

