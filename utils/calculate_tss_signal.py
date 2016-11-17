import metaseq
import pybedtools
import numpy as np
import multiprocessing



from matplotlib import pyplot as plt

tss_annotation = str(sys.argv[1])
bam_file = str(sys.argv[2])
slop_region = str(sys.argv[3])


bins = str(sys.argv[4])



#read in tss annotation file in bed format
tss_bed = pybedtools.BedTool(tss_annotation)

# extend by 1000 bp up/downstream
tss_slop = tss_bed.slop(b=1000, genome = annotationDir + 'chromLength')


bam_gsignal = metaseq.genomic_signal(bam_file, 'bam')

# the region +/-500bp around each TSS will be split into a total of 100 bins,
# change as needed

x = np.linspace(-1000, 1000, bins)

# most of the work happens here
test1_tss = test1_bam.array(tss_slop, bins = bins, processes = cpus)
test2_tss = test2_bam.array(tss_slop, bins = bins, processes = cpus)
bc1_tss = bc1_bam.array(tss_slop, bins = bins, processes = cpus)
bc2_tss = bc2_bam.array(tss_slop, bins = bins, processes = cpus)
mg_tss = mg_bam.array(tss_slop, bins = bins, processes = cpus)
tp_tss = tp_bam.array(tss_slop, bins = bins, processes = cpus)
