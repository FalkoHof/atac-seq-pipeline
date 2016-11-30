import scipy.stats as stats

def read_file(f):
    fin = open(f)
    lines = fin.readlines
    d=process_lines(lines)
    return d

def process_lines(lines,col_idx):
    histoDict = dict()
    for line in lines:
        cols = line.split('\t')
        tf = col[col_idx]
    return histoDict, len(lines)

#add +1 to a dict holding key, int value, othervise initilizea and set value =1
def addCountToDict(d, key):
    if fragmentSize in d:
        d[key] += 1
    else:
        d[key] = 1
    return d

def hypergeom_test(sample_dict, sample_size, background_dict, background_size):

    M = background_size
    N = sample_size
    for sample in sample_dict:
        tf_sample = sample_dict[sample]
        tf_background = background_dict[sample]

        oddsratio,p_value = stats.fisher_exact(([tf_sample, sample_size-tf_sample],
                            [tf_background,background_size-tf_background]),
                            'greater')
#
#             tf_x                                all_tfs
# open_chrom  tf_sample                           sample_size-tf_sample
# background  tf_background                       background_size-tf_background
#


sample_dict, sample_size =read_file(sample)
bg_dict, background_size = read_file(background)
