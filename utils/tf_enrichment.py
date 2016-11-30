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
        n = sample_dict[sample]

        stats.fisher_exact([8,2])




        Atlantic  Indian
whales     8        2
sharks     1        5


        open_chrom  genome
whales     8        2
sharks     1        5


We use this table to find the p-value:
>>>

>>> import scipy.stats as stats
>>> oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])


>>> pvalue
0.0349...







sample_dict, sample_size =read_file(sample)
bg_dict, background_size = read_file(background)
