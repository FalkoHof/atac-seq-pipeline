import scipy.stats as stats

def read_file(f):
    print ("Processing file:" + f)
    fin = open(f)
    lines = fin.readlines
    dist,size = process_lines(lines)
    print ("Processing file:" + f + " -Done")
    return d

def process_lines(lines,col_idx):
    d = dict()
    for line in lines:
        cols = line.split('\t')
        key = col[col_idx]
        d = addCountToDict(d,key)
    return d, len(lines)

#add +1 to a dict holding key, int value, othervise initilizea and set value =1
def addCountToDict(d, key):
    if fragmentSize in d:
        d[key] += 1
    else:
        d[key] = 1
    return d

def hypergeom_test(sample_dict, sample_size, background_dict, background_size):
    header=['ID', 'p_value','oddsratio','#obs_sample', 'adj_sample_size', \
            '#obs_backgroud','adj_bg_size']
    results=[]
    for sample in sample_dict:
        tf_sample = sample_dict[sample]
        tf_background = background_dict[sample]

        table = ([[tf_sample, sample_size-tf_sample],
                 [tf_background, background_size-tf_background]])
# matrix for testing
#             tf_x                                all_tfs
# open_chrom  tf_sample                           sample_size-tf_sample
# background  tf_background                       background_size-tf_background
#
        oddsratio, p_value = stats.fisher_exact(table, 'greater')

        values = [v for r in table for v in r]
        one_res = [sample, p_value, oddsratio]
        one_res.append(values)

    return header, results


def writeToFile(header,data,f):
    fout = open(f,'w')
    fout.write('\t'.join(header))

    for d in data:
        fout.write('\t'.join(map(str, d)))

    fout.close()

sample_file = sys.argv[1]
bg_file = sys.argv[2]
fout = sys.argv[3]

sample_dist, sample_dist = read_file(sample_file)
bg_dist, bg_size = read_file(bg_file)


header, results = hypergeom_test(sample_dist,
                                 sample_size,
                                 background_dict,
                                 background_size)

writeToFile(header,results,fout)
