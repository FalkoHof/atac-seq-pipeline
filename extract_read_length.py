#script to extract the lenght of each sequenced mate pair
#works only if the file is sorted according to the mate pair id
#filenameInput = 'test_1_mate_sorted.bed'

#script to extract the lenght of each sequenced mate pair
#works only if the file is sorted according to the mate pair id
#first argument input file, second argument is output file
import sys
import os.path
import argparse
from sets import Set

parser = argparse.ArgumentParser(description='Calculates distance between between paired end reads of a BED file, sorted on read ids')
parser.add_argument("-g", "--histogram", action="store_true", help="export fragment counts table",default=False)
parser.add_argument("-d", "--detailed", action="store_true", help="export detailed data in bed like format for each read",default=False)
parser.add_argument("-v", "--verbose",action="store_true", help="information about the run is printed to stdout",default=False)
parser.add_argument("-f", "--fragmentSize", help="set threshold for maximal fragment length (default: 2000 bp)",default=2000,type=int)
parser.add_argument("-o", "--output", dest='output', help="set output directory", required=True, type=str)
#parser.add_argument("-n", "--ncbi", action="store_true", help="specifies that the dataset is from ncbi)",default=False)
parser.add_argument("input", help="BED file of the regions you want to generate calculate the distance for")
args = parser.parse_args()

histogram = args.histogram
detailedInfo = args.detailed
fragmentThreshold = args.fragmentSize
filenameInput = args.input
outputDir = args.output
#ncbi_file = args.ncbi

blacklistChr = Set(['Ath_chrc','Ath_chrm'])
histoDict = dict()

#inputFolder = '/Volumes/nodine/lab/members/falko/ATAC-SEQ/500cells_human/'
#filenameInput = inputFolder + '34769.filtered.sorted.bed'
#filenameOutput = inputFolder + 'read_length/' + 'readLength_34769.tab'

if not os.path.isfile(filenameInput):
    sys.exit('Input file does not exist!')

#generate output filenames. Pattern is name.histogram or name.detailed_info
def generateFilenames(filename,outputDir):
    basename = os.path.splitext(filename)[0]
    outputPath = os.path.join(outputDir,basename)
    filenameHisto = outputPath + '.histogram'
    filenameDetail =  outputPath + '.detailed_info'
    return (filenameHisto, filenameDetail)

#add +1 to a dict holding key, int value, othervise initilizea and set value =1
def addCountToDict(d, fragmentSize):
    if fragmentSize in d:
        d[fragmentSize] += 1
    else:
        d[fragmentSize] = 1
    return d

#write dictionary key/value int pairs within a certain range to tab delimed file
def writeHistogramToFile(d,filename, fragmentSize):
    #open file and write header
    fout = open(filename,'w')
    fout.write('FragmentSize' + '\tCounts\n')
    #loop over defined range, if no entry is present write out 0
    for i in range(1,fragmentSize+1):
        count = 0
        if i in d:
            count = d[i]
        fout.write(str(i) + '\t' + str(count) + '\n')
    fout.close()

#get return tulple with max, min positional values of 2 bed file entries
def getMaxMinValues(read1, read2):
    positionList = []
    positionList.append(int(read1[1]))
    positionList.append(int(read1[2]))
    positionList.append(int(read2[1]))
    positionList.append(int(read2[2]))
    maxPosition = max(positionList)
    minPosition = min(positionList)
    return (maxPosition, minPosition)

#calculates the distance between two bed file features in bp
def getReadLength(read1, read2):
    #get max and min position
    maxPosition, minPosition = getMaxMinValues(read1, read2)
    #caltulate read length and return
    readLength = maxPosition - minPosition
    return (readLength)

#writes a line of detailed info to file, if verbose prints info
def writeDetailedInfo(fout,id, chr, length, pos1, pos2):
    fout.write(id + '\t' + str(length) + '\t' + chr + '\t' + str(pos1) + '\t' + str(pos2) +'\n')

#generate filenames for output
filenameHisto, filenameDetail = generateFilenames(filenameInput, outputDir)

#open input file
fin = open(filenameInput,'r')
print('InputFile: ' + filenameInput)

#conditional output file
if detailedInfo:
    detailFout = open(filenameDetail,'w')

#some counters for transparency
lineCounter = 0
mateCounter = 0
blacklistCounter = 0

#loop ove the input file as long as there are lines
while True:
    mate1 = fin.readline().rstrip('\n')
    mate2 = fin.readline().rstrip('\n')

    #breaks the loop when file is read
    if not mate2 or not mate1:
        if not mate2:
            #counts the additional last line when uneven number of mates
            lineCounter+=1
        if args.verbose:
            print('Lines read: '+ str(lineCounter))
        break

    lineCounter+=2

    mate1Cols = mate1.split('\t')
    mate2Cols = mate2.split('\t')

    #if ncbi:
    #mate1Id = (mate1Cols[3])
    #mate2Id = (mate2Cols[3])

    #else:
    mate1Id = (mate1Cols[3])[:-2]
    mate2Id = (mate2Cols[3])[:-2]


    #loop that skips 1 line as long ids are different, in case only 1 mate is present
    while not mate1Id == mate2Id:
        if args.verbose:
            print('Warning mate1 and 2 are not identical!\t' + 'mate1Id\t' + 'mate2Id\t' + 'Line Number: ' + str((lineCounter-1))+ '-' + str(lineCounter) + '\tline skipped...')
        mate1 = mate2
        mate2 = fin.readline().rstrip('\n')

        if not mate2:
            if args.verbose:
                print('Lines read: '+ str(lineCounter))
            break

        lineCounter+=1

        mate1Cols = mate1.split('\t')
        mate2Cols = mate2.split('\t')

        #if ncbi:
        mate1Id = (mate1Cols[3])
        mate2Id = (mate2Cols[3])
        #else:
        #    mate1Id = (mate1Cols[3])[:-2]
        #    mate2Id = (mate2Cols[3])[:-2]

    chrMate1 = mate1Cols[0]
    chrMate2 = mate2Cols[0]
        #check if chromosome has been blacklisted
    chromosomes = [chrMate1, chrMate2]
        #check if any of the chromosomes has been blacklisted
    if not any(x in blacklistChr for x in chromosomes):
        #check if both chromosomes are the same
        if chrMate1 == chrMate2:
            mateCounter+=1
            readLength = getReadLength(mate1Cols,mate2Cols)
            #print some warnings of length is longer than threshold
            if readLength > fragmentThreshold and args.verbose:
                print('Warning: ' + mate1Id + 'read is longer than ' + str(fragmentThreshold) +'bp.')
            #make histogram dict, if desired
            if histogram:
                histoDict = addCountToDict(histoDict,readLength)
                #generate detailed output, if desired
            if detailedInfo:
                maxVal, minVal = getMaxMinValues(mate1Cols,mate2Cols)
                writeDetailedInfo(detailFout,chrMate1,mate1Id,readLength,minVal,maxVal)
            #writeDetailedInfo(detailFout,mateChr,mate1Id[:-2],readLength,minVal,maxVal)
            #if chromosomes are not the same, skip. probably error in mapping
        else:
            print('Warning: ' + mate1Id + ' reads map on differen chromosomes! ' + \
            str(mate1Cols[0]) +'\t' + mate2Cols[0]+ 'Skipping...')
    else:
        blacklistCounter+=1

#write histogram to file, give some overview statistics and close everything
if histogram:
    writeHistogramToFile(histoDict,filenameHisto, fragmentThreshold)
#substract one becuase of reading asymetry
lineCounter-=1
print('Number of mate pairs: ' + str(mateCounter) + '\tNumber of reads: ' + str(lineCounter))
print('Number of blacklisted pairs: ' + str(blacklistCounter) + '\tNumber of reads: ' + str(lineCounter))
print('Percentage of considered pairs: ' + str(((mateCounter*1.0)/(lineCounter/2.0)*100)))
print('Number of pairs not accounted by blacklist: ' + str((lineCounter/2.0)-mateCounter-blacklistCounter))
fin.close()
if detailedInfo:
    detailFout.close()
