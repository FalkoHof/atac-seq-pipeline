# script to offset reads aligning to the + strand by +4 bp, and all reads aligning to the - strand by - 5 bp.
# provide input file as first argument, file output name as second, the script takes bed files as input
# change the values in dictionary genome_size according to the reference genome chromosome names and sizes used for mapping.
# to log output of the file simply pipe std out to a txt file
import sys
import os.path

filenameInput = str(sys.argv[1])
filenameOutput = str(sys.argv[2])
#needs to be set according to reference genome that was used for mapping
genome_size ={'Ath_chr1':30427671,'Ath_chr2':19698289,'Ath_chr3':23459830,
              'Ath_chr4':18585056,'Ath_chr5':26975502,'Ath_chrc':154478,'Ath_chrm':366924}

if not os.path.isfile(filenameInput):
    sys.exit('Input file does not exist!')

if filenameInput == filenameOutput:
    sys.exit('Input file cannot be output file!')

fin = open(filenameInput,'r')
fout = open(filenameOutput, 'w')

while True:
    line = fin.readline()

    #breaks the loop when file is read
    if not line:
        break

    col = line.split('\t')

    start = int(col[1])
    end = int(col[2])
    strand = col[5]

    #adds strand specific offset and warn if position would exceed chromosome size
    if '+' in strand:
        start+=4
        end+=4
        if genome_size[col[0]] < start:
            start = genome_size[col[0]]
            print('Warning: ' + col[3] + 'start position ' + str(start) + ' > '  + col[0] + str(genome_size[col[0]]))
        if genome_size[col[0]] < end:
            end = genome_size[col[0]]
            print('Warning: ' + col[3] + 'end position ' + str(end) + ' > ' + col[0] + str(genome_size[col[0]]))

    else:
        start-=5
        end-=5
        if start < 0:
            start=0
            print('Warning: ' + col[3] + 'start position ' + str(start) + '< 0'  + col[0])
        if end < 0:
            end = 0
            print('Warning: ' + col[3] + 'end position ' + str(end) + '< 0'  + col[0])


    fout.write(col[0] + '\t' + str(start) + '\t' + str(end) + '\t' + col[3] + '\t' + col[4] + '\t' + col[5])

fin.close()
fout.close()
print('Done!')
