import sys
from sys import argv
#from sys import stderr
import re

# gene island makes sense only in heterochromatin regions
# filter < 500bp GIs
# split by chr 


hetfile,gifile,outfile,filter = argv[1],argv[2],argv[3],int(argv[4])
#print filter

#print argv[1:3]
#exit()

#print hetfile,gifile,"\n"
#exit()
try:
    fh_het = open(hetfile,"r")
except IOError:
    print >> sys.stderr, "err open file", hetfile
else:
    print >> sys.stderr, "read for file ", fh_het.name

'''
#use while loop
line= fh_het.readline()
while line :

    print line,
    line = fh_het.readline()

'''

#####sub#########

def readBed(file):
#for genearl bed 4 file
    fh = open(file,'r')
    bed = []
    for line in fh:
        line = line.strip()
        if len(line) == 0 or line.startswith("#"):
            continue
        chr,s,e,id = re.split("[\t ]+",line)
        s = int(s)
        e = int(e)
        #print filter
        if e - s + 1 >= filter:
            bed.append([chr,s,e,id])
    fh.close()
    return bed



########main ############

#use for loop directly to readin het regions
region = {}
for line in fh_het.readlines():
    line = line.strip()
    #print line,
    if len(line) == 0 or line.startswith('#'):
        #print "empty line",line
        continue
    chr,s,e = re.split("[\t ]+",line)
    s=int(s)
    e=int(e)
    #print chr,s,e
    region[chr] = (s,e)
    #print chr,region[chr]

#print region
#print region["chr04"][2]
fh_het.close()


#readin GI bed, filter and output
beds = readBed(gifile)
fh_out = open(outfile,"w")
for i in beds:
    #print i
    chr,s,e,id = i
    if e < region[chr][0] or s > region[chr][1]:
        pass
    else:
        fh_out.write("%s\t%s\t%i\t%s\n" % (chr,s,e,id) )

fh_out.flush()
fh_out.close()

print >> sys.stderr, "done"


