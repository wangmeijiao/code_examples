

import pybedtools
import sys
import copy



#a = pybedtools.example_bedtool('a.bed')

afile = sys.argv[1]
bfile = sys.argv[2]
step = sys.argv[3]
step = int (step)

a = pybedtools.BedTool(afile)
b = pybedtools.BedTool(bfile)

#a = pybedtools.BedTool("gene.bed")
#b = pybedtools.BedTool("tigr7_RM.DNAt_CACTA-EnSpm.bed")

a_interval = a.as_intervalfile()
b_interval = b.as_intervalfile()

#query = a[0]
#b_interval.all_hits(query)
#b_interval.count_hits(query)



## the left flanking of the left part
result1 = [0,0,0,0,0,0,0,0,0,0]
for i in a:
    #print 'Search for region: ' + i.chrom + '\t' + ('%d' %i.start) + '\t' +\
    #('%d' %i.stop) + '\t' + i.name  + '\t' + i.score + '\t' + i.strand
    region = copy.copy(i)
    #print 'Search for region: ' + region.chrom + '\t' + ('%d' %region.start) + '\t' +\
    #('%d' %region.stop) + '\t' + region.name  + '\t' + region.score + '\t' + region.strand
    for j in range(1,11):
        region.start = i.start - j*step
        if(region.start <= 0):
            region.start = 0
        region.end = i.start - (j-1)*step  
        cnt = b_interval.count_hits(region)
        result1[j-1]+=cnt
        #print 'Found %d hits at %d round when region.start is %d region.end is %d ; i.start is %d i.end is %d' %(cnt,j,region.start,region.end,i.start,i.end)



result1.reverse()
for i in result1:
    print i


## the right flanking of the left part
result2 = [0,0,0,0,0,0,0,0,0,0]
for i in a:
    #print 'Search for region: ' + i.chrom + '\t' + ('%d' %i.start) + '\t' +\
    #('%d' %i.stop) + '\t' + i.name  + '\t' + i.score + '\t' + i.strand
    region = copy.copy(i)
    #print 'Search for region: ' + region.chrom + '\t' + ('%d' %region.start) + '\t' +\
    #('%d' %region.stop) + '\t' + region.name  + '\t' + region.score + '\t' + region.strand
    for j in range(1,11):
        region.start = i.start + (j-1)*step
        #if(region.start <= 0):
        #    region.start = 0
        region.end = i.start + j*step
        cnt = b_interval.count_hits(region)
        result2[j-1]+=cnt
        #print 'Found %d hits at %d round when region.start is %d region.end is %d ; i.start is %d i.end is %d' %(cnt,j,region.start,region.end,i.start,i.end)


for i in result2:
    print i


print "0"
#print "0"
#print "0"


## the left flanking of the right part
result3 = [0,0,0,0,0,0,0,0,0,0]
for i in a:
    #print 'Search for region: ' + i.chrom + '\t' + ('%d' %i.start) + '\t' +\
    #('%d' %i.stop) + '\t' + i.name  + '\t' + i.score + '\t' + i.strand
    region = copy.copy(i)
    #print 'Search for region: ' + region.chrom + '\t' + ('%d' %region.start) + '\t' +\
    #('%d' %region.stop) + '\t' + region.name  + '\t' + region.score + '\t' + region.strand
    for j in range(1,11):
        region.start = i.end - j*step
        if(region.start <= 0):
            region.start = 0
        region.end = i.end - (j-1)*step  
        cnt = b_interval.count_hits(region)
        result3[j-1]+=cnt
        #print 'Found %d hits at %d round when region.start is %d region.end is %d ; i.start is %d i.end is %d' %(cnt,j,region.start,region.end,i.start,i.end)



result3.reverse()
for i in result3:
    print i



## the right flanking of the right part
result4 = [0,0,0,0,0,0,0,0,0,0]
for i in a:
    #print 'Search for region: ' + i.chrom + '\t' + ('%d' %i.start) + '\t' +\
    #('%d' %i.stop) + '\t' + i.name  + '\t' + i.score + '\t' + i.strand
    region = copy.copy(i)
    #print 'Search for region: ' + region.chrom + '\t' + ('%d' %region.start) + '\t' +\
    #('%d' %region.stop) + '\t' + region.name  + '\t' + region.score + '\t' + region.strand
    for j in range(1,11):
        region.start = i.end + (j-1)*step
        #if(region.start <= 0):
        #    region.start = 0
        region.end = i.end + j*step
        cnt = b_interval.count_hits(region)
        result4[j-1]+=cnt
        #print 'Found %d hits at %d round when region.start is %d region.end is %d ; i.start is %d i.end is %d' %(cnt,j,region.start,region.end,i.start,i.end)

#print "\n"

for i in result4:
    print i



'''
    hits = b_interval.all_hits(i)
    cnt = 0
    for j in hits:
        cnt+=1
        print ('Found %d hits:\n' %cnt) + j.chrom + '\t' + ('%d' %j.start) +'\t' + ('%d' %j.end)
'''


#slice = a[1:3]
#print( slice)

#print (b - a).count()

#for i in a:
#    print i



