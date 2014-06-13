__author__ = 'seidz'
#!/usr/bin/env python

import sys

current_kmer = None
kmer = None
infile=open('sampleoutput0lasts.txt','r')
infile=sorted(infile)
outfile = open('sampleoutput2lasts.txt','w')
# input comes from STDIN
for line in infile:
    # remove leading and trailing whitespace
    line = line.strip()
    #print line
    # parse the input we got from mapper.py
    kmer, v1 = line.split('\t')

    # convert count (currently a string) to int
    try:
        v1 = int(v1)
    except ValueError:
        # count was not a number, so silently
        # ignore/discard this line
        print "error"
        continue

    # this IF-switch only works because Hadoop sorts map output
    # by key (here: word) before it is passed to the reducer
    if current_kmer == kmer:
        current_v1 += v1
    else:
        if current_kmer:
            # write result to STDOUT
            outfile.write('%s\t%s\n' % (current_kmer, current_v1))
        current_kmer = kmer
        current_v1 = v1



# do not forget to output the last word if needed!
if current_kmer == kmer:
    outfile.write('%s\t%s\n' % (current_kmer, current_v1))
