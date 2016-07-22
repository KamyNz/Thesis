# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 20:52:19 2015

@author: CamilaMV
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import sys

## reading genbank file

genome_name=sys.argv[1]

genbank_file=SeqIO.parse(open(sys.argv[2],"r"),"genbank")

#output file where the rnas are going to be written

output_file=open(genome_name,"w")


#creatign RNA list

rRNAs=[]

for genome in genbank_file:
    
    for feature in genome.features:
        
        if(feature.type == "rRNA"):
            
            #print feature
            
            start = feature.location.start
            
            end = feature.location.end
            
            location = str(start) + ":" + str(end)
            
            
            desc = feature.qualifiers['locus_tag'][0]            
            
            product = feature.qualifiers['product'][0]
            
            seq = feature.extract(genome.seq)
            
#            location=[start,end]
#            
#            print location
#            
#            type(location)
            
            #print seq
            
            record = SeqRecord(seq,id=desc,description=product + " " + location)
#            
#            # in record variable ID, desc and seq are together
            rRNAs.append(record)
#            
#            output_file.write(desc)
#
#            output_file.write("\t")
#
#            output_file.write(product)
#
#            output_file.write("\t")
#
#            output_file.write(location)
#            
#            output_file.write("\n")

            
SeqIO.write(rRNAs,output_file,"fasta")
output_file.close()


