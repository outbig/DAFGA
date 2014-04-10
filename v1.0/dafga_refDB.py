#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
# Author:   Yongkyu Kim, PhD
            Max Planck Institute for terrestrial microbiology
# 
# Date:     2014-04-02
# Version:  1.0
# 
# Try 'dafga_refDB.py -h' for more information and See manual
# 
# Purpose:   
# Bugs: Please report to https://github.com/outbig/DAFGA/issues?state=open
"""
from Bio import SeqIO, Entrez
from argparse import ArgumentParser
import os

def make_directory(dirname):
    try:
        os.makedirs(dirname)
    except OSError:
        raise
    return os.path.abspath(dirname)

def delete_contig_line_from_gbk(gp): # delete unreconizable lines by bioptyhon 
    processed = os.path.join(o_dir,"processed.gp")
    with open(processed,"w") as new_gp:
        for line in open(gp,"r"):
            if line.startswith("CONTIG"):
                pass
            elif line.startswith("     SecStr"):
                pass
            elif line.startswith("     Het"):
                pass
            elif line.startswith("     Bond"):
                pass
            else:
                new_gp.write(line)
    return processed
  
def parsing_gp(gp, LENGTH):
    print '\n   ...Parsing reference sequnece information in gp format...\n'
    source_strain, ref_taxa_xref = {}, {}
    with open(prefix+"_ref_seqs.fasta","w") as refseq:
        for record in SeqIO.parse(gp,"genbank"):
            source = record.features[0].qualifiers
            if "culture" not in record.description:
                source_strain_info = ""
                taxa_xref = [x for x in source["db_xref"] if "taxon:" in x]
                xref = taxa_xref[0].split("taxon:")[-1]
                ref_taxa_xref[record.name] = xref 
                refseq.write(">"+record.name+"\t"+record.description.strip(".")+"\t"+xref+"\n")
                refseq.write(str(record.seq)+"\n")
                if len(record.seq) >= int(LENGTH):
                    if source.get("strain"):
                        if source["strain"][0] in source["organism"][0]:
                            source_strain_info = source["organism"][0]
                        else:
                            source_strain_info = source["organism"][0]+" "+source["strain"][0].split(" ")[-1]
                        source_strain[record.name] = [source_strain_info, record.seq]
    print 'Done\n'
    return ref_taxa_xref, source_strain, refseq

def efetch_from_taxonomy(xref): # retreive full taxonomic lineage from taxonomy database
    with open(prefix+"_ID_to_taxonomy.txt","w") as outfile:
        handle_taxon = Entrez.efetch(db="taxonomy", id=xref.values(), retmode="xml")
        records = Entrez.read(handle_taxon)
        handle_taxon.close()
        print "No. of reference sequences: {0}\n".format(len(records))
        print "   ...Retrieve taxonomic lineage information...   "
        outfile.write("ReferenceID\tNCBI_taxanomy\n")
        for i, seq in enumerate(xref.keys()):
            outfile.write(seq+'\t'+records[i]['Lineage']+'\n')
        print "\nDone\n"
        
def searching_nt(source_strain):
    print "  ...Searching nt database and fetching 16S rRNA sequences...\n  "
    print "16S ID\t\tFunGene ID\tDescription\t\t16S length"
    source_list = []
    with open(prefix+"_strain_16S_rRNAs.fasta","w") as rRNAs:
        for seqID, strain in source_strain.items():
            search = strain[0]+"[Organism] AND 16S NOT genome"  # key word to search 16S rRNA of source organism
            handle_esearch = Entrez.esearch(db="nucleotide",term=search)
            records = Entrez.read(handle_esearch)
            handle_esearch.close()
            if len(records["IdList"]) > 0:
                source_list.append(efetch_from_nt_list(seqID, strain, records["IdList"], rRNAs))
    print "Done\n"
    return source_list
    
def efetch_from_nt_list(seqID, strain, nt_list, rRNA_out): 
    handle_efetch = Entrez.efetch(db="nucleotide", id=nt_list, retmode="xml", validate=False)
    records = Entrez.read(handle_efetch, validate=False)
    handle_efetch.close()
    temp = None
    for record in records:
        HIT_DEF = record["GBSeq_definition"]
        if '16' in HIT_DEF and 1200 < int(record["GBSeq_length"]) < 1600 and 'inter' not in HIT_DEF:
            if temp == None:
                if record.get('GBSeq_sequence'):
                    rRNA_out.write('>'+record['GBSeq_locus']+'\t'+seqID+'\t'+HIT_DEF+'\n'+record['GBSeq_sequence']+'\n')
                    temp = seqID
                    print record['GBSeq_locus']+'\t'+seqID+'\t'+HIT_DEF+'\t'+str(len(record['GBSeq_sequence']))
    return  temp
    
def delete_redundancy(rRNA_sources, refseq):    
    with open(prefix+'_strain.fasta','w') as fg_strain:
        count = 0 
        for item in rRNA_sources:
            if item in refseq.keys():
                fg_strain.write('>'+item+'\n'+str(refseq[item][1])+'\n')
                count+=1
        print "No. of refseqs with strain level description: {0}\n".format(count)
            
    
if __name__ == "__main__":    
    parser = ArgumentParser()
    parser.add_argument('-gp',dest='gp',required=True, help='publically available sequences retrieved from NCBI protein database in gp format')
    parser.add_argument('--email',dest='email',required=True, help='to inform NCBI who you are')
    parser.add_argument('-o',dest='o_dir',required=True, help='the directory where output files will be saved')
    parser.add_argument('-l',dest='length', default = 50, help='minimum length of reference sequences. Default: 50')
    args=parser.parse_args()
    Entrez.email = args.email
    global o_dir, prefix
    o_dir = make_directory(args.o_dir)
    prefix = os.path.join(o_dir,os.path.basename(args.gp).split(".")[0])
    gp = delete_contig_line_from_gbk(args.gp)
    ref_xref, source_strain, refseq = parsing_gp(gp, args.length)   
    efetch_from_taxonomy(ref_xref)
    source_list = searching_nt(source_strain)
    delete_redundancy(source_list, source_strain)
    
