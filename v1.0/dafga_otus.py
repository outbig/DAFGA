#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
# Author:   Yongkyu Kim, PhD
            Max Planck Institute for terrestrial microbiology
# 
# Date:     2014-04-02
# Version:  1.0
# 
# Try 'dafga_otus.py -h' for more information and See manual
# 
# 
# Bugs: Please report to https://github.com/outbig/DAFGA/issues?state=open
"""

import subprocess, os, re, glob
from argparse import ArgumentParser
from Bio import SeqIO
import pandas as pd

def make_subdirectory(dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    return os.path.abspath(dirname)
    
def preclustering(filename,nt_id): # running usearch to generate consensus seq of preclustering
    print "\n   ...Pre-Clustering of amplicon sequences(nt) with %s identity...\n" % nt_id 
    subdir = make_subdirectory(o_dir+"/preclustering")    
    name_head = subdir+"/preclstr_"+nt_id
    prog = subprocess.Popen(["usearch","-cluster_fast",filename,"-id",nt_id,"-centroids",name_head+"_centroids.fasta",
                             "-consout",name_head+"_consensus.fasta","-uc",name_head+".uc"])
    subprocess.Popen.wait(prog)
    print "\nDone\n"
    return name_head+"_consensus.fasta", name_head+".uc"
    
def translate_consensus(cons_seqs, gen_code):
    codon = {1:"Standard Code",2:"Vertebrate Mitochondrial Code",3:"Yeast Mitochondrial Code",
             4:"Mold, Protozoan, and Coelenterate Mitrochondrial Code and Mycoplasma/Spiroplasma Code",
             5:"Invertebrate Mitochondrial Code",6:"Ciliate, Dasyclatdacean and Hexamita Nuclear Code",
             9:"Echinoderm and Flatworm Mitochondrial Code",10:"Euplotid Nuclear Code",
             11:"Bacterial, Archaeal and Plant Plastide Code",12:"Alternative Yeast Nuclear Code",
             13:"Ascidian Mitochondrial Code",14:"Alternative Flatworm Mitochondrial Code",
             15:"Blepharisma Nuclear Code",16:"Chlorophycean Mitochondrial Code"}
    print "   ...Translation into protein sequence..."
    print "Representative consensus sequences are translated with {0}\n".format(codon[gen_code])
    aa_dir = make_subdirectory(o_dir+"/clustering")
    with open(aa_dir+"/translated_consensus.fasta","w") as faa, open(cons_seqs,"r") as cons:    
        for count, record in enumerate(SeqIO.parse(cons,"fasta")):
            faa.write(">precluster"+str(count)+"\t"+str(record.id)+"\n")
            seq = trim_seqs(record.seq)
            faa.write(str(seq.translate(table=gen_code))+"\n")
    print 'Done\n'
    return aa_dir+"/translated_consensus.fasta"

def rep_aa_seqs(rep_otu, faa):
    print '\n   ...Pick representative seq of OTUs...\n'
    new_fasta = os.path.join(os.path.dirname(faa),"../rep_seqs.fasta")
    with open(faa,"r") as in_clstr, open(new_fasta,"w") as new_rep:
        for record in SeqIO.parse(in_clstr,"fasta"):
            if record.id in rep_otu.keys():
                rep_preclstr = record.description.split("\t")[0]
                seq_id = re.search("centroid=(.+?);seq",record.description).group(1)
                new_rep.write(">"+rep_otu[record.id]+"\t"+rep_preclstr+"\t"+seq_id+"\n")
                new_rep.write(str(record.seq)+"\n")
    print 'Done\n'
    return new_fasta
    
def trim_seqs(seq):
    if tl == '+1':
        seq = end_trim(seq) 
    elif tl == '+2':
        seq = end_trim(seq[1:])
    elif tl == '+3':
        seq = end_trim(seq[2:])
    return seq 
    
def end_trim(seq):
    if len(seq)%3 == 0:
        pass
    elif len(seq)%3 == 1:
        seq = seq[:-1]
    elif len(seq)%3 == 2:
        seq = seq[:-2]
    return seq

def clustering_aa(faa):
    print '\n   ...Clustering of translated representative sequences...'
    threshold = pd.read_csv(corr, delimiter='\t',index_col='Unnamed: 0')
    iden_cut, siml_cut = float(threshold[taxa][1])/100, float(threshold[taxa][2])/100
    print '     with {0} identity using usearch\n'.format(iden_cut)
    clstr, clstr_rep = run_usearch(faa,iden_cut)
    clstr_otu,rep_otu = parse_usearch(clstr,"cluster")
    print 'Done\n'
    return rep_otu, clstr_otu

def run_usearch(filename, aa_id):
    name_head = os.path.dirname(filename)   
    cen, uc = name_head+"/cluster_{0}_centroids.fasta".format(aa_id), name_head+"/clusters_{0}.uc".format(aa_id)
    prog = subprocess.Popen(["usearch","-cluster_fast",filename,"-id",str(aa_id),"-centroids", cen,"-uc",uc])
    subprocess.Popen.wait(prog)
    return uc, cen

def parse_usearch(clstr, clstr_level):
    base = os.path.splitext(clstr)[0]
    with open(base+"_OTUs.txt","w") as otu, open(clstr,"r") as cluster:
        clstr_dict,rep_dict = {}, {}
        for line in cluster:            
            temp = line.split("\t")
            if temp[0] == "S":
                clstr_dict[temp[-2]]=[temp[-2]]
            elif temp[0] == "H":
                clstr_dict[temp[-1].rstrip("\n")].append(temp[-2])
            elif temp[0] == "C":
                otu.write(clstr_level+temp[1])
                for seq in clstr_dict[temp[-2]]:
                    otu.write("\t"+seq)
                    rep_dict[temp[-2]]="cluster"+temp[1]
                otu.write("\n")
        return base+"_OTUs.txt",rep_dict

def OTU_table(preclstr,clstr):
    print '\n   ...Create OTU table...\n'
    p_dir = os.path.join(os.path.dirname(preclstr),"..")
    preclstr_dict, clstr_dict = {}, {}
    with open(preclstr,"r") as preOTU, open(clstr,"r") as OTU, open(os.path.join(p_dir,"OTU_mapping.txt"),"w") as seqOTU:
        for line in preOTU:
            temp = line.strip("\n").split("\t",1)    
            preclstr_dict[temp[0]] = temp[1]
        for line in OTU:
            temp = line.strip("\n").split("\t")
            clstr_dict[temp[0]] = temp[1:]
        for key in clstr_dict.keys():
            seqOTU.write(key)
            for value in clstr_dict[key]:
                seqOTU.write("\t"+preclstr_dict[value])
            seqOTU.write("\n")
    print 'Done\n'
    

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i",dest="fasta",metavar="FASTA",required=True, help="input amplicon sequences resulting from split_library.py")
    parser.add_argument("-o",dest="o_dir",metavar="output directory",required=True, help="the output directory where output files are saved")
    parser.add_argument("-t",dest="tl",metavar='translation',required=True, choices=['+1','+2','+3'],help='start of translation in your amplicons')
    parser.add_argument('-c',dest='corr',metavar='correlation',required=True, help='the output directory of the correlation_plot.py')    
    parser.add_argument("--nt_id",dest="nt_id",metavar="NT_SEQ_IDENTITY",default="0.97",help="sequence identity threshold for preclustering using nucleotide sequences [0-1]. Default: 0.97")
    parser.add_argument('--taxa',dest='taxa',default='species',choices=['phylum','class','order','family','genus','species','strain'],help='selected taxonomic rank for reliable assignment. Default: species')
    parser.add_argument("--gen_code", dest="gen_code",metavar="CODON",type=int, choices=range(1,7)+range(9,17),default="1",help="NCBI codon table for translation [1-6 and 9-16]. Default: 11 for bacteria and archaea")
    args=parser.parse_args()
    
    global o_dir,tl, taxa, corr
    o_dir = make_subdirectory(args.o_dir)
    tl, taxa = args.tl, args.taxa
    corr = glob.glob(os.path.join(args.corr,'*_to_16S.txt'))[0]
    cons_seqs, preclstr = preclustering(args.fasta, args.nt_id)  
    preclstr_otu,pre_rep_otu = parse_usearch(preclstr,"precluster")
    faa = translate_consensus(cons_seqs, args.gen_code)
    rep_otu, clstr_otu = clustering_aa(faa)
    rep_aa_seqs(rep_otu,faa)
    OTU_table(preclstr_otu,clstr_otu)
