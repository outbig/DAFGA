#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
# Author:   Yongkyu Kim, PhD
            Max Planck Institute for terrestrial microbiology
# 
# Date:     2014-04-02
# Version:  1.0
# 
# Try 'taxa_assignment.py -h' for more information and See manual
# 
# Purpose:   
# Bugs: Please report to https://github.com/outbig/DAFGA/issues?state=open

"""

from argparse import ArgumentParser
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO, SeqIO
import os, subprocess, errno, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as grd
import pandas as pd

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    return path
    
def makeBLASTdb(ref_dir):
    print '\n  ...Construct BLAST DB with reference seqs...'
    ref_seqs =  glob.glob(os.path.join(os.path.abspath(ref_dir),'*_ref_seqs.fasta'))[0]
    db_name = os.path.basename(ref_seqs).split('_seqs')[0]
    db_out = os.path.join(o_dir,db_name)
    prog = subprocess.Popen(["makeblastdb","-in",ref_seqs,"-dbtype","prot","-out",db_out])
    subprocess.Popen.wait(prog)
    print '\nDone\n'
    return db_out
    
def runBLASTP(ref_db, thread,evalue):
    print '\n   ...BLASTP of representative seqs to reference DB...'
    base = os.path.basename(ref_db).split('ref')[0]
    xml_out = os.path.join(o_dir,base+'BLASTP.xml')
    prog = subprocess.Popen(["blastp","-query",rep_seqs,"-db",ref_db,"-outfmt","5","-out", xml_out,"-evalue",str(evalue),"-num_threads",str(thread),"-max_target_seqs","1"])
    subprocess.Popen.wait(prog)
    print '\nDone\n'
    return xml_out
    
def parsing_BLAST(blast,ref_dir):
    print '   ...Parsing BLAST output...\n'
    id_to_taxa = glob.glob(os.path.join(os.path.abspath(ref_dir),'*_ID_to_taxonomy.txt'))[0]
    base = os.path.basename(blast).split('BLASTP')[0]
    with open(os.path.join(o_dir,base+"taxa_assignment.txt"),"w") as taxa_assignment, open(blast,"r") as blast_xml, open(id_to_taxa,"r") as lineage:
        NO_HIT = ("\t"+"No blast hit")+"\n"    
        taxa_lineage,blast_out = {}, []        
        taxa_assignment.write("OTU_name\tNCBI_taxonomy\tcoverage(%)\talignment_identity(%)\talignment_similarity(%)\tReference_sequence_ID\n")
        for line in lineage:
            temp = line.strip().split("\t")
            taxa_lineage[temp[0]]=temp[1]
        for record in NCBIXML.parse(blast_xml):
            if record.alignments:
                ref_id = str(record.descriptions[0].title.split(" ")[1])
                identity = float(record.alignments[0].hsps[0].identities)
                positive = float(record.alignments[0].hsps[0].positives)
                align = float(record.alignments[0].hsps[0].align_length)
                cov, iden, siml = format(align*100/record.query_length,'.2f'), format(identity*100/align,'.2f'), format(positive*100/align,'.2f')
                ncbi_lin = ['noble']*8
                for i, taxa_name in enumerate(taxa_lineage[ref_id].split(';')[1:8]):
                    ncbi_lin[i] = taxa_name
                full_taxa = pd.Series(ncbi_lin,index=['domain','phylum','class','order','family','genus','species','strain'])
                selected_taxa = ";".join(list(full_taxa[:rank]))
                blast_out.append([record.query,selected_taxa,float(cov),float(iden),float(siml)])         
                taxa_assignment.write(str(record.query)+"\t"+taxa_lineage[ref_id]+"\t"+str(cov)+"\t"+str(iden)+"\t"+str(siml)+"\t"+str(ref_id)+"\n") 
            else: 
                taxa_assignment.write(str(record.query)+NO_HIT)
        df = pd.DataFrame(blast_out,columns=['OTU_name','NCBI_taxonomy','coverage(%)','identity(%)','similarity(%)'])
    print '\nDone\n'
    return df, os.path.join(o_dir,base+"taxa_assignment.txt")

def plot_histogram(df):
    fig = plt.figure(figsize=(15,7),dpi=300)
    gs = grd.GridSpec(1,2,wspace=0.2,left=0.05, right=0.95)
#   histogram1 showing alignment coverage
    ax0 = plt.subplot(gs[0])
    n,bins,patches = plt.hist(df['coverage(%)'],bins=10,range=(0,100),rwidth=0.4,edgecolor='white', color='m')
    plt.xlim(0,100)
    plt.xticks(np.arange(0,110,10,dtype=int))
    plt.ylim(0,max(n)*1.01)
    plt.xlabel('Query coverage(%) in alignment')
    plt.ylabel('Number of alignments')
#   Histogram2 showing identity and similarity of alignments
    ax1 = plt.subplot(gs[1])
    n,bins,patches = plt.hist([df['identity(%)'],df['similarity(%)']],bins=10,range=(0,100),rwidth=0.9,color={'c':(0.0,0.75),'y':(0.75,0.75,0.0)},edgecolor='white',\
                              label=['identity','similarity'])
    lg = plt.legend(fontsize=12)
    lg.draw_frame(False)
    max_count = max(max(n[0]),max(n[1]))
    plt.xlim(0,100)
    plt.xticks(np.arange(0,110,10,dtype=int))
    plt.ylim(0,max_count*1.01)
    plt.xlabel('Identity/Similarity(%) in local alignments')
    plt.ylabel('Number of alignments')
    plt.savefig(os.path.join(o_dir,'align_prof.pdf'))

def taxonomy_assignment(df, similarity, coverage, qt):
    print '\n   ...Taxonomy is assigned at {0} level...\n'.format(rank)
    siml_table = glob.glob(os.path.join(os.path.abspath(similarity),'*similarity_to_16S.txt'))[0]
    taxa_iden = os.path.join(o_dir,'taxa_at_{0}_iden.txt'.format(rank))
    taxa_siml = os.path.join(o_dir,'taxa_at_{0}_siml.txt'.format(rank))
    with open(siml_table,"r") as correspondence, open(taxa_iden,"w") as taid, open(taxa_siml,"w") as tasi:
        threshold = pd.read_csv(correspondence, delimiter='\t',index_col='Unnamed: 0')
        iden_cut, siml_cut = float(threshold[rank][1]),float(threshold[rank][2])
        df_iden = df[(df['coverage(%)'] > float(coverage)) & (df['identity(%)'] > iden_cut)]
        df_siml = df[(df['coverage(%)'] > float(coverage)) & (df['similarity(%)'] > siml_cut)]
        print 'No of OTUs taxonomically classified using identiy:     {0}'.format(df_iden['OTU_name'].count())
        print 'No of OTUs taxonomically classified using similarity:  {0}'.format(df_siml['OTU_name'].count())
        df_iden.to_csv(taid, header=False, index=False, sep='\t')
        df_siml.to_csv(tasi, header=False, index=False, sep='\t')
        
def MSA_muscle():
    base = os.path.splitext(rep_seqs)[0]
    cline = ["muscle","-in",rep_seqs,"-out",base+"_MSA.aln","-quiet"]
    MSA = subprocess.Popen(cline)
    subprocess.Popen.wait(MSA)
    return base+"_MSA.aln"
    
def phylo_FastTree(msa,method):
    ofile = os.path.join(o_dir,"rep_phyltree.tre")
    print '\n   ...Tree building using FastTree...\n'
    with open(ofile,'w') as tree:
        try:
            if method == "fasttree":
                cline = ["FastTree","-quiet","-nopr",msa]
            elif method == "fasttreemp":
                cline = ["FastTreeMP","-quiet","-nopr",msa]
        except ValueError:
            print 'Not a valid method'
        phylo = subprocess.Popen(cline,stdout=tree)
        subprocess.Popen.wait(phylo)
        print '\nDone\n'   
    
    
if __name__ == "__main__": 
    parser = ArgumentParser(description='''Taxonomy of representative protein sequences (1_otu_clustering.py) will be assigned 
                            to the best homologs of reference sequences retrieved from NCBI protein database (2_ref_DB.py)''')
    parser.add_argument('-g',dest='otus', required=True, help='the output directory of the OTU_clustering.py')
    parser.add_argument('-r',dest='ref', required=True, help='the output directory of  the reference_db.py')
    parser.add_argument('-c',dest='corr', required=True, help='the output directory of the correlation_plot.py')
    parser.add_argument('-o',dest='o_dir',required=True, help='the directory where output files will be saved')
    parser.add_argument('--num_threads',dest='num_threads',default=1, help='Number of threads (CPUs) to use in the BLAST search. Default: 2')
    parser.add_argument('--evalue',dest='evalue',default=10,help='Expectation value (E) threshold for saving hits. Default: 10')
    parser.add_argument('--taxa',dest='taxa',default='species',help='selected taxonomic rank for reliable assignment. Default: species')
    parser.add_argument('--cov',dest='coverage',default=30,help='minum query coverage in sequence alignments. Default: 30')
    parser.add_argument('--phylo',dest='method',default='fasttree',choices=['fasttree','fasttreemp'],help='FastTree methods for tree building. Default: FastTree' )
    args = parser.parse_args() 
    global o_dir, rep_seqs, otus, rank
    o_dir = make_sure_path_exists(args.o_dir)
    rank = args.taxa 
    otus = glob.glob(os.path.join(os.path.abspath(args.otus),'OTU_mapping.txt'))[0]
    rep_seqs = glob.glob(os.path.join(os.path.abspath(args.otus),'rep_seqs.fasta'))[0]
    db = makeBLASTdb(args.ref)
    xml = runBLASTP(db, args.num_threads, args.evalue)
    df,taxa = parsing_BLAST(xml, args.ref)
    plot_histogram(df)
    qt = otu_table()
    iden_taxa, siml_taxa = taxonomy_assignment(df, args.corr, args.coverage, qt)
    otu_taxa(qt,taxa,iden_taxa,siml_taxa)
    msa = MSA_muscle()
    phylo_FastTree(msa,args.method)
