#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
# Author:   Yongkyu Kim, PhD
            Max Planck Institute for terrestrial microbiology
# 
# Date:     2014-04-02
# Version:  v1.0
# 
# Try 'dafga_correlation.py -h' for more information and See manual
# 
# 
# Bugs: Please report to https://github.com/outbig/DAFGA/issues?state=open
"""

from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Emboss.Applications import WaterCommandline
import os, glob, linecache
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import stats

def pairing_fungene_with_rRNA(rRNA_seqs, fungene_seqs):
    base_fg, base_ssu = os.path.splitext(fungene_seqs)[0], os.path.splitext(rRNA_seqs)[0]
    fg, ssu = {}, {}
    with open(fungene_seqs,"r") as fungene, open(rRNA_seqs,"r") as rRNA, open(base_fg+"_pair.fasta","w") as fungene_pair, open(base_ssu+"_pair.fasta","w") as rRNA_pair:             
        for record_fg in SeqIO.parse(fungene,"fasta"):
            fg[record_fg.id] = record_fg
        for record_ssu in SeqIO.parse(rRNA,"fasta"):
            ssu[record_ssu.description.split("\t")[1]] = record_ssu
        for seqs in set(fg.keys())&set(ssu.keys()):
            SeqIO.write(fg[seqs],fungene_pair,"fasta")
            SeqIO.write(ssu[seqs],rRNA_pair,"fasta")
    return base_fg+"_pair.fasta", base_ssu+"_pair.fasta"

def input_files(in_dir):
    fg = glob.glob(os.path.join(in_dir,'*_strain.fasta'))[0]
    ssu = glob.glob(os.path.join(in_dir,'*_strain_16S_rRNAs.fasta'))[0]
    return fg, ssu

def batch_iterator(iterator, batch_size) :
    entry = True
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                break
            batch.append(entry)
        if batch :
            yield batch
       
def make_directory(dirname):
    try:
        os.makedirs(dirname)
    except OSError:
        raise
    return os.path.abspath(dirname)
    
def split_into_pieces(fasta, seq_type):
    split_dir = make_directory(os.path.join(o_dir,seq_type+"_split"))
    record_iter = SeqIO.parse(open(fasta,"r"),"fasta")
    for i, batch in enumerate(batch_iterator(record_iter,1)):
        filename = split_dir+"/query_%i.fa" % (i+1)
        handle = open(filename,"w")
        SeqIO.write(batch, handle, "fasta")
        handle.close()
    return split_dir
    
def emboss_local_pairwise_alignment(query_dir, seq_type):
    if seq_type == 'fg':
        print '\n   ...pairwise comparison of functional gene sequences...\n' 
    elif seq_type == 'ssu':
        print '\n   ...pairwise comparison of SSU rRNA sequences...\n'
    water_cline = WaterCommandline()
    water_cline.gapopen=10
    water_cline.gapextend=0.5
    query_list = [query for query in sorted(glob.glob(query_dir+"/*.fa"))]
    for i, a_seq in enumerate(query_list): 
        water_cline.asequence=str(a_seq)
        for j, b_seq in enumerate(query_list[i:]):
            water_cline.bsequence=str(b_seq)
            align_out = query_dir+"/pairwise_"+str(i+1)+"_"+str(i+j+1)+".aln"
            water_cline.outfile=str(align_out)
            water_cline()
    print 'Done\n'
    return query_dir+"/*.aln"
            
def identity_similarity(alignments):
    identity, similarity, aln_out = [], [], []
    for alignment in sorted(glob.glob(alignments)):
        iden, siml = linecache.getline(alignment,24), linecache.getline(alignment,25)
        aln_out.append([float(iden.split("(")[1].split("%")[0]),float(siml.split("(")[1].split("%")[0])])
    df = pd.DataFrame(aln_out, columns=('identity(%)','similarity(%)'))
    return df

def intercept(b):
    if b > 0:
        b = "+%.2f" % abs(b)
    elif b < 0:
        b = "-%.2f" % abs(b)
    else:
        b = " "
    return b
 
def scatter_plot(ssu_df, fg_df):
    ssu_iden, fg_iden, fg_siml = ssu_df['identity(%)'], fg_df['identity(%)'], fg_df['similarity(%)']
    fig = plt.figure(figsize=(15,7),dpi=300)
    gs = gridspec.GridSpec(1,2,wspace=0.2,left=0.05, right=0.95)
#   correlation plot of 16S rRNA identity versus funtional gene identity
    ax0 = plt.subplot(gs[0])
    plt.scatter(ssu_iden,fg_iden,color='blue',s=1)
    iden_func = stats.linregress(ssu_iden,fg_iden)
    x_rg = range(int(min(ssu_iden)),int(max(ssu_iden))+1)
    y_rg = np.polyval([iden_func[0],iden_func[1]],x_rg)
    plt.text(5,95, r'$y =  %.2f x  %s $' % (iden_func[0],intercept(iden_func[1])), fontsize=15)
    plt.text(5,90, r'$R^2=%.4f$' % (iden_func[2]**2))
    plt.text(5,85, r'$P-value=%.2e$' % (iden_func[3]))
    plt.text(5,80, r'$StdErr=%.4f$' % (iden_func[4]))
    plt.title('16S rRNA identity vs. Funtional gene identity')
    plt.plot(x_rg,y_rg,'r--',label='line 1')   
    plt.xlabel('16S rRNA gene identity (%)')
    plt.ylabel('Funtional gene identity (%)')
    plt.ylim(0,100)
    plt.xlim(0,100)
#   correlation plot of 16S rRNA identity versus funtional gene similarity
    ax1 = plt.subplot(gs[1])
    plt.scatter(ssu_iden,fg_siml,color='green',s=1)
    siml_func = stats.linregress(ssu_iden,fg_siml)
    x_rg = range(int(min(ssu_iden)),int(max(ssu_iden))+1)
    y_rg = np.polyval([siml_func[0],siml_func[1]],x_rg)
    plt.text(5,95, r'$y =  %.2f x  %s $' % (siml_func[0],intercept(siml_func[1])), fontsize=15)
    plt.text(5,90, r'$R^2=%.4f$' % (siml_func[2]**2))
    plt.text(5,85, r'$P-value=%.2e$' % (siml_func[3]))
    plt.text(5,80, r'$StdErr=%.4f$' % (siml_func[4]))
    plt.title('16S rRNA identity vs. Funtional gene similarity')
    (m,b) = np.polyfit(ssu_iden,fg_siml, 1)
    x_rg = range(int(min(ssu_iden)),int(max(ssu_iden))+1)
    y_rg = np.polyval([m,b],x_rg)
    plt.plot(x_rg,y_rg,'r--')
    plt.xlabel('16S rRNA gene identity (%)')
    plt.ylabel('Funtional gene similarity (%)')
    plt.ylim(0,100)
    plt.xlim(0,100)
    plt.savefig(o_dir+'/correlation_plot.pdf')
    return iden_func, siml_func

def corresponding_similarity(iden_eq,siml_eq):
    base = os.path.splitext(o_dir)[0]
    rRNA, fg_iden, fg_siml = [99,97,95,93,90,85,80], [], []
    for item in rRNA:
        fg_iden.append(np.polyval([iden_eq[0],iden_eq[1]],item))
        fg_siml.append(np.polyval([siml_eq[0],siml_eq[1]],item))
    df_corr = pd.DataFrame([rRNA,fg_iden,fg_siml],dtype=int,columns=['strain','species','genus','family','order','class','phylum'],index=['16S_rRNA_iden(%)','Functional_gene_iden(%)','Functional_gene_siml(%)'])
    with open(base+"/corresponding_similarity_to_16S.txt","w") as correspondence:
        df_corr.to_csv(correspondence,sep="\t")
    print df_corr
    
if __name__ == "__main__":
    parser =  ArgumentParser()
    parser.add_argument('-r', dest='ref', required=True, help='the output directory of reference_db.py which includes *_strain.fasta and *_strain_SSU_rRNA.fasta')    
    parser.add_argument('-o', dest='o_dir', required=True, help='the directory where output files will be saved')
    args = parser.parse_args()
    global o_dir
    o_dir = os.path.abspath(args.o_dir)
    fg, ssu = input_files(args.ref) 
    fg_dir = split_into_pieces(fg,"fg")
    ssu_dir = split_into_pieces(ssu,"ssu")
    fg_aln = emboss_local_pairwise_alignment(fg_dir,"fg")
    ssu_aln = emboss_local_pairwise_alignment(ssu_dir,"ssu")
    ssu_df = identity_similarity(ssu_aln)   
    fg_df = identity_similarity(fg_aln)
    iden_func, siml_func = scatter_plot(ssu_df, fg_df)
    corresponding_similarity(iden_func,siml_func)
