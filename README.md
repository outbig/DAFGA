DAFGA
=====

Diversity analysis of Functional Gene Amplicons

Diversity analysis of functional marker genes provides physiological insights into microbial guilds that perform an ecologically relevant process. However, it is challenging to group functional gene sequences to valid taxonomic units, primarily due to varying evolutionary rates of the individual gene sequences and possible horizontal gene transfer events. We developed a python script package named DAFGA, which estimates the evolutionary rate of a particular functional gene in a standardized manner by relating its sequence divergence to that of the 16S rRNA gene. As a result, DAFGA provides gene-specific parameter sets for OTU clustering and taxonomic assignment at desired rank, and it can be implemented into the diversity measurements offered by QIIME.


Dependencies
------------

- python:	The toolkit is primarily written in python 2.7.2. 
- biopython:	A set of freely available tools for biological computation (Download from http://biopython.org/wiki/Download)	
- numpy, scipy, pandas, matplotlib:	Python modules used in the toolkit
- usearch:	Search and clustering algorithms used for OTU clustering
- emboss:	Smith-Waterman algorithm is used to calculate the local alignment of two sequences (Download from http://emboss.sourceforge.net/download/)
- FastTree:	Construction of approximately-maximum-likelihood phylogenetic trees from protein sequences alignment (Download from http://www.microbesonline.org/fasttree/#Install)
- Muscle:	Multiple sequence alignment software for protein sequences (Download from http://www.drive5.com/muscle/downloads.htm)
- BLAST+:	A new suite of BLAST tools that utilizes the NCBI C++ Toolkit (Download from  ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) 

INSTALL
-------
- DAFGA is developed and tested under linux system. The package installation and function works in linux terminal. 
- Download DAFGA scripts from https://github.com/outbig/DAFGA. You can find a tab of "Download ZIP" on the right panel
- > unzip DAFGA-mater.zip 
- > cd v1.0 (Go to the folder of the lastest version of DAFGA scripts)
- > sudo python setup.py install


QUESTIONS & BUG REPORTS
-----------------------

If you have questions or find a bug please report it at
https://github.com/outbig/DAFGA/issues?state=open


COPYRIGHT AND LICENSE
---------------------

Copyrights (C) 2014 by Yongkyu Kim All rights reserved

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.
