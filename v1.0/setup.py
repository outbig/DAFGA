#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
# Author:   Yongkyu Kim, PhD
            Max Planck Institute for terrestrial microbiology
# 
# Date:     2014-04-02
# Version:  1.0
# 
# 
# Purpose:   
# Bugs: Please report to https://github.com/outbig/DAFGA/issues?state=open

"""

from distutils.core import setup

setup(
    name='DAFGA',
    version='1.0',
    description='Diversity analysis of functional gene amplicons',
    author='Yongkyu Kim',
    author_email='outbig@gmail.com ',
    url = "https://github.com/outbig/DAFGA",
    scripts=['dafga_refDB.py','dafga_correlation.py','dafga_otus.py','dafga_taxonomy.py']
    )
