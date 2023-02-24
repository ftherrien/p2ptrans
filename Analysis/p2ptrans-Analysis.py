#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 12:22:01 2022

@author: tim
"""

import pandas as pd
r1=pd.read_csv('results1.csv', sep=' ', header=None)
r2=pd.read_csv('results2.csv',sep='    ', header=None)
r3=pd.read_csv('results3.csv',sep='    ', header=None)
r4=pd.read_csv('results4.csv', header=None)
r5=pd.read_csv('results5.csv', sep=' ', header=None)
r6=pd.read_csv('results6.csv', sep=' ', header=None)
r7=pd.read_csv('results7.csv', sep=' ', header=None)
r8=pd.read_csv('results8.csv', sep=' ', header=None)
r9=pd.read_csv('results9.csv', sep=' ', header=None)

r1.columns = ['location', 'dist_old', 'dist_new', 'Difference', 'Percent']
r2.columns = ['closest hkl']
r3.columns = ['closest hkl']
r4.columns = ['Transition structures']
r5.columns = ['Mapped atoms']
r6.columns = ['Stretching factor']
r7.columns = ['TC size']
r8.columns = ['Total distance', '0', '1']
r9.columns = ['1','e1', '2','e2', '3','e3']



dataframes = [r1, r2, r3, r4, r5, r6, r7, r8['Total distance'], r9['e1'], r9['e2'], r9['e3']]
dataframes = pd.concat(dataframes, axis=1)

dataframes.to_csv('Results.csv', index=False, sep='\t')
