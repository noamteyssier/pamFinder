#!/usr/bin/env python


"""
A sketchpad to use w/ jupyter
"""

import os, sys, itertools
import numpy as np
import pandas as pd
from itertools import product


def zipList(a,b, delim='-'):
    """function to zip two lists by a delim"""
    if len(a) == len(b):
        return [delim.join([str(a[i]),str(b[i])]) for i in range(len(a))]
    else:
        sys.exit("zipList requires iterators of same length")



lPos = [7, 24, 35, 48, 64, 138, 150, 162, 177]
lID = ['CC', 'GG','CC', 'GG','CC', 'GG','CC', 'GG', 'GG']
rPos = [16, 66, 86, 87, 103, 136, 150, 167, 207]
rID = ['CC', 'GG','CC', 'GG','CC', 'GG','CC', 'GG', 'GG']


start, end, oStart, oEnd = 114472, 114559, 114259, 114772

trueLeft = [l + oStart for l in lPos]
trueRight = [r + end for r in rPos]

# left = zipList(trueLeft, lID)
# right = zipList(trueRight, rID)
left = zipList(lPos, lID)
right = zipList(rPos, rID)


## convert to df
combinations = pd.DataFrame(list(product(left, right)), columns = ['left','right'])
leftSide = combinations['left'].str.split('-', expand=True).rename(columns = {0 : 'left', 1 : 'lPAM'})
rightSide = combinations['right'].str.split('-', expand=True).rename(columns = {0 : 'right', 1 : 'rPAM'})

combinations = leftSide.join(rightSide, how = 'outer')
combinations['left'] = pd.to_numeric(combinations['left'])
combinations['right'] = pd.to_numeric(combinations['right'])

combinations['left'] = np.where(combinations['lID' == 'CC', ])


## calculate size of interval
combinations['size'] = combinations['right'] - combinations['left']
## return intervals less than 300
intervals = combinations.loc[combinations['size'] <= 300].reset_index(drop=True)

leftSide['side'] = 'l'
rightSide['side'] = 'r'
l = leftSide.rename(columns = {'left' : 'pos', 'lPAM' : 'PAM'})
r = rightSide.rename(columns = {'right' : 'pos', 'rPAM' : 'PAM'})

l.append(r)
