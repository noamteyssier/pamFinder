#!/usr/bin/env python

import sys, argparse
import pandas as pd
from itertools import product

class Region:
    regions = list()
    chromDict = dict()
    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.size = self.end - self.start
        self.offset = 300 - self.size # range to look on either side of region

        self.offsetStart = self.start - self.offset
        self.offsetEnd = self.end + self.offset

        self.posInfo = {'l' : dict(), 'r' : dict() } # {side : CC/GG : index}
        self.leftPositions = self.find_PAM(self.offsetSequence(l = True), 'l')
        self.rightPositions = self.find_PAM(self.offsetSequence(r = True), 'r')

        # combinations of left and right gRNA positions where size < 300
        self.intervals = None

        self.find_intervals()
        self.pull_sequences()
        sys.exit()


        Region.regions.append(self)
    def add_position(self, side, pam, index):
        if pam not in self.posInfo[side]:
            self.posInfo[side][pam] = list()
        self.posInfo[side][pam].append(index)
    def pull_positions(self, side, pam):
        """return list of positions for side and pam"""
        try:
            return ['+'.join([str(s), pam]) for s in self.posInfo[side][pam]]
        except KeyError:
            return []
    def offsetSequence(self, l=False, r=False):
        """return the full sequence"""
        if (l == False) and (r == False):
            return Region.chromDict[self.chrom][self.offsetStart: self.offsetEnd]
        elif (l == True) and (r == False):
            return Region.chromDict[self.chrom][self.offsetStart : self.start]
        elif (l == False) and (r == True):
            return Region.chromDict[self.chrom][self.end : self.offsetEnd]
        else:
            return [Region.chromDict[self.chrom][self.offsetStart : self.start], Region.chromDict[self.chrom][self.end : self.offsetEnd]]
    def find_PAM(self, sequence, side):
        """method to find positions of GG or CC in sequence"""
        ggcc = list()
        for i in range(len(sequence)):
            try:
                # search through sequence two positions at a time to search for GG or CC
                twoMer = ''.join([sequence[i], sequence[i+1]])
                if twoMer == 'GG' or twoMer == 'CC':
                    # actual cut site will be 4 positions before GG
                    if twoMer == 'GG':
                        index = i - 4
                    # and 4 positions before CC
                    else:
                        index = i + 4
                    self.add_position(side, twoMer, index)
                    ggcc.append(index)
            except IndexError:
                break
        return ggcc
    def find_intervals(self):
        """find intervals with a size less than 300"""

        # pull left and right indexes with corresponding PAM
        left = self.pull_positions('l','GG') + self.pull_positions('l','CC')
        right = self.pull_positions('r','GG') + self.pull_positions('r','CC')

        # matrix multiplication of left and right for all combinations
        combinations = pd.DataFrame(list(product(left, right)), columns = ['left','right'])

        # split strings on delimeter to separate index from pam into individual cols
        leftSide = combinations['left'].str.split('+', expand=True).rename(columns = {0 : 'left', 1 : 'lPAM'})
        rightSide = combinations['right'].str.split('+', expand=True).rename(columns = {0 : 'right', 1 : 'rPAM'})

        # only perform if dataframe is not empty
        if combinations.empty == False:
            # join left and right sides
            combinations = leftSide.join(rightSide, how = 'outer')

            # balance index to chromosomal position
            combinations['left'] = pd.to_numeric(combinations['left']) + self.offsetStart
            combinations['right'] = pd.to_numeric(combinations['right']) + self.end

            ## calculate size of interval
            combinations['size'] = combinations['right'] - combinations['left']

            ## return intervals less than 300
            self.intervals = combinations.loc[combinations['size'] <= 300].reset_index(drop=True)
        else:
            self.intervals = None
    def pull_sequences(self):
        print self.intervals
        pass







def parse_fasta(genome):
    """parse each line of file"""
    with open(genome) as f:
        while True:
            try:
                yield next(f).strip('\n')
            except StopIteration:
                break
def parse_bed(bed):
    """read lines of bed and split on tabs"""
    with open(bed) as f:
        while True:
            try:
                yield next(f).strip('\n').split('\t')
            except StopIteration:
                break
def join_lines(chromDict):
    """"join lines of sequence for one long string"""
    for c in chromDict:
        chromDict[c] = ''.join(chromDict[c])
    return chromDict
def make_chrom_dict(genome_fn):
    """assign each line to its respective chromosome"""
    chromDict = dict()
    for line in parse_fasta(genome_fn):
        if '>' in line:
            currentChrom = line.strip('>').split('|')[0].strip(' ')
            chromDict[currentChrom] = list()
            continue
        chromDict[currentChrom].append(line)
    return join_lines(chromDict)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", '--input', help="input genome", required=True)
    p.add_argument('-b', '--bed', help='input bed file', required=True)
    args = p.parse_args()

    Region.chromDict = make_chrom_dict(args.input)
    [Region(b[0], b[1], b[2]) for b in parse_bed(args.bed)]


'''
statistics to gather :
- number of nGG/nCC in gRNA
- GC content of gRNA
- number of times gRNA found in Genome
'''




if __name__ == '__main__':
    main()
