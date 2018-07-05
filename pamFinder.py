#!/usr/bin/env python

import sys, argparse

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

        self.find_PAM()
        Region.regions.append(self)

    def offsetSequence(self):
        """return the full sequence"""
        return Region.chromDict[self.chrom][self.offsetStart: self.offsetEnd]
    def find_PAM(self):
        seq = self.offsetSequence()
        # for (a, b) in zip(seq[0::2], seq[1::2]):
        #     print a, b
        # for a in seq:
        #     print a




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






if __name__ == '__main__':
    main()
