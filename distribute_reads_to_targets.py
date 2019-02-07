#!/usr/bin/env python

from __future__ import division     # for correct progress bar output in both Python 2 and 3

import sys,os,errno
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

"""
This script is part of a pipeline to extract phylogenetically-useful sequences from
Illumina data using the targeted (liquid-phase) sequence enrichment approach.

After a BLASTx search of the raw reads against the target sequences, the reads need to be
sorted according to the successful hits. This script takes the BLASTx output (tabular)
and the raw read files, and distributes the reads into FASTA files ready for assembly.

If there are multiple BLAST results (for example, one for each read direction),
concatenate them prior to sorting.
"""

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

class ProgressIndicator:
    """Class for output of ASCII progress indicator"""

    def __init__(self, max, width):
        self.max = max
        self.width = width
        self.current = None
        self.prev_n = None

    def show(self, value):
        maximum = self.max
        n = int(self.width * value/maximum)
        if n != self.prev_n:
            sys.stderr.write('\r[' +  n*'=' + (self.width - n)*' ' + ']')
            sys.stderr.flush()
            self.prev_n = n

    def set(self, value):
        self.prev_n = None      # force redraw
        self.show(value)

    def finish(self):
        sys.stderr.write("\n")

def read_sorting(blastfilename):
    read_hit_dict = {}
    blastfile = open(blastfilename)
    for line in blastfile:
        line = line.split()
        readID = line[0]
        target = line[1].split("-")[-1]
        if readID in read_hit_dict:
            if target not in read_hit_dict[readID]:
                read_hit_dict[readID].append(target)
        else:
            read_hit_dict[readID] = [target]
    return read_hit_dict

def write_paired_seqs(target,ID1,Seq1,ID2,Seq2,single=True):
    mkdir_p(target)
    if single:
        outfile = open(os.path.join(target,"{}_interleaved.fasta".format(target)),'a')
        outfile.write(">{}\n{}\n".format(ID1,Seq1))
        outfile.write(">{}\n{}\n".format(ID2,Seq2))
        outfile.close()
    else:
        outfile1 = open(os.path.join(target,"{}_1.fasta".format(target)),'a')
        outfile1.write(">{}\n{}\n".format(ID1,Seq1))
        outfile2 = open(os.path.join(target,"{}_2.fasta".format(target)),'a')
        outfile2.write(">{}\n{}\n".format(ID2,Seq2))
        outfile1.close()
        outfile2.close()

def write_single_seqs(target,ID1,Seq1):
    """Distributing targets from single-end sequencing"""
    mkdir_p(target)
    outfile = open(os.path.join(target,"{}_unpaired.fasta".format(target)),'a')
    outfile.write(">{}\n{}\n".format(ID1,Seq1))
    outfile.close()

def distribute_reads(readfiles, read_hit_dict, single=True):
    if len(readfiles) != 2:
        print("\nError: distribute_reads expected len(readfiles)=2")
        sys.exit(1)

    num_reads_to_write = len(read_hit_dict)
    iterator1 = FastqGeneralIterator(open(readfiles[0]))
    iterator2 = FastqGeneralIterator(open(readfiles[1]))

    reads_written = 0
    sys.stderr.write("Read distributing progress:\n")

    prog = ProgressIndicator(num_reads_to_write, 40)
    for ID1_long, Seq1, Qual1 in iterator1:

        ID2_long, Seq2, Qual2 = next(iterator2)

        ID1 = ID1_long.split()[0]
        if ID1.endswith("/1") or ID1.endswith("/2"):
            ID1 = ID1[:-2]

        ID2 = ID2_long.split()[0]
        if ID2.endswith("/1") or ID2.endswith("/2"):
            ID2 = ID2[:-2]

        if ID1 in read_hit_dict:
            for target in read_hit_dict[ID1]:
                write_paired_seqs(target,ID1,Seq1,ID2,Seq2)
            reads_written += 1
            prog.show(reads_written)
        elif ID2 in read_hit_dict:
            for target in read_hit_dict[ID2]:
                write_paired_seqs(target,ID1,Seq1,ID2,Seq2)
            reads_written += 1
            prog.show(reads_written)
    prog.finish()

def main():
    blastfilename = sys.argv[1]
    readfiles = sys.argv[2:]

    print("\nDistributing reads")
    print("  blastfilename = {}".format(blastfilename))
    print("  readfiles     = {}".format("\n                  ".join(readfiles)))
    read_hit_dict = read_sorting(blastfilename)
    print("\nUnique reads with hits: {}".format(len(read_hit_dict)))
    distribute_reads(readfiles, read_hit_dict, single=True)

if __name__ == "__main__":
	main()
