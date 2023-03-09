#!/usr/bin/env python

# Author: Justine Macalindong justinem@uoregon.edu


'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''
__version__ = "0.6"         
                            
import re

DNA_bases = set('ATGCNatcgn')
RNA_bases = set('AUGCNaucgn')

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter) - 33

def qual_score(phred_score: str) -> float:
    """Calculates the average quality score of an entire string"""
    sum = 0 
    length = 0 
    for letter in phred_score:
        length += 1
        sum += convert_phred(letter)
    return(sum/length)

def validate_base_seq(seq,RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA: str) -> float:
    '''This function calculates the GC content of a DNA sequence.'''
    DNA = DNA.upper()
    return (DNA.count("G") + DNA.count("C")) / len(DNA)

def oneline_fasta(file):
    '''This function takes all sequence lines of a FASTA file and puts them all on one line.'''
    # make dict with headers as keys and sequences as values
    seq_dict = {}
    with open(file, 'r') as fa:
        line_count = 0
        for line in fa:
            line_count +=1
            line = line.strip('\n')
            # only get header lines
            if line[0] == '>':
                header_line = line
            # populate dict with seq lines (non-header lines)
            else:
                if header_line not in seq_dict:
                    seq_dict[header_line] = line
                else:
                    seq_dict[header_line] += line
    # write out to file
    fa_one_line = open('one_line.fa', 'w')
    for keys,vals in seq_dict.items():
        fa_one_line.write(str(keys) + '\n' + str(vals) + '\n')
    fa_one_line.close()
    # return len(seq_dict)

def reverse_complement(seq: str) -> str:
    '''This function takes a sequence and returns the reverse complement of it.'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement[i] for i in reversed(seq))

def strand_flag(flag: str) -> str:
    '''This function takes a a bitwise flag and returns its strandedness.'''
    if ((flag & 16) == 16):
	    return "-"
    else:
        return "+"

def position_adjust(pos: int,cigar: str,strand: str) -> int:
    '''This function adjusts position based on cigar string and strandedness.'''
    if strand == "+":
        soft = re.findall(r'^([0-9]+)S', cigar)
        if soft:
            return pos - int(soft[0])
        else:
            return pos
    if strand == "-":
        deletion = 0
        skip = 0
        match = 0
        soft = re.findall(r'([0-9]+)S$', cigar)
        D = re.findall(r'([0-9]+)D', cigar)
        N = re.findall(r'([0-9]+)N', cigar)
        M = re.findall(r'([0-9]+)M', cigar)
        for x in D:
            deletion += int(x)
        for y in N:
            skip += int(y)
        for z in M:
            match += int(z)
        if soft:
            total = pos + int(soft[0]) + deletion + skip + match
            return total
        else:
            total = pos + deletion + skip + match
            return total
        

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calcluated the correct average phred score")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATGCAT") == 0.5
    print("correctly calculated GC content")