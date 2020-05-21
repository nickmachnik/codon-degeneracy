#!/usr/bin/env python

import kb


def truncate(seq: str):
    """Removes the prefix before the first ATG.

    Args:
        seq (str): Nucleotide sequence to truncate

    Returns:
        str or None: the truncated sequence or None if
                     no start codon was found.
    """
    seq = seq.upper()
    for start in range(0, len(seq) - 2, 3):
        codon = seq[start:start+3]
        if codon == "ATG":
            return seq[start:]


def transcribe(seq: str):
    """Transcribe a DNA sequence into a protein.

    Args:
        seq (str): DNA sequence

    Returns:
        str: Protein sequence
    """
    seq = seq.upper()
    translated = ""
    for start in range(0, len(seq) - 2, 3):
        codon = seq[start:start+3]
        try:
            aa = kb.codon_table[codon]
        except KeyError:
            print("Codon {} is not in the codon table!".format(codon))
            exit(1)
        if aa == '*':
            return seq
        else:
            translated += aa
