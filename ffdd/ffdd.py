#!/usr/bin/env python


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