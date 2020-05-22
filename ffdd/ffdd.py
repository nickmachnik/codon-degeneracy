#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import unambiguous_dna_by_name
from skbio.alignment import local_pairwise_align_ssw
from skbio.sequence import Protein
from skbio.alignment._pairwise import blosum50


def _truncate(seq: str):
    """Removes the prefix before the first ATG.

    Args:
        seq (str): Nucleotide sequence to truncate

    Returns:
        str or None: the truncated sequence or None if
                     no start codon was found.
    """
    seq = seq.upper()
    for start in range(0, len(seq) - 2):
        codon = seq[start:start+3]
        if codon == "ATG":
            return seq[start:]


def _translate(seq: str, table="Standard"):
    """Translates a given DNA sequence into a
    protein sequence, using a specified codon table.

    Args:
        seq (str): DNA sequence. Expects the coding strand and
                    a start with a leading ATG codon.
        table (str, optional): NCBI table name as used in Bio.Data.CodonTable

    Returns:
        skbio.sequence.Protein: Protein object with the translated sequence.
    """
    return Protein(str(Seq(seq, IUPAC.unambiguous_dna).translate(table=table)))


def _align(a: Protein, b: Protein):
    """Wraps the skbio pairwise ssw alilgner.

    Args:
        a (str): sequence a
        b (str): sequence b

    Returns:
        skbio.alignment.TabularMSA: skbio alignment table
    """
    return local_pairwise_align_ssw(a, b, substitution_matrix=blosum50)


def _x_fold_degenerate_from_codon_table(x: int, table: str):
    """Extracts amino acids that are encoded by x
    different codons from an NCBI codon table.

    Args:
        x (int): fold degeneration to filter by.
        table (str, optional): NCBI table name as used in Bio.Data.CodonTable

    Returns:
        dict: amino acids (keys) and their x codons (values)
    """
    codon_table = unambiguous_dna_by_name[str].forward_table
    reverse_table = {}
    for codon, aa in codon_table.items():
        reverse_table.setdefault(aa, [])
        reverse_table[aa].append[codon]
    x_fold_table = {}
    for aa, codons in reverse_table.items():
        if len(codons) == x:
            x_fold_table[aa] = codons
    return x_fold_table


def _triplet(seq: str, i: int):
    """Return the ith triplet in a sequence.

    Args:
        seq (str): Sequence
        i (int): 0-based index of the target triplet.

    Returns:
        str: The target triplet.
    """
    start = i * 3
    return seq[start:start+3]


def _get_aligned_ffds(a: str, b: str, table_a="Standard", table_b="Standard"):
    proteins = []
    for s, table in zip([a, table_a], [b, table_b]):
        ts = _truncate(s)
        if not ts:
            raise ValueError("DNA sequence without ATG codon provided!")
        else:
            proteins.append(_translate(ts, table))
    alignment = _align(*proteins)
    unaligned_a_start = alignment[2][0][0]
    unaligned_b_start = alignment[2][1][0]

