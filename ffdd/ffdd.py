#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import unambiguous_dna_by_name
from skbio.alignment import local_pairwise_align_ssw
from skbio.sequence import Protein
from skbio.alignment._pairwise import blosum50
from itertools import combinations


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


def _hamming_distance(a, b):
    distance = 0
    differential_sites = []
    for p in range(len(a)):
        if a[p] != b[p]:
            distance += 1
            differential_sites.append(p)
    return distance, differential_sites


def _site_degeneracy(codons):
    """Group codons by site for which they are degenerate.

    Args:
        codons (array_like): array of base triplets, e.g. codons for 
                             a particular amino acid

    Returns:
        dict: mapping between each degenerate site and the groups of 
              codons which are degenerate at that site.
    """
    def resolve_pairs(pairs):
        codons = set()
        for a, b in pairs:
            codons.add(a)
            codons.add(b)

        resolved = [set([codons.pop()])]
        for a in codons:
            added = False
            for group in resolved:
                if any((a, b) in pairs or (b, a) in pairs for b in group):
                    group.add(a)
                    added = True
                    break
            if not added:
                resolved.append(set([a]))
        return resolved

    pairs2sites = {}
    sites2pairs = {}
    for a, b in combinations(codons, 2):
        distance, sites = _hamming_distance(a, b)
        if distance == 1:
            pairs2sites[(a, b)] = sites[0]
            sites2pairs.setdefault(sites[0], set())
            sites2pairs[sites[0]].add((a, b))

    return {site: resolve_pairs(pairs)
            for site, pairs in sites2pairs.items()}


def _x_fold_degenerate_aa_from_codon_table(x: int, table: str):
    """Extracts amino acids that are encoded by x
    different codons from an NCBI codon table.

    Args:
        x (int): fold degeneration to filter by.
        table (str, optional): NCBI table name as used in Bio.Data.CodonTable

    Returns:
        dict: amino acids (keys) and their x codons (values)
    """
    codon_table = unambiguous_dna_by_name[table].forward_table
    reverse_table = {}
    for codon, aa in codon_table.items():
        reverse_table.setdefault(aa, [])
        reverse_table[aa].append(codon)
    x_fold_table = {}
    for aa, codons in reverse_table.items():
        if len(codons) >= x:
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


def _aligned_ffds(a: str, b: str, table_a="Standard", table_b="Standard"):
    """Extracts the four-fold degenerate codons at conserved sites
    from a protein sequence alignment of two coding DNA sequences.

    Args:
        a (str): coding DNA sequence a
        b (str): coding DNA sequence b
        table_a (str, optional): NCBI table name as used in Bio.Data.CodonTable
        table_b (str, optional): NCBI table name as used in Bio.Data.CodonTable

    Yields:
        tuple: aligned four-fold degenerate codons in sequence a and b

    Raises:
        ValueError: if any of the input sequences do not contain
                    an ATG start codon.
    """
    proteins = []
    truncated = []
    for s, table in [[a, table_a], [b, table_b]]:
        ts = _truncate(s)
        if not ts:
            raise ValueError("DNA sequence without ATG codon provided!")
        else:
            truncated.append(ts)
            proteins.append(_translate(ts, table))
    alignment = _align(*proteins)

    # shorten the input sequences to the aligned subsequences
    a = truncated[0][alignment[2][0][0]*3:]
    b = truncated[1][alignment[2][1][0]*3:]

    degenerate_aa = set(
        _x_fold_degenerate_aa_from_codon_table(4, table_a)).intersection(
            set(_x_fold_degenerate_aa_from_codon_table(4, table_b)))

    # iterate over the aligned sequences
    for i, (ca, cb) in enumerate(zip(
        str(alignment[0][0]), str(alignment[0][1]))
    ):
        if ca == cb and ca in degenerate_aa:
            yield _triplet(a, i), _triplet(b, i)


def substitution_rate_at_ffds(
    a: str, b: str, table_a="Standard", table_b="Standard"
) -> ((int, int), [str, str]):
    """Estimates the rate of neutral substitutions by counting
    the number of substitutions at four-fold degenerate sites.
    sites.

    Args:
        a (str): coding DNA sequence a
        b (str): coding DNA sequence b
        table_a (str, optional): NCBI table name as used in Bio.Data.CodonTable
        table_b (str, optional): NCBI table name as used in Bio.Data.CodonTable

    Returns:
        (int, int): number of substitutions, number of sites
        [str, str]: the selected ORFs of the input sequences
    """
    proteins = []
    truncated = []
    for s, table in [[a, table_a], [b, table_b]]:
        ts = _truncate(s)
        if not ts:
            raise ValueError("DNA sequence without ATG codon provided!")
        else:
            truncated.append(ts)
            proteins.append(_translate(ts, table))
    alignment = _align(*proteins)

    # shorten the input sequences to the aligned subsequences
    a = truncated[0][alignment[2][0][0]*3:]
    b = truncated[1][alignment[2][1][0]*3:]

    degenerate_aa = set(
        _x_fold_degenerate_aa_from_codon_table(4, table_a)).intersection(
            set(_x_fold_degenerate_aa_from_codon_table(4, table_b)))

    n_sites = 0
    n_sub = 0
    # iterate over the aligned sequences
    for i, (ca, cb) in enumerate(zip(
        str(alignment[0][0]), str(alignment[0][1]))
    ):
        if ca == cb and ca in degenerate_aa:
            tra, trb = _triplet(a, i), _triplet(b, i)
            n_sites += 1
            if tra[2] != trb[2]:
                n_sub += 1

    return (n_sub, n_sites), truncated
