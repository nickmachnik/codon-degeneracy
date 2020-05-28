#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Data.CodonTable import unambiguous_dna_by_name
from skbio.alignment import local_pairwise_align_ssw
from skbio.sequence import Protein
from skbio.alignment._pairwise import blosum50
from itertools import combinations


def _truncate(seq: str):
    """Removes the prefix before the first ATG.
    and any trailing nucleotides that are not
    part of a full codon.

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
            res = seq[start:]
            trailing_nuc = len(res) % 3
            if trailing_nuc > 0:
                return res[:-trailing_nuc]
            else:
                return res


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
    return Protein(str(Seq(seq).translate(table=table)))


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


def _site_degeneracy(codons, min_x=None):
    """Group codons by site for which they are degenerate.

    Args:
        codons (array_like): array of base triplets, e.g. codons for
                             a particular amino acid.
        min_x (None, optional): Filters out groups of degenerate codons that
                                have less than min_x elements.

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

    sites2pairs = {}
    for a, b in combinations(codons, 2):
        distance, sites = _hamming_distance(a, b)
        if distance == 1:
            sites2pairs.setdefault(sites[0], set())
            sites2pairs[sites[0]].add((a, b))

    if min_x is None:
        res = {site: resolve_pairs(pairs)
               for site, pairs in sites2pairs.items()}
    else:
        res = {}
        for site, pairs in sites2pairs.items():
            resolved = [s for s in resolve_pairs(pairs) if len(s) >= min_x]
            if resolved:
                res[site] = resolved
    return res


def _x_fold_degenerate_aa_from_codon_table(x: int, table: str):
    """Extracts amino acids that are encoded by x
    different codons from an NCBI codon table.

    Args:
        x (int): fold degeneration to filter by.
        table (str, optional): NCBI table name as used in Bio.Data.CodonTable

    Returns:
        dict: maps each amino acid (keys) to a dict that specifies
              the degenerate codons for each site.

    """
    codon_table = unambiguous_dna_by_name[table].forward_table
    reverse_table = {}
    for codon, aa in codon_table.items():
        reverse_table.setdefault(aa, [])
        reverse_table[aa].append(codon)
    x_fold_table = {}
    for aa, codons in reverse_table.items():
        if len(codons) >= x:
            degeneracy = _site_degeneracy(codons, x)
            if degeneracy:
                x_fold_table[aa] = degeneracy
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


def _triplet_with_context(seq: str, i: int):
    """Return the ith triplet in a sequence with the
    previous and following nucleotide

    Args:
        seq (str): Sequence
        i (int): 0-based index of the target triplet.

    Returns:
        str: The target triplet.
    """
    start = i * 3
    return seq[start-1:start+4]


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


def substitutions_per_ffds(
    a: str, b: str, table_a="Standard", table_b="Standard"
) -> ((int, int), [str, str]):
    """Estimates the numbers of neutral substitutions per site by counting
    the number of substitutions at four-fold degenerate sites sites.

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

    ffdc_a = _x_fold_degenerate_aa_from_codon_table(4, table_a)
    ffdc_b = _x_fold_degenerate_aa_from_codon_table(4, table_b)

    common_ffds = {}
    for aa in set(ffdc_a).intersection(set(ffdc_b)):
        for site, codons in ffdc_a[aa].items():
            if site in ffdc_b[aa]:
                common_codons = []
                for codon_set in codons:
                    if codon_set in ffdc_b[aa][site]:
                        common_codons.append(codon_set)
                if common_codons:
                    common_ffds.setdefault(aa, dict())
                    common_ffds[aa][site] = common_codons

    n_sites = 0
    n_sub = 0
    offset_a, offset_b = 0, 0
    # iterate over the aligned sequences
    for i, (ca, cb) in enumerate(zip(
        str(alignment[0][0]), str(alignment[0][1]))
    ):
        if ca == "-":
            offset_a += 1
        if cb == "-":
            offset_b += 1
        if ca == cb and ca in common_ffds:
            triplet_a, triplet_b = _triplet(a, i - offset_a), _triplet(b, i - offset_b)
            _, locs = _hamming_distance(triplet_a, triplet_b)
            n_sites += len(common_ffds[ca])
            for loc in locs:
                if loc in common_ffds[ca]:
                    n_sub += 1

    return (n_sub, n_sites), truncated


def substitutions_per_ffds_by_cpg_context(
    a: str, b: str, table_a="Standard", table_b="Standard"
):
    """Estimates the numbers of neutral substitutions per site by counting
    the number of substitutions at four fold degenerate sites sites.
    Differentiates between four fold degenerate sites in different potential
    CpG contexts: preceded by C and not followed by G (nonCpG), preceded by C
    but not followed by G (postC), followed by G but not preceded by C (preG),
    and preceded by C and followed by G (postCpreG).

    Note: the number of sites considered here may be substantially lower than
          in substitutions_per_ffds, as this function requires the sites
          preceeding and following site of a four fold degenerate site
          to be identical in the two aligned sequences.

    See also:
        Kondrashov FA, Ogurtsov AY, Kondrashov AS. Selection in favor
        of nucleotidesG and C diversifies evolution rates and levels of
        polymorphism at mammalian synonymous sites.
        J Theor Biol. 2006;240(4):616‐626. doi:10.1016/j.jtbi.2005.10.020

    Args:
        a (str): coding DNA sequence a
        b (str): coding DNA sequence b
        table_a (str, optional): NCBI table name as used in Bio.Data.CodonTable
        table_b (str, optional): NCBI table name as used in Bio.Data.CodonTable
    """

    def categorize_ffds(qintuplet_a, qintuplet_b, loc):
        if qintuplet_a[loc - 1] != qintuplet_b[loc - 1]:
            return None
        elif qintuplet_a[loc + 1] != qintuplet_b[loc + 1]:
            return None
        else:
            pre, post = qintuplet_a[loc - 1], qintuplet_a[loc + 1]

        if pre == 'C' and post == 'G':
            return "postCpreG"
        elif pre == 'C':
            return "postC"
        elif post == 'G':
            return "preG"
        else:
            return "nonCpG"

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
    alignment_len = len(str(alignment[0][1]))

    ffdc_a = _x_fold_degenerate_aa_from_codon_table(4, table_a)
    ffdc_b = _x_fold_degenerate_aa_from_codon_table(4, table_b)

    common_ffds = {}
    for aa in set(ffdc_a).intersection(set(ffdc_b)):
        for site, codons in ffdc_a[aa].items():
            if site in ffdc_b[aa]:
                common_codons = []
                for codon_set in codons:
                    if codon_set in ffdc_b[aa][site]:
                        common_codons.append(codon_set)
                if common_codons:
                    common_ffds.setdefault(aa, dict())
                    common_ffds[aa][site] = common_codons

    subs_and_sites = {
        "nonCpG": [0, 0],
        "postC": [0, 0],
        "preG": [0, 0],
        "postCpreG": [0, 0]
    }
    offset_a, offset_b = 0, 0

    # iterate over the aligned sequences
    for i, (ca, cb) in enumerate(zip(
        str(alignment[0][0]), str(alignment[0][1]))
    ):
        if ca == "-":
            offset_a += 1
        if cb == "-":
            offset_b += 1
        if i == 0 or i == alignment_len - 1:
            continue
        if ca == cb and ca in common_ffds:
            qintuplet_a = _triplet_with_context(a, i - offset_a)
            qintuplet_b = _triplet_with_context(b, i - offset_b)
            _, locs = _hamming_distance(qintuplet_a[1:4], qintuplet_b[1:4])
            for loc in common_ffds[ca]:
                category = categorize_ffds(qintuplet_a, qintuplet_b, loc + 1)
                if category is None:
                    continue
                else:
                    subs_and_sites[category][1] += 1
                    subs_and_sites[category][0] += (loc in locs)

    return subs_and_sites
