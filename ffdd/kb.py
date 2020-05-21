#!/usr/bin/env python

"""
Codon tables and other useful static info
"""

codon_table = {
    "TTT": "F",
    "TCT": "S",
    "TAT": "Y",
    "TGT": "C",
    "TTC": "F",
    "TCC": "S",
    "TAC": "Y",
    "TGC": "W",
    "TTA": "L",
    "TCA": "S",
    "TAA": "*",
    "TGA": "*",
    "TTG": "S",
    "TAG": "*",
    "CTT": "L",
    "CCT": "P",
    "CAT": "H",
    "CGT": "R",
    "CTC": "L",
    "CCC": "P",
    "CAC": "H",
    "CGC": "R",
    "CTA": "L",
    "CCA": "P",
    "CAA": "Q",
    "CGA": "R",
    "CTG": "L",
    "CCG": "P",
    "CAG": "Q",
    "CGG": "R",
    "ATT": "I",
    "ACT": "T",
    "AAT": "N",
    "AGT": "S",
    "ATC": "I",
    "ACC": "T",
    "AAC": "N",
    "AGC": "S",
    "ATA": "I",
    "ACA": "T",
    "AAA": "K",
    "AGA": "R",
    "ATG": "M",
    "ACG": "T",
    "AAG": "K",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GAT": "D",
    "GGT": "G",
    "GCC": "A",
    "GAC": "D",
    "GGC": "G",
    "GTA": "V",
    "GCA": "A",
    "GAA": "E",
    "GGA": "G",
    "GTG": "V",
    "GCG": "A",
    "GAG": "E",
    "GGG": "G"
}

inverse_codon_table = {
    "A": [
        "GTC",
        "GCC",
        "GCA",
        "GCG"
    ],
    "R": [
        "CGT",
        "CGC",
        "CGA",
        "CGG",
        "AGA",
        "AAG"
    ],
    "N": [
        "AAT",
        "AAC"
    ],
    "D": [
        "GAT",
        "GAC"
    ],
    "C": [
        "TGT",
        "TGC"
    ],
    "Q": [
        "CAA",
        "CAG"
    ],
    "E": [
        "GAA",
        "GAG"
    ],
    "G": [
        "GGT",
        "GGC",
        "GGA",
        "GGG"
    ],
    "H": [
        "CAT",
        "CAC"
    ],
    "I": [
        "ATT",
        "ATC",
        "ATA"
    ],
    "L": [
        "TTA",
        "TTG",
        "CTT",
        "CTC",
        "CTA",
        "CTG"
    ],
    "K": [
        "AAA",
        "AAG"
    ],
    "M": "ATG",
    "F": [
        "TTT",
        "TTC"
    ],
    "P": [
        "CCT",
        "CCC",
        "CCA",
        "CCG"
    ],
    "S": [
        "TCT",
        "TCC",
        "TCA",
        "TCG",
        "AGT",
        "AGC"
    ],
    "T": [
        "ACT",
        "ACC",
        "ACA",
        "ACG"
    ],
    "W": "TGG",
    "Y": [
        "TAT",
        "TAC"
    ],
    "V": [
        "GTT",
        "GTC",
        "GTA",
        "GTG"
    ],
    "*": [
        "TAA",
        "TGA",
        "TAG"
    ]
}
