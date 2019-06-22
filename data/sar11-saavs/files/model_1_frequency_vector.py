#!/usr/bin/env python
# -*- coding: utf-8
import sys
import numpy as np
import pandas as pd

# r is assumed 0
kappa = float(sys.argv[1]) # transition vs transversion bias

from anvio.sequence import Codon
from anvio.constants import codons, codon_to_AA, amino_acids

dist = Codon().get_codon_to_codon_dist_dictionary()

genes = [int(gene.strip()) for gene in open('List-core-799genes.txt').readlines()]
df = pd.read_csv('codon_frequencies.txt', sep='\t')
vector = df.loc[df['gene_callers_id'].isin(genes), :].iloc[:, 1:].sum(axis=0)
vector /= vector.sum()

exchange_rate = {}

for X in codons:
    for Y in codons:
        d = dist[X][Y][0]
        amino_acid_X = codon_to_AA[X]
        amino_acid_Y = codon_to_AA[Y]
        AAST = ''.join(sorted([amino_acid_X, amino_acid_Y]))

        if d != 1:
            continue
        if codons.index(Y) <= codons.index(X):
            continue
        if amino_acid_X == amino_acid_Y:
            continue

        t = 'transversion' if dist[X][Y][1] else 'transition'
        P = kappa / (kappa + 2) if t == 'transition' else 1 / (kappa + 2)

        fX = vector[X]
        fY = vector[Y]
        score = (fX + fY)/2 * P

        if AAST not in exchange_rate:
            exchange_rate[AAST] = score
        else:
            exchange_rate[AAST] += score

top_twenty_five = ["IleVal",
                   "AspGlu",
                   "AsnAsp",
                   "AsnLys",
                   "AsnSer",
                   "ArgLys",
                   "GluLys",
                   "SerThr",
                   "IleLeu",
                   "AlaThr",
                   "IleThr",
                   "GlnLys",
                   "AlaSer",
                   "LeuPhe",
                   "AlaVal",
                   "IleMet",
                   "GlnGlu",
                   "PheTyr",
                   "LeuVal",
                   "GlySer",
                   "HisTyr",
                   "LeuSer",
                   "LysThr",
                   "AsnThr",
                   "ProSer"]

# normalize rates
exchange_rate = {aast: rate/sum(exchange_rate.values()) for aast, rate in exchange_rate.items()}

x = top_twenty_five
y = [exchange_rate[aast] for aast in top_twenty_five]

df = pd.DataFrame(list(zip(*[x,y])), columns=('AAST', 'fraction'))
df.to_csv('neutral_freq_dist_model_1.txt', sep='\t')
