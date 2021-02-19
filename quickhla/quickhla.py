#!/usr/bin/python
# -*- coding: UTF-8 -*-


"""
Author: Krisian Ullrich
date: February 2021
email: ullrich@evolbio.mpg.de
License: MIT

The MIT License (MIT)

Copyright (c) 2021 Kristian Ullrich

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


import os
import sys
import argparse
import numpy as np
from Bio import SeqIO


def build_taxid_tree(hla_0d, hla_2d, hla_4d, hla_6d, hla_8d):
    t = 1
    tree = {}
    for x_idx, y in enumerate(hla_0d):
        id_0d = y
        id_2d = id_0d + '*' + hla_2d[x_idx]
        id_4d = id_2d + ':' + hla_4d[x_idx]
        id_6d = id_4d + ':' + hla_6d[x_idx]
        id_8d = id_6d + ':' + hla_8d[x_idx]
        if id_0d not in tree:
            parent_t = 1
            t += 1
            t_type = 'kingdom'
            tree[id_0d] = [t, parent_t, t_type, id_0d]
        if id_2d not in tree:
            parent_t = tree[id_0d][0]
            t += 1
            t_type = 'phylum'
            tree[id_2d] = [t, parent_t, t_type, id_2d]
        if id_4d not in tree:
            parent_t = tree[id_2d][0]
            t += 1
            t_type = 'class'
            tree[id_4d] = [t, parent_t, t_type, id_4d]
        if id_6d not in tree:
            parent_t = tree[id_4d][0]
            t += 1
            t_type = 'order'
            tree[id_6d] = [t, parent_t, t_type, id_6d]
        if id_8d not in tree:
            parent_t = tree[id_6d][0]
            t += 1
            t_type = 'family'
            tree[id_8d] = [t, parent_t, t_type, id_8d]
    return tree


def parse_sam(samfile):
    hit_dict = {}
    with open(samfile, 'r') as inhandle:
        for line in inhandle:
            if line[0] != '@':
                t = line.strip().split('\t')[2].split('|')[0]
                if t == '*':
                    continue
                id_0d = t.split('*')[0]
                id_2d = id_0d + '*' + t.split('*')[1].split(':')[0]
                id_4d = id_2d + ':' + t.split('*')[1].split(':')[1]
                id_6d = id_4d + ':' + t.split('*')[1].split(':')[2]
                id_8d = id_6d + ':' + t.split('*')[1].split(':')[3]
                if id_0d in hit_dict:
                    hit_dict[id_0d][0] += 1
                    if id_2d in hit_dict[id_0d][1]:
                        hit_dict[id_0d][1][id_2d][0] += 1
                        if id_4d in hit_dict[id_0d][1][id_2d][1]:
                            hit_dict[id_0d][1][id_2d][1][id_4d][0] += 1
                            if id_6d in hit_dict[id_0d][1][id_2d][1][id_4d][1]:
                                hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][0] += 1
                                if id_8d in hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1]:
                                    hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1][id_8d][0] += 1
                                if id_8d not in hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1]:
                                    hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1][id_8d] = [1, {}]
                            if id_6d not in hit_dict[id_0d][1][id_2d][1][id_4d][1]:
                                hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d] = [1, {}]
                                hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1][id_8d] = [1, {}]
                        if id_4d not in hit_dict[id_0d][1][id_2d][1]:
                            hit_dict[id_0d][1][id_2d][1][id_4d] = [1, {}]
                            hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d] = [1, {}]
                            hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1][id_8d] = [1, {}]
                    if id_2d not in hit_dict[id_0d][1]:
                        hit_dict[id_0d][1][id_2d] = [1, {}]
                        hit_dict[id_0d][1][id_2d][1][id_4d] = [1, {}]
                        hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d] = [1, {}]
                        hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1][id_8d] = [1, {}]
                if id_0d not in hit_dict:
                    hit_dict[id_0d] = [1, {}]
                    hit_dict[id_0d][1][id_2d] = [1, {}]
                    hit_dict[id_0d][1][id_2d][1][id_4d] = [1, {}]
                    hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d] = [1, {}]
                    hit_dict[id_0d][1][id_2d][1][id_4d][1][id_6d][1][id_8d] = [1, {}]
    return hit_dict


def weight_hit_dict(hit_dict, weights):
    for id_0d_k in hit_dict.keys():
        hit_dict[id_0d_k][0] = hit_dict[id_0d_k][0] * weights[id_0d_k]
        for id_2d_k in hit_dict[id_0d_k][1]:
            hit_dict[id_0d_k][1][id_2d_k][0] = hit_dict[id_0d_k][1][id_2d_k][0] * weights[id_2d_k]
            for id_4d_k in hit_dict[id_0d_k][1][id_2d_k][1]:
                hit_dict[id_0d_k][1][id_2d_k][1][id_4d_k][0] = hit_dict[id_0d_k][1][id_2d_k][1][id_4d_k][0] * weights[
                    id_4d_k]
                for id_6d_k in hit_dict[id_0d_k][1][id_2d_k][1][id_4d_k][1]:
                    hit_dict[id_0d_k][1][id_2d_k][1][id_4d_k][1][id_6d_k][0] = hit_dict[
                        id_0d_k][1][id_2d_k][1][id_4d_k][1][id_6d_k][0] * weights[id_6d_k]
    return hit_dict


def write_hit_dict(hit_dict, outfile, n):
    with open(outfile, 'w') as outhandle:
        for id_0d_k in sorted(list(hit_dict.keys())):
            # define empty hit arrays for each digit
            id_2d_hits = []
            id_4d_hits = []
            id_6d_hits = []
            id_8d_hits = []
            outhandle.write('%s\t%i\n' % (id_0d_k, hit_dict[id_0d_k][0]))
            # print('%s\t%i\n' % (id_0d_k, hit_dict[id_0d_k][0]))
            id_2d_values = [[x[0], x[1][0]] for x in
                            sorted(hit_dict[id_0d_k][1].items(), key=lambda x: x[1][0], reverse=True)]
            id_2d_hits.extend(sorted(id_2d_values, key=lambda x: x[1], reverse=True))
            for id_2d_k, id_2d_n in id_2d_values:
                id_4d_values = [[x[0], x[1][0]] for x in
                                sorted(hit_dict[id_0d_k][1][id_2d_k][1].items(), key=lambda x: x[1][0], reverse=True)]
                id_2d_k_id_4d_hits = sorted(id_4d_values, key=lambda x: x[1], reverse=True)
                id_4d_hits.extend(id_2d_k_id_4d_hits)
                for id_4d_k, id_4d_n in id_4d_values:
                    id_6d_values = [[x[0], x[1][0]] for x in
                                    sorted(hit_dict[id_0d_k][1][id_2d_k][1][id_4d_k][1].items(), key=lambda x: x[1][0],
                                           reverse=True)]
                    id_2d_k_id_4d_k_id_6d_hits = sorted(id_6d_values, key=lambda x: x[1], reverse=True)
                    id_6d_hits.extend(id_2d_k_id_4d_k_id_6d_hits)
                    for id_6d_k, id_6d_n in id_6d_values:
                        id_8d_values = [[x[0], x[1][0]] for x in
                                        sorted(hit_dict[id_0d_k][1][id_2d_k][1][id_4d_k][1][id_6d_k][1].items(),
                                               key=lambda x: x[1][0], reverse=True)]
                        id_2d_k_id_4d_k_id_6d_k_id_8d_hits = sorted(id_8d_values, key=lambda x: x[1], reverse=True)
                        id_8d_hits.extend(id_2d_k_id_4d_k_id_6d_k_id_8d_hits)
            id_2d_hits = sorted(id_2d_hits, key=lambda x: x[1], reverse=True)[:n]
            for id_2d_h, id_2d_n in id_2d_hits:
                outhandle.write('\t%s\t%i\n' % (id_2d_h, id_2d_n))
                # print('\t%s\t%i\n' % (id_2d_h, id_2d_n))
            id_4d_hits = sorted(id_4d_hits, key=lambda x: x[1], reverse=True)[:n]
            for id_4d_h, id_4d_n in id_4d_hits:
                outhandle.write('\t%s\t%i\n' % (id_4d_h, id_4d_n))
                # print('\t%s\t%i\n' % (id_4d_h, id_4d_n))
            id_6d_hits = sorted(id_6d_hits, key=lambda x: x[1], reverse=True)[:n]
            for id_6d_h, id_6d_n in id_6d_hits:
                outhandle.write('\t%s\t%i\n' % (id_6d_h, id_6d_n))
                # print('\t%s\t%i\n' % (id_6d_h, id_6d_n))
            id_8d_hits = sorted(id_8d_hits, key=lambda x: x[1], reverse=True)[:n]
            for id_8d_h, id_8d_n in id_8d_hits:
                outhandle.write('\t%s\t%i\n' % (id_8d_h, id_8d_n))
                # print('\t%s\t%i\n' % (id_8d_h, id_8d_n))


def tree_to_names(tree, foo):
    with open(foo, 'w') as outhandle:
        outhandle.write('1\t|\thlaroot\t|\thlaroot\t|\tscientific name\t|\n')
        for i in tree.keys():
            outhandle.write(str(tree[i][0]) + '\t|\t' +
                            str(tree[i][3]) + '\t|\t' +
                            str(tree[i][3]) + '\t|\tscientific name\t|\n')


def tree_to_nodes(tree, foo):
    with open(foo, 'w') as outhandle:
        outhandle.write('1\t|\t0\t|\tno rank\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n')
        for i in tree.keys():
            outhandle.write(str(tree[i][0]) + '\t|\t' +
                            str(tree[i][1]) + '\t|\t' +
                            str(tree[i][2]) + '\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\t\t|\n')


def create_ids(hla_0d, hla_2d, hla_4d, hla_6d, hla_8d, tree, seq_array, level):
    for x_idx, y in enumerate(hla_0d):
        id_0d = y
        id_2d = y + '*' + hla_2d[x_idx]
        id_4d = y + '*' + hla_2d[x_idx] + ':' + hla_4d[x_idx]
        id_6d = y + '*' + hla_2d[x_idx] + ':' + hla_4d[x_idx] + ':' + hla_6d[x_idx]
        id_8d = y + '*' + hla_2d[x_idx] + ':' + hla_4d[x_idx] + ':' + hla_6d[x_idx] + ':' + hla_8d[x_idx]
        if level == '0d':
            seq_array[x_idx].id = id_8d + '|kraken:taxid|' + str(tree[id_0d][0]) + '|' + id_0d
            seq_array[x_idx].name = id_8d + '|kraken:taxid|' + str(tree[id_0d][0]) + '|' + id_0d
            seq_array[x_idx].description = id_8d + '|kraken:taxid|' + str(tree[id_0d][0]) + '|' + id_0d
        if level == '2d':
            seq_array[x_idx].id = id_8d + '|kraken:taxid|' + str(tree[id_2d][0]) + '|' + id_2d
            seq_array[x_idx].name = id_8d + '|kraken:taxid|' + str(tree[id_2d][0]) + '|' + id_2d
            seq_array[x_idx].description = id_8d + '|kraken:taxid|' + str(tree[id_2d][0]) + '|' + id_2d
        if level == '4d':
            seq_array[x_idx].id = id_8d + '|kraken:taxid|' + str(tree[id_4d][0]) + '|' + id_4d
            seq_array[x_idx].name = id_8d + '|kraken:taxid|' + str(tree[id_4d][0]) + '|' + id_4d
            seq_array[x_idx].description = id_8d + '|kraken:taxid|' + str(tree[id_4d][0]) + '|' + id_4d
        if level == '6d':
            seq_array[x_idx].id = id_8d + '|kraken:taxid|' + str(tree[id_6d][0]) + '|' + id_6d
            seq_array[x_idx].name = id_8d + '|kraken:taxid|' + str(tree[id_6d][0]) + '|' + id_6d
            seq_array[x_idx].description = id_8d + '|kraken:taxid|' + str(tree[id_6d][0]) + '|' + id_6d
        if level == '8d':
            seq_array[x_idx].id = id_8d + '|kraken:taxid|' + str(tree[id_8d][0]) + '|' + id_8d
            seq_array[x_idx].name = id_8d + '|kraken:taxid|' + str(tree[id_8d][0]) + '|' + id_8d
            seq_array[x_idx].description = id_8d + '|kraken:taxid|' + str(tree[id_8d][0]) + '|' + id_8d
    return seq_array


def get_weights(seq, hla_0d, hla_2d, hla_4d, hla_6d, hla_8d):
    # combine
    hla_0d_ids = hla_0d
    hla_2d_ids = [x + '*' + y for x, y in zip(hla_0d_ids, hla_2d)]
    hla_4d_ids = [x + ':' + y for x, y in zip(hla_2d_ids, hla_4d)]
    hla_6d_ids = [x + ':' + y for x, y in zip(hla_4d_ids, hla_6d)]
    hla_8d_ids = [x + ':' + y for x, y in zip(hla_6d_ids, hla_8d)]
    # get total seq length
    seq_len = np.sum([len(x) for x in seq])
    # get hla_0d seq length
    hla_0d_seq_len = {}
    hla_2d_seq_len = {}
    hla_4d_seq_len = {}
    hla_6d_seq_len = {}
    hla_8d_seq_len = {}
    for x, y in zip(seq, hla_8d_ids):
        y_0d = y.split('*')[0]
        y_2d = y.split(':')[0]
        y_4d = y_2d + ':' + y.split(':')[1]
        y_6d = y_4d + ':' + y.split(':')[2]
        y_8d = y
        if y_0d in hla_0d_seq_len:
            hla_0d_seq_len[y_0d] += len(x)
        if y_0d not in hla_0d_seq_len:
            hla_0d_seq_len[y_0d] = len(x)
        if y_2d in hla_2d_seq_len:
            hla_2d_seq_len[y_2d] += len(x)
        if y_2d not in hla_2d_seq_len:
            hla_2d_seq_len[y_2d] = len(x)
        if y_4d in hla_4d_seq_len:
            hla_4d_seq_len[y_4d] += len(x)
        if y_4d not in hla_4d_seq_len:
            hla_4d_seq_len[y_4d] = len(x)
        if y_6d in hla_6d_seq_len:
            hla_6d_seq_len[y_6d] += len(x)
        if y_6d not in hla_6d_seq_len:
            hla_6d_seq_len[y_6d] = len(x)
        if y_8d in hla_8d_seq_len:
            hla_8d_seq_len[y_8d] += len(x)
        if y_8d not in hla_8d_seq_len:
            hla_8d_seq_len[y_8d] = len(x)
    weights = {}
    # 0d
    for k, v in hla_0d_seq_len.items():
        weights[k] = seq_len/v
    # 2d
    for k, v in hla_2d_seq_len.items():
        y_0d = k.split('*')[0]
        weights[k] = hla_0d_seq_len[y_0d]/v
    # 4d
    for k, v in hla_4d_seq_len.items():
        y_2d = k.split(':')[0]
        weights[k] = hla_2d_seq_len[y_2d]/v
    # 6d
    for k, v in hla_6d_seq_len.items():
        y_2d = k.split(':')[0]
        y_4d = y_2d + ':' + k.split(':')[1]
        weights[k] = hla_4d_seq_len[y_4d]/v
    # 8d
    for k, v in hla_8d_seq_len.items():
        y_2d = k.split(':')[0]
        y_4d = y_2d + ':' + k.split(':')[1]
        y_6d = y_4d + ':' + k.split(':')[2]
        weights[k] = hla_6d_seq_len[y_6d]/v
    return weights


def write_weights(weights, foo):
    with open(foo, 'w') as outhandle:
        for k in weights:
            outhandle.write('%s\t%s\n' % (k, str(weights[k])))


def read_weights(foo):
    weights = {}
    with open(foo, 'r') as inhandle:
        for line in inhandle:
            k = line.strip().split('\t')[0]
            v = line.strip().split('\t')[1]
            weights[k] = float(v)
    return weights


def build_kraken(args, parser):
    """
    creates kraken database files

    :param args:
    :param parser:
    :return:
    """
    if not args.gen:
        parser.print_help()
        sys.exit('\nPlease specify input gen file ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_gen.fasta')
    if not args.nuc:
        parser.print_help()
        sys.exit('\nPlease specify input gen file ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/hla_nuc.fasta')
    if not args.kb:
        args.kb = 'kraken2-build'
    if not args.bb:
        args.bb = 'bowtie2-build'
    if not args.hb:
        args.hb = 'hisat2-build'
    print(args)
    for kl in args.kl:
        os.system('mkdir -p ' + args.o + '/hla.gen.0d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.gen.2d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.gen.4d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.gen.6d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.gen.8d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.nuc.0d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.nuc.2d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.nuc.4d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.nuc.6d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.nuc.8d.' + str(kl))
        os.system('mkdir -p ' + args.o + '/hla.gen.0d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.gen.2d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.gen.4d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.gen.6d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.gen.8d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.nuc.0d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.nuc.2d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.nuc.4d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.nuc.6d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.nuc.8d.' + str(kl) + '/data')
        os.system('mkdir -p ' + args.o + '/hla.gen.0d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.gen.2d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.gen.4d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.gen.6d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.gen.8d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.nuc.0d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.nuc.2d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.nuc.4d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.nuc.6d.' + str(kl) + '/taxonomy')
        os.system('mkdir -p ' + args.o + '/hla.nuc.8d.' + str(kl) + '/taxonomy')
    gen = list(SeqIO.parse(args.gen, 'fasta'))
    nuc = list(SeqIO.parse(args.nuc, 'fasta'))
    # separate record names into 0d, 2d, 4d, 6d, 8d
    gen_hla = [x.description.split(' ')[1] for x in gen]
    nuc_hla = [x.description.split(' ')[1] for x in nuc]
    # 0d
    gen_hla_0d = [x.split('*')[0] for x in gen_hla]
    nuc_hla_0d = [x.split('*')[0] for x in nuc_hla]
    # fill NA to 8d for all entries
    gen_hla_sub = [x.split('*')[1] for x in gen_hla]
    nuc_hla_sub = [x.split('*')[1] for x in nuc_hla]
    gen_hla_sub_len = [len(x.split(':')) for x in gen_hla_sub]
    nuc_hla_sub_len = [len(x.split(':')) for x in nuc_hla_sub]
    for x_idx, y_len in enumerate(gen_hla_sub_len):
        if y_len != 4:
            gen_hla_sub[x_idx] += ''.join(np.repeat(':NA', 4 - y_len))
    for x_idx, y_len in enumerate(nuc_hla_sub_len):
        if y_len != 4:
            nuc_hla_sub[x_idx] += ''.join(np.repeat(':NA', 4 - y_len))
    # 2d
    gen_hla_2d = [x.split(':')[0] for x in gen_hla_sub]
    nuc_hla_2d = [x.split(':')[0] for x in nuc_hla_sub]
    # 4d
    gen_hla_4d = [x.split(':')[1] for x in gen_hla_sub]
    nuc_hla_4d = [x.split(':')[1] for x in nuc_hla_sub]
    # 6d
    gen_hla_6d = [x.split(':')[2] for x in gen_hla_sub]
    nuc_hla_6d = [x.split(':')[2] for x in nuc_hla_sub]
    # 8d
    gen_hla_8d = [x.split(':')[3] for x in gen_hla_sub]
    nuc_hla_8d = [x.split(':')[3] for x in nuc_hla_sub]
    # build taxid tree
    gen_tree = build_taxid_tree(gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d)
    nuc_tree = build_taxid_tree(nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d)
    # get weights
    gen_weights = get_weights(gen, gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d)
    nuc_weights = get_weights(nuc, nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d)
    # create kraken database
    for kl, ml, ms in zip(args.kl, args.ml, args.ms):
        # 0d
        # write names.dmp and nodes.dmp
        tree_to_names(gen_tree, args.o + '/hla.gen.0d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_names(nuc_tree, args.o + '/hla.nuc.0d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_nodes(gen_tree, args.o + '/hla.gen.0d.' + str(kl) + '/taxonomy/nodes.dmp')
        tree_to_nodes(nuc_tree, args.o + '/hla.nuc.0d.' + str(kl) + '/taxonomy/nodes.dmp')
        # change fasta file id and description
        gen_0d = list(SeqIO.parse(args.gen, 'fasta'))
        nuc_0d = list(SeqIO.parse(args.nuc, 'fasta'))
        gen_0d = create_ids(gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d, gen_tree, gen_0d, '0d')
        nuc_0d = create_ids(nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d, nuc_tree, nuc_0d, '0d')
        # write fasta file
        SeqIO.write(gen_0d, args.o + '/hla.gen.0d.' + str(kl) + '/data/data.fasta', 'fasta')
        SeqIO.write(nuc_0d, args.o + '/hla.nuc.0d.' + str(kl) + '/data/data.fasta', 'fasta')
        write_weights(gen_weights, args.o + '/hla.gen.0d.' + str(kl) + '/data/data.weights')
        write_weights(nuc_weights, args.o + '/hla.nuc.0d.' + str(kl) + '/data/data.weights')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.0d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.0d.' + str(kl) + '/data/data')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.0d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.0d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.0d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.0d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.0d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.0d.' + str(kl) + '/data/data')
        # kraken-build add library
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.gen.0d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.gen.0d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.nuc.0d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.nuc.0d.' + str(kl))
        # kraken-build build
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.gen.0d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.nuc.0d.' + str(kl))
        # 2d
        # write names.dmp and nodes.dmp
        tree_to_names(gen_tree, args.o + '/hla.gen.2d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_names(nuc_tree, args.o + '/hla.nuc.2d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_nodes(gen_tree, args.o + '/hla.gen.2d.' + str(kl) + '/taxonomy/nodes.dmp')
        tree_to_nodes(nuc_tree, args.o + '/hla.nuc.2d.' + str(kl) + '/taxonomy/nodes.dmp')
        # change fasta file id and description
        gen_2d = list(SeqIO.parse(args.gen, 'fasta'))
        nuc_2d = list(SeqIO.parse(args.nuc, 'fasta'))
        gen_2d = create_ids(gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d, gen_tree, gen_2d, '2d')
        nuc_2d = create_ids(nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d, nuc_tree, nuc_2d, '2d')
        # write fasta file
        SeqIO.write(gen_2d, args.o + '/hla.gen.2d.' + str(kl) + '/data/data.fasta', 'fasta')
        SeqIO.write(nuc_2d, args.o + '/hla.nuc.2d.' + str(kl) + '/data/data.fasta', 'fasta')
        write_weights(gen_weights, args.o + '/hla.gen.2d.' + str(kl) + '/data/data.weights')
        write_weights(nuc_weights, args.o + '/hla.nuc.2d.' + str(kl) + '/data/data.weights')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.2d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.2d.' + str(kl) + '/data/data')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.2d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.2d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.2d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.2d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.2d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.2d.' + str(kl) + '/data/data')
        # kraken-build add library
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.gen.2d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.gen.2d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.nuc.2d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.nuc.2d.' + str(kl))
        # kraken-build build
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.gen.2d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.nuc.2d.' + str(kl))
        # 4d
        # write names.dmp and nodes.dmp
        tree_to_names(gen_tree, args.o + '/hla.gen.4d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_names(nuc_tree, args.o + '/hla.nuc.4d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_nodes(gen_tree, args.o + '/hla.gen.4d.' + str(kl) + '/taxonomy/nodes.dmp')
        tree_to_nodes(nuc_tree, args.o + '/hla.nuc.4d.' + str(kl) + '/taxonomy/nodes.dmp')
        # change fasta file id and description
        gen_4d = list(SeqIO.parse(args.gen, 'fasta'))
        nuc_4d = list(SeqIO.parse(args.nuc, 'fasta'))
        gen_4d = create_ids(gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d, gen_tree, gen_4d, '4d')
        nuc_4d = create_ids(nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d, nuc_tree, nuc_4d, '4d')
        # write fasta file
        SeqIO.write(gen_4d, args.o + '/hla.gen.4d.' + str(kl) + '/data/data.fasta', 'fasta')
        SeqIO.write(nuc_4d, args.o + '/hla.nuc.4d.' + str(kl) + '/data/data.fasta', 'fasta')
        write_weights(gen_weights, args.o + '/hla.gen.4d.' + str(kl) + '/data/data.weights')
        write_weights(nuc_weights, args.o + '/hla.nuc.4d.' + str(kl) + '/data/data.weights')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.4d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.4d.' + str(kl) + '/data/data')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.4d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.4d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.4d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.4d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.4d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.4d.' + str(kl) + '/data/data')
        # kraken-build add library
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.gen.4d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.gen.4d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.nuc.4d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.nuc.4d.' + str(kl))
        # kraken-build build
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.gen.4d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.nuc.4d.' + str(kl))
        # 6d
        # write names.dmp and nodes.dmp
        tree_to_names(gen_tree, args.o + '/hla.gen.6d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_names(nuc_tree, args.o + '/hla.nuc.6d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_nodes(gen_tree, args.o + '/hla.gen.6d.' + str(kl) + '/taxonomy/nodes.dmp')
        tree_to_nodes(nuc_tree, args.o + '/hla.nuc.6d.' + str(kl) + '/taxonomy/nodes.dmp')
        # change fasta file id and description
        gen_6d = list(SeqIO.parse(args.gen, 'fasta'))
        nuc_6d = list(SeqIO.parse(args.nuc, 'fasta'))
        gen_6d = create_ids(gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d, gen_tree, gen_6d, '6d')
        nuc_6d = create_ids(nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d, nuc_tree, nuc_6d, '6d')
        # write fasta file
        SeqIO.write(gen_6d, args.o + '/hla.gen.6d.' + str(kl) + '/data/data.fasta', 'fasta')
        SeqIO.write(nuc_6d, args.o + '/hla.nuc.6d.' + str(kl) + '/data/data.fasta', 'fasta')
        write_weights(gen_weights, args.o + '/hla.gen.6d.' + str(kl) + '/data/data.weights')
        write_weights(nuc_weights, args.o + '/hla.nuc.6d.' + str(kl) + '/data/data.weights')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.6d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.6d.' + str(kl) + '/data/data')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.6d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.6d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.6d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.6d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.6d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.6d.' + str(kl) + '/data/data')
        # kraken-build add library
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.gen.6d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.gen.6d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.nuc.6d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.nuc.6d.' + str(kl))
        # kraken-build build
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.gen.6d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.nuc.6d.' + str(kl))
        # 8d
        tree_to_names(gen_tree, args.o + '/hla.gen.8d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_names(nuc_tree, args.o + '/hla.nuc.8d.' + str(kl) + '/taxonomy/names.dmp')
        tree_to_nodes(gen_tree, args.o + '/hla.gen.8d.' + str(kl) + '/taxonomy/nodes.dmp')
        tree_to_nodes(nuc_tree, args.o + '/hla.nuc.8d.' + str(kl) + '/taxonomy/nodes.dmp')
        # change fasta file id and description
        gen_8d = list(SeqIO.parse(args.gen, 'fasta'))
        nuc_8d = list(SeqIO.parse(args.nuc, 'fasta'))
        gen_8d = create_ids(gen_hla_0d, gen_hla_2d, gen_hla_4d, gen_hla_6d, gen_hla_8d, gen_tree, gen_8d, '8d')
        nuc_8d = create_ids(nuc_hla_0d, nuc_hla_2d, nuc_hla_4d, nuc_hla_6d, nuc_hla_8d, nuc_tree, nuc_8d, '8d')
        # write fasta file
        SeqIO.write(gen_8d, args.o + '/hla.gen.8d.' + str(kl) + '/data/data.fasta', 'fasta')
        SeqIO.write(nuc_8d, args.o + '/hla.nuc.8d.' + str(kl) + '/data/data.fasta', 'fasta')
        write_weights(gen_weights, args.o + '/hla.gen.8d.' + str(kl) + '/data/data.weights')
        write_weights(nuc_weights, args.o + '/hla.nuc.8d.' + str(kl) + '/data/data.weights')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.8d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.8d.' + str(kl) + '/data/data')
        os.system(args.bb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.8d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.8d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.gen.8d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.gen.8d.' + str(kl) + '/data/data')
        os.system(args.hb + ' --threads ' + str(args.t) + ' ' +
                  args.o + '/hla.nuc.8d.' + str(kl) + '/data/data.fasta ' +
                  args.o + '/hla.nuc.8d.' + str(kl) + '/data/data')
        # kraken-build add library
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.gen.8d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.gen.8d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --add-to-library ' + args.o + '/hla.nuc.8d.' +
                  str(kl) + '/data/data.fasta' + ' --db ' + args.o + '/hla.nuc.8d.' + str(kl))
        # kraken-build build
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.gen.8d.' + str(kl))
        os.system(args.kb + ' --kmer-len ' + str(kl) + ' --threads ' + str(args.t) +
                  ' --minimizer-len ' + str(ml) + ' --minimizer-spaces ' +
                  str(ms) + ' --build --db ' + args.o + '/hla.nuc.8d.' + str(kl))


def classify(args, parser):
    """
    classify fastq files

    :param args:
    :param parser:
    :return:
    """
    if not args.f:
        parser.print_help()
        sys.exit('\nPlease specify fastq file')
    if not args.db:
        parser.print_help()
        sys.exit('\nPlease specify database to use for classification')
    if not args.kb:
        args.kb = 'kraken2'
    if not args.bb:
        args.bb = 'bowtie2'
    if not args.hb:
        args.hb = 'hisat2'
    print(args)
    # kraken2 search
    os.system(
        args.kb +
        ' --use-names' +
        ' --threads ' + str(args.t) +
        ' --db ' + args.d + '/' + args.db +
        ' --classified-out ' + args.o + '.' + args.db + '#.fq' +
        ' --report ' + args.o + '.' + args.db + '.report ' +
        ' --output - ' +
        ' --paired ' + args.f + ' ' + args.r)
    print('kraken2 search finished')
    if args.alg == 'none':
        print('quickhla finished')
    else:
        # hisat2 or bowtie2 search
        if args.alg == 'hisat2':
            os.system(
                args.hb + ' ' +
                args.ho + ' ' +
                ' --threads ' + str(args.t) +
                ' -x ' + args.d + '/' + args.db + '/data/data' +
                ' -1 ' + args.o + '.' + args.db + '_1.fq' +
                ' -2 ' + args.o + '.' + args.db + '_2.fq' +
                ' -S ' + args.o + '.' + args.db + '.sam')
        if args.alg == 'bowtie2':
            os.system(
                args.bb + ' ' +
                args.bo + ' ' +
                ' --threads ' + str(args.t) +
                ' -x ' + args.d + '/' + args.db + '/data/data' +
                ' -1 ' + args.o + '.' + args.db + '_1.fq' +
                ' -2 ' + args.o + '.' + args.db + '_2.fq' +
                ' -S ' + args.o + '.' + args.db + '.sam')
        # create hit dictionary
        hit_dict = parse_sam(args.o + '.' + args.db + '.sam')
        # get weights
        if args.w:
            weights = read_weights(args.d + '/' + args.db + '/data/data.weights')
            hit_dict = weight_hit_dict(hit_dict, weights)
        write_hit_dict(hit_dict, args.o + '.' + args.db + '.class.txt', args.n)


def subparser(subparsers):
    # build; parser
    parser_build = subparsers.add_parser('build', help='build help')
    parser_build.add_argument('-gen', help='input gen file [mandatory]')
    parser_build.add_argument('-nuc', help='input nuc file [mandatory]')
    parser_build.add_argument('-o', help='output kraken2 database directory [default: hla.db]', default='hla.db')
    parser_build.add_argument('-kb', help='specify kraken2-build binary [if not given assumes to be in PATH]')
    parser_build.add_argument('-bb', help='specify bowtie2-build binary [if not given assumes to be in PATH]')
    parser_build.add_argument('-hb', help='specify hisat2-build binary [if not given assumes to be in PATH]')
    parser_build.add_argument('-kl',
                              help='kraken2 kmer-len array [default: 35]'
                                   ' can be an array of values and needs to be same length as -ml and -ms',
                              default=[35])
    parser_build.add_argument('-ml', help='kraken2 minimizer-len array [default: 31]'
                                          ' if array needs to be same length as -kl and -ms',
                              default=[31])
    parser_build.add_argument('-ms', help='kraken2 minimizer-space array [default: 7]'
                                          ' if array needs to be same length as -kl and -ml',
                              default=[7])
    parser_build.add_argument('-t', help='specify number threads [default: 1]', default=1, type=int)
    parser_build.set_defaults(func=build_kraken)
    # classify; parser
    parser_classify = subparsers.add_parser('classify', help='classify help')
    parser_classify.add_argument('-f', help='specify forward fastq [mandatory]')
    parser_classify.add_argument('-r', help='specify reverse fastq [mandatory]')
    parser_classify.add_argument('-d', help='specify db directory [mandatory]')
    parser_classify.add_argument('-db', help='specify db for classification [mandatory]')
    parser_classify.add_argument('-o', help='specify output prefix [default: out]', default='out')
    parser_classify.add_argument('-t', help='specify number threads [default: 1]', default=1, type=int)
    parser_classify.add_argument('-kb', help='specify kraken2 binary [if not given assumes to be in PATH]')
    parser_classify.add_argument('-alg', help='specify aligner [default: hisat2]', default='hisat2')
    parser_classify.add_argument('-w', help='apply weights on read counts [default: False]', action='store_true')
    parser_classify.add_argument('-bb', help='specify bowtie2 binary [if not given assumes to be in PATH]')
    parser_classify.add_argument('-bo', help='specify bowtie2 options '
                                             '[default: --very-fast --no-unal -k 1000]',
                                 default='--very-fast --no-unal -k 1000')
    parser_classify.add_argument('-hb', help='specify hisat2 binary [if not given assumes to be in PATH]')
    parser_classify.add_argument('-ho', help='specify hisat2 options '
                                             '[default: --fast -k 1000]',
                                 default='--fast -k 1000')
    parser_classify.add_argument('-n', help='specify number of top hits to report [default: show all]', type=int)
    parser_classify.set_defaults(func=classify)


def main():
    # top-level parser
    parser = argparse.ArgumentParser(prog='quickhla', usage='%(prog)s <sub-script> [options] [<arguments>...]',
                                     description='classifies HLA genes from fastq files')
    subparsers = parser.add_subparsers(title='sub-scripts', description='valid sub-scripts', help='sub-scripts help',
                                       dest='cmd')
    # sub-level parser
    subparser(subparsers)
    # get args
    args = parser.parse_args()
    # call function
    try:
        args.func(args, parser)
    except AttributeError:
        parser.print_help()
        parser.exit()


if __name__ == '__main__':
    main()
