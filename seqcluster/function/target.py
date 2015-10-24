"""Get enrichment of targets doing permutations for significance"""
import os
import os.path as op
from collections import defaultdict
import gzip

from bcbio import utils


def targets_enrichment(args):
    if args.sps not in ["hsa"]:
        raise ValueError("Species not supported yet.")

    target, mirna_map, mirna_id = _get_files(args.annotation)
    mirna_id_dict = _get_mirbase_id(mirna_id, args.sps)
    mirna_map_dict = _get_target_id(mirna_map, args.sps)
    mirna_list = _get_mirna_input(args.input)
    res = defaultdict(list)
    collect = defaultdict(list)
    with open(target) as in_handle:
        in_handle.next()
        for line in in_handle:
            cols = line.split()
            mit_id = mirna_map_dict[cols[2]]
            for mit in mit_id:
                if mit in mirna_id_dict:
                    mirbase_id = mirna_id_dict[mit]
                    if mirbase_id in mirna_list:
                        if cols[-3] == "NULL":
                            cols[-3] = 0
                        collect[mirbase_id].append({'info': line, 'score': float(cols[-3])})
    for mir in collect:
        sorted_targets = sorted(collect[mir], key=lambda k: k['score'])
        for e in range(0, 100):
            line = sorted_targets[e]['info']
            cols = line.split()
            name = cols[0].split(".")[0]
            res[(name, 'info')] = line.strip()
            res[(name, 'mirs')].append(mir)

    with open(op.join(args.out, "pairs.tsv"), 'w') as out_handle:
       for gene, field in res:
            if field == "info":
                counts = len(res[(gene, 'mirs')])
                mirs = ",".join(res[(gene, 'mirs')])
                print >>out_handle, "%s\t%s\t%s\t%s" %  (mirs, counts,  gene, res[(gene, 'info')])


def _get_mirna_input(fn):
    need = set()
    with open(fn) as in_handle:
        for line in in_handle:
            if len(line.split("-")) > 2:
                need.add(line.strip())
    return need


def _get_target_id(fn,sps):
    mapping = defaultdict(list)
    with open(fn) as in_handle:
        for line in in_handle:
            cols = line.split()
            mapping[cols[1]].append(cols[-1].strip())
    return mapping


def _get_mirbase_id(fn, sps):
    mapping = dict()
    with gzip.open(fn, 'rb') as in_handle:
        for line in in_handle:
            cols = line.split()
            if cols[1].startswith(sps):
                cols = line.split()
                mapping[cols[3]] = cols[1]
    return mapping


def _get_files(path):
    target = op.join(path, "Summary_Counts_sps.txt")
    mirna_map = op.join(path, "miR_Family_Info_sps.txt")
    mirna_id = op.join(path, "mirna_mature.txt.gz")
    return target, mirna_map, mirna_id
