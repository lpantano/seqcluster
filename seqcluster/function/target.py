"""Get enrichment of targets doing permutations for significance"""
from __future__ import print_function
import os.path as op
from collections import defaultdict
import gzip

import seqcluster.libs.logger as mylog

logger = mylog.getLogger(__name__)

def targets_enrichment(args):
    if args.sps not in ["hsa", "mmu"]:
        raise ValueError("Species not supported yet.")

    target, mirna_map, mirna_id = _get_files(args.annotation, args.sps)
    mirna_id_dict = _get_mirbase_id(mirna_id, args.sps)
    mirna_map_dict = _get_target_id(mirna_map, args.sps)
    mirna_list = _get_mirna_input(args.input)
    res = defaultdict(list)
    collect = defaultdict(list)
    with open(target) as in_handle:
        in_handle.readline()
        for line in in_handle:
            cols = line.split()
            mit_id = mirna_map_dict[cols[2]]
            for mit in mit_id:
                if mit in mirna_id_dict:
                    mirbase_id = mirna_id_dict[mit]
                    if mirbase_id in mirna_list:
                        info = map(str, [cols[0], cols[1], cols[4], cols[-3], cols[-1]])
                        keep = "\t".join(info)
                        if cols[-3] == "NULL":
                            cols[-3] = 0
                        if cols[-1] == "NULL":
                            cols[-1] = -1
                        score = float(cols[-3])
                        pt = float(cols[-1]) > 0.1 or cols[-1] == -1
                        if int(cols[4]) > 0 and pt:
                            if score < -0.2:
                                collect[mirbase_id].append({'info': keep, 'score': float(cols[-3])})
    with open(op.join(args.out, "matrix.tsv"), 'w') as ma_handle:
        for mir in collect:
            sorted_targets = sorted(collect[mir], key=lambda k: k['score'])
            for e in range(0, len(sorted_targets)):
                line = sorted_targets[e]['info']
                cols = line.split()
                name = cols[0].split(".")[0]
                res[(name, 'info')] = line.strip()
                res[(name, 'mirs')].append(mir)
                print("%s\t%s\t%s" % (name, mir, sorted_targets[e]['score']), file=ma_handle, end="")

    with open(op.join(args.out, "pairs.tsv"), 'w') as out_handle:
        print("mirs\tcounts\tgene\ttranscript\tgene\tsites\tcontext_score\tAggr_PCT",file=out_handle, end="")
        for gene, field in res:
            if field == "info":
                counts = len(res[(gene, 'mirs')])
                mirs = ",".join(list(set(res[(gene, 'mirs')])))
                print("%s\t%s\t%s\t%s" %  (mirs, counts,  gene, res[(gene, 'info')]), file=out_handle, end="")


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


def _get_files(path, sps):
    out = []
    fns = ["Summary_Counts.default_predictions.txt", "miR_Family_Info.txt", "mirna_mature.txt.gz"]
    for fn in fns:
        fn = op.join(path, fn)
        out.append(fn)
        if not op.exists(fn):
            raise IOError("Files not installed, please use this command:"
                          "seqcluster install --data path --genomes mm10")
    return out
