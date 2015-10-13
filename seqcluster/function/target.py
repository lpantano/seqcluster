"""Get enrichment of targets doing permutations for significance"""
import os
import gzip

from bcbio import utils


def targets_enrichment(args):
    if args.species not in ["human"]:
        raise ValueError("Species not supported yet.")

    target, mirna_map, mirna_id = _get_file(args.annotation, args.species)
    mirna_id_dict = _get_mirbase_id(mirna_id, args.species)
    mirna_map_dict = _get_target_id(mirna_map)

def _get_target_id(fn,sps):
    mapping = dict()
    with gzip.open(fn, 'rb') as in_handle:
        for line in in_handle:
            cols = line.split()
            mapping[cols[0]] = cols[-1].strip()
    return mapping


def _get_mirbase_id(fn,sps):
    mapping = dict()
    with gzip.open(fn, 'rb') as in_handle:
        for line in in_handle:
            if cols[1].startswith(sps):
                cols = line.split()
                mapping[cols[3]] = cols[2]
    return mapping


def _get_files(path):
    target = op.join(path, "Conserved_Family_Conserved_Targets_Info.txt.zip")
    mirna_map = op.join(path, "miR_Family_info.txt.zip")
    mirna_id = op.join(path, "mirna_mature.txt.gz")
    return target, mirna_map, mirna_id
