"""
Modules to support sqlite implementation
"""
import os.path as op
import sqlite3 as lite
import json

from bcbio.utils import safe_makedir


def _create_db(name):
    """
    creater connection to sqlite
    """
    con = lite.connect(name)
    return con


def _get_description(string):
    """
    Parse annotation to get nice description
    """
    ann = set()
    if not string:
        return "This cluster is inter-genic."
    for item in string:
        for db in item:
            ann = ann.union(set(item[db]))
    return "annotated as: %s ..." % ",".join(list(ann)[:3])

def _get_sequences(cluster):
    seqs = [s.values()[0] for s in cluster['seqs']]
    freqs = [f.values()[0] for f in cluster['freq']]
    data = []
    for s, f in zip(seqs, freqs):
        data.append({'name': s, 'freq': f})
    return data

def _insert_data(con, data):
    """
    insert line for each cluster
    """
    with con:
        cur = con.cursor()
        cur.execute("DROP TABLE IF EXISTS clusters;")
        cur.execute("CREATE TABLE clusters(Id INT, Description TEXT, Locus TEXT, Annotation TEXT, Sequences TEXT, Profile TXT)")
        for c in data[0]:
            locus = json.dumps(data[0][c]['loci'])
            annotation = json.dumps(data[0][c]['ann'])
            description = _get_description(data[0][c]['ann'])
            sequences = json.dumps(_get_sequences(data[0][c]))
            keys = data[0][c]['freq'][0].values()[0].keys()
            profile = "Not available."
            if data[0][c]['profile']:
                data[0][c]['profile']['num_lines'] = len(data[0][c]['profile'])
                profile = json.dumps(data[0][c]['profile'])
            cur.execute("INSERT INTO clusters VALUES(%s, '%s', '%s', '%s', '%s', '%s')" % (c, description, locus, annotation, sequences, profile))


def _close(con):
    if con:
        con.close()


def make_database(data, name="seqcluster.db", out_dir="database"):
    out_dir = safe_makedir(out_dir)
    op.abspath(out_dir)
    con = _create_db(op.join(out_dir, name))
    _insert_data(con, data)
    _close(con)

