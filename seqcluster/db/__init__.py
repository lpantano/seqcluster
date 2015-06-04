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


def _insert_data(con, data):
    """
    insert line for each cluster
    """
    with con:
        cur = con.cursor()
        cur.execute("DROP TABLE IF EXISTS clusters;")
        cur.execute("CREATE TABLE clusters(Id INT, Locus TEXT, Annotation TEXT, Sequences TEXT)")
        for c in data[0]:
            locus = json.dumps(data[0][c]['loci'])
            annotation = json.dumps(data[0][c]['ann'])
            sequences = json.dumps(data[0][c]['seqs'])
            cur.execute("INSERT INTO clusters VALUES(%s, '%s', '%s', '%s')" % (c, locus, annotation, sequences))


def _close(con):
    if con:
        con.close()


def make_database(data, name="seqcluster.db", out_dir="database"):
    out_dir = safe_makedir(out_dir)
    op.abspath(out_dir)
    con = _create_db(op.join(out_dir, name))
    _insert_data(con, data)
    _close(con)

