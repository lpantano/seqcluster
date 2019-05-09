"""
Modules to support sqlite implementation
"""
import os.path as op
import sqlite3 as lite
import logging
import json

from collections import defaultdict

from seqcluster.libs.utils import safe_dirs

logger = logging.getLogger('report')

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
    seqs = [list(s.values())[0] for s in cluster['seqs']]
    freqs = [list(f.values())[0] for f in cluster['freq']]
    data = []
    total_freq = {}

    for s, f in zip(seqs, freqs):
        fix = dict(zip(list(f.keys()), list(f.values())))
        data.append({'name': s, 'freq': fix})
        total_freq[s] = 1.0 * sum(list(fix.values())) / len(list(fix.values()))

    if len(total_freq) > 100:
        counts_50 = sorted(list(total_freq.values()))[-100]
        data = [e for e in data if 1.0 * sum(e['freq'].values()) / len(list(e['freq'].values())) > counts_50]
    return data

def _take_closest(num,collection):
    return min(collection, key=lambda x:abs(x-num))

def _get_closer(dat, pos):
    if pos in dat:
        return pos
    else:
        closest_pos = _take_closest(pos, dat.keys())
        if abs(closest_pos - pos) < 3:
            return closest_pos

def _set_format(profile):
    """
    Prepare dict to list of y values with same x
    """
    x = set()
    for sample in profile:
        x = x.union(set(profile[sample].keys()))
    if not x:
        return ''
    end, start = max(x), min(x)
    x = range(start, end, 4)
    scaled_profile = defaultdict(list)
    for pos in x:
        for sample in profile:
            y = _get_closer(profile[sample], pos)
            if y:
                scaled_profile[sample].append(profile[sample][y])
            else:
                scaled_profile[sample].append(0)
    return {'x': list(x), 'y': scaled_profile, 'names': list(scaled_profile.keys())}

def _insert_data(con, data):
    """
    insert line for each cluster
    """
    n = 0
    with con:
        cur = con.cursor()
        cur.execute("DROP TABLE IF EXISTS clusters;")
        cur.execute("CREATE TABLE clusters(Id INT, Description TEXT, Locus TEXT, Annotation TEXT, Sequences TEXT, Profile TXT, Precursor TXT)")
        for c in data[0]:
            n += 1
            locus = json.dumps(data[0][c]['loci'])
            annotation = json.dumps(data[0][c]['ann'])
            description = _get_description(data[0][c]['ann'])
            sequences = json.dumps(_get_sequences(data[0][c]))
            profile = "Not available."
            if 'profile' in data[0][c]:
                profile = json.dumps(_set_format(data[0][c]['profile']))
            precursor = json.dumps(data[0][c].get('precursor'))
            statement = "INSERT INTO clusters VALUES (%s, '%s', '%s', '%s', '%s', '%s', '%s')" % (c, description, locus, annotation, sequences, profile, precursor)
            cur.execute(statement)

    logger.info("Clusters inserted: %s" % n)

def _close(con):
    if con:
        con.close()

def make_database(data, name="seqcluster.db", out_dir="database"):
    out_dir = safe_dirs(out_dir)
    op.abspath(out_dir)
    con = _create_db(op.join(out_dir, name))
    _insert_data(con, data)
    _close(con)
