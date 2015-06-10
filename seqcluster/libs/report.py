import logging
import string
import os
# from collections import Counter

import matplotlib
matplotlib.use('Agg')
import pylab
pylab.rcParams['figure.figsize'] = (25.0, 10.0)
from collections import Counter, defaultdict
import pandas as pd

from read import map_to_precursors
from utils import safe_dirs
from progressbar import ProgressBar

from bcbio.utils import file_exists

from seqcluster.html import HTML
from seqcluster import templates

logger = logging.getLogger('html')


def _get_link(c):
    """Gives html link tag for cluster link information"""
    return "<a href=%s/maps.html>%s</a>" % (c, c)


def _get_ann(dbs, features):
    """
    Gives format to annotation for html table output
    """
    value = ""
    for db, feature in zip(dbs, features):
        value += db + ":" + feature
    return value


def make_profile(data, out_dir, args):
    """
    Make html for each cluster
    """
    main_table = []
    header = ['id', 'ann']
    html_file = os.path.join(out_dir, "index.html")
    n = len(data[0])
    with ProgressBar(maxval=n, redirect_stdout=True) as p:
        for itern, c in enumerate(data[0]):
            p.update(itern)
            logger.debug("creating cluser: {}".format(c))
            safe_dirs(os.path.join(out_dir, c))
            valid, ann = _single_cluster(c, data, os.path.join(out_dir, c, "maps.tsv"), args)
            if valid:
                main_table.append([_get_link(c), _get_ann(valid, ann)])

    main_html = HTML.table(main_table, header_row=header, attribs={'id': 'keywords'})
    html_template = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(templates.__file__)), "main"))
    content = open(html_template).read()
    data = {'main': main_html}
    out_content = string.Template(content).safe_substitute(data)
    with open(html_file, 'w') as out_handle:
        print >>out_handle, out_content


def _expand(dat, counts, start, end):
    """
    expand the same counts from start to end
    """
    for pos in range(start, end):
        dat[pos] += counts
    return dat


def _convert_to_df(in_file):
    """
    convert data frame into table with pandas
    """
    dat = Counter()
    with open(in_file) as in_handle:
        for line in in_handle:
            cols = line.strip().split(" ")
            counts = float(cols[1].replace("cx", ""))
            dat = _expand(dat, counts, int(cols[4]), int(cols[5]))
    dat = {'positions': dat.keys(), 'expression': dat.values()}
    df = pd.DataFrame(data=dat, columns=dat.keys())
    df.set_index('positions', inplace=True)
    return df


def _make_html(c, html_file, figure_file, prefix):
    """
    create html from template, adding figure,
    annotation and sequences counts
    """
    ann = defaultdict(list)
    seqs_table = []

    src_img = "<img src=\"%s\" width=\"800\" height=\"350\" />" % os.path.basename(figure_file)
    coor_list = [" ".join(map(str, l)) for l in c['loci']]

    for pos in c['ann']:
        for db in pos:
            ann[db] += list(pos[db])
    logger.debug(ann)

    valid = [l for l in c['valid']]
    ann_list = [", ".join(list(set(ann[feature]))) for feature in ann if feature in valid]

    seqs = [s.values()[0] for s in c['seqs']]
    freq = [map(float, s.values()[0].values()) for s in c['freq']]
    header = ['seq'] + c['freq'][0].values()[0].keys()
    for s, f in zip(seqs, freq):
        f = map(round, f)
        seqs_table.append([s] + map(str, f))
    # seqs_html = seqs_html.replace("TABLE", "TABLE id=\"keywords\"")
    if not file_exists(html_file):
        coor_html = HTML.list(coor_list)
        ann_html = HTML.list(ann_list)
        seqs_html = HTML.table(seqs_table,
                               header_row=header, attribs={'id': 'keywords'})
        html_template = os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(templates.__file__)), "cluster"))
        content = open(html_template).read()
        data = {'profile': src_img,
                'loci': coor_html,
                'annotation': ann_html,
                'table': seqs_html}
        out_content = string.Template(content).safe_substitute(data)
        with open(html_file, 'w') as out_handle:
            print >>out_handle, out_content

    return valid, ann_list


def _single_cluster(c, data, out_file, args):
    """
    Map sequences on precursors and create
    expression profile
    """
    valid, ann = 0, 0
    figure_file = out_file.replace(".tsv", ".png")
    html_file = out_file.replace(".tsv", ".html")
    prefix = os.path.dirname(out_file)
    names = [round(sum(s.values()[0].values())) for s in data[0][c]['freq']]
    seqs = [s.values()[0] for s in data[0][c]['seqs']]

    loci = data[0][c]['loci']

    if loci[0][3] - loci[0][2] > 500:
        logger.info("locus bigger > 500 nt, skipping: %s" % loci)
        return valid, ann
    if not file_exists(out_file):
        logger.debug("map all sequences to all loci %s " % loci)
        map_to_precursors(seqs, names, {loci[0][0]: [loci[0][0:5]]}, out_file, args)
    # map_sequences_w_bowtie(sequences, precursors)

    logger.debug("plot sequences on loci")
    df = _convert_to_df(out_file)
    if not df.empty:
        if not file_exists(figure_file):
            plot = df.plot()
            plot.set_ylabel('Normalized expression', fontsize=25)
            plot.set_xlabel('Position', fontsize=25)
            plot.tick_params(axis='both', which='major', labelsize=20)
            plot.tick_params(axis='both', which='minor', labelsize=20)
            plot.get_figure().savefig(figure_file)
        valid, ann = _make_html(data[0][c], html_file, figure_file, prefix)

    return valid, ann
