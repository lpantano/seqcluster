import logger as mylog
import os

from seqcluster.libs.classes import annotation, dbannotation

logger = mylog.getLogger("run")


def read_gtf_line(cols):
    """parse gtf line to get class/name information"""
    try:
        group = cols[2]
        attrs = cols[8].split(";")
        name = [attr.strip().split(" ")[1] for attr in attrs if attr.strip().split(" ")[0].lower().endswith("name")]
        biotype = [attr.strip().split(" ")[1] for attr in attrs if attr.strip().split(" ")[0].lower().endswith("biotype")]
        if biotype:
            group = biotype[0]
        c = cols[0]
        s = int(cols[3])
        e = int(cols[4])
        st = cols[6]
        return [c, s, e, st, group, name[0]]
    except Exception, e:
        logger.error(cols)
        logger.error("File is not in correct format")
        logger.error("Expect chr source feature start end . strand attributes")
        logger.error("Attributes are 'gene_name SNCA; gene_id ENSG; '")
        logger.error("The 3rd column is used as type of small RNA (like miRNA)")
        logger.error("at least should contains '; *name NAME; '")
        logger.error(e)
        raise


def _position_in_feature(pos_a, pos_b):
    """return distance to 3' and 5' end of the feature"""
    strd = "-"
    if pos_a[2] in pos_b[2]:
        strd = "+"
    if pos_a[2] in "+" and pos_b[2] in "+":
        lento5 = pos_a[0] - pos_b[1] + 1
        lento3 = pos_a[1] - pos_b[1] + 1
    if pos_a[2] in "+" and pos_b[2] in "-":
        lento5 = pos_a[1] - pos_b[0] + 1
        lento3 = pos_a[0] - pos_b[1] + 1
    if pos_a[2] in "-" and pos_b[2] in "+":
        lento5 = pos_a[0] - pos_b[1] + 1
        lento3 = pos_a[1] - pos_b[0] + 1
    if pos_a[2] in "-" and pos_b[2] in "-":
        lento3 = pos_a[0] - pos_b[0] + 1
        lento5 = pos_a[1] - pos_b[1] + 1
    else:
        lento5 = pos_a[0] - pos_b[0] + 1
        lento3 = pos_a[1] - pos_b[1] + 1
    return lento5, lento3, strd


def anncluster(c, clus_obj, db, type_ann):
    """intersect transcription position with annotation files"""
    id_sa, id_ea, id_id, id_idl, id_sta = 1, 2, 3, 4, 5
    if type_ann == "bed":
        id_sb = 7
        id_eb = 8
        id_stb = 11
        id_tag = 9
    ida = 0
    clus_id = clus_obj.clus
    loci_id = clus_obj.loci
    db = os.path.splitext(db)[0]
    logger.debug("Type:%s\n" % type_ann)
    for cols in c.features():
        if type_ann == "gtf":
            cb, sb, eb, stb, db, tag = read_gtf_line(cols[6:])
        else:
            sb = int(cols[id_sb])
            eb = int(cols[id_eb])
            stb = cols[id_stb]
            tag = cols[id_tag]
        id = int(cols[id_id])
        idl = int(cols[id_idl])
        if (id in clus_id):
            clus = clus_id[id]
            sa = int(cols[id_sa])
            ea = int(cols[id_ea])
            ida += 1
            lento5, lento3, strd = _position_in_feature([sa, ea, cols[id_sta]], [sb, eb, stb])
            if db in loci_id[idl].db_ann:
                ann = annotation(db, tag, strd, lento5, lento3)
                tdb = loci_id[idl].db_ann[db]
                tdb.add_db_ann(ida, ann)
                loci_id[idl].add_db(db, tdb)
            else:
                ann = annotation(db, tag, strd, lento5, lento3)
                tdb = dbannotation(1)
                tdb.add_db_ann(ida, ann)
                loci_id[idl].add_db(db, tdb)
            clus_id[id] = clus
    clus_obj.clus = clus_id
    clus_obj.loci = loci_id
    return clus_obj


