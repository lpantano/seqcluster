
def read_gtf_line(cols):
    """parse gtf line to get class/name information"""
    try:
        group = cols[2]
        attrs = cols[8].split("; ")
        name = [attr.split(" ")[1] for attr in attrs if attr.split(" ")[0].endswith("name")]
        c = cols[0]
        s = int(cols[3])
        e = int(cols[4])
        st = cols[6]
        return [c, s, e, st, group, name[0]]
    except Exception, e:
        print "File is not in correct format"
        print "Expect chr source feature start end . strand attributes"
        print "Attributes are 'gene_name SNCA; gene_id ENSG; '"
        print "at least should contains '; *name NAME; '"
        print e
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
