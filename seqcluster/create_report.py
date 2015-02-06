import os
import logging
from libs.read import load_data
from libs.report import make_profile
from libs.utils import safe_dirs


logger = logging.getLogger('report')


def report(args):
    """
    Create report in html format
    """
    logger.info("reading sequeces")
    data = load_data(args.json)
    out_dir = os.path.join(args.out, "html")
    safe_dirs(out_dir)

    logger.info("create profile")
    make_profile(data, out_dir, args)
    # get_sequences_from_cluster()
    # get_matrix_position()
    # plot_sequences()
    logger.info("Done")
