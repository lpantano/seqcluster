#!/usr/bin/env bash
# bash script inspared from bcbio-nextgen: https://github.com/chapmanb/bcbio-nextgen/blob/master/tests/run_tests.sh
# It allows to run specific test using the name like
# ./run_test.sh test_align
set -e

readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1; }
SC_DIR=$(dirname "$(readlinkf `which seqcluster`)")
unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1
echo $SC_DIR
"$SC_DIR/nosetests" -v -s -a "$@"
