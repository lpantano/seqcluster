#!/usr/bin/env bash
# bash script inspared from bcbio-nextgen: https://github.com/chapmanb/bcbio-nextgen/blob/master/tests/run_tests.sh
# It allows to run specific test using the name like
# ./run_test.sh test_align
set -o pipefail  # trace ERR through pipes
set -o errtrace  # trace ERR through 'time command' and other functions
set -o nounset   ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1; }
SC_DIR=$(dirname "$(readlinkf `which seqcluster`)")
unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1
echo $SC_DIR

echo "Run module test."
"$SC_DIR/nosetests" -v -s -a "$@"
