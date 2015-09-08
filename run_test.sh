readlinkf(){ perl -MCwd -e 'print Cwd::abs_path shift' $1; }
SC_DIR=$(dirname "$(readlinkf `which seqcluster`)")
unset PYTHONHOME
unset PYTHONPATH
export PYTHONNOUSERSITE=1
echo $SC_DIR
"$SC_DIR/nosetests" -v -s -a "$@"
