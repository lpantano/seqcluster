set -o pipefail

bcbio_nextgen.py -w template template.yaml mirqc_bcbio.csv raw/* --only-metadata
bcbio_nextgen.py -w template template.yaml non_mirqc_bcbio.csv raw/* --only-metadata


cd mirqc_bcbio/work
bcbio_nextgen.py ../config/mirqc_bcbio.yaml

cd non_mirqc_bcbio/work
bcbio_nextgen.py ../config/non_mirqc_bcbio.yaml
