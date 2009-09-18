#!/bin/sh

. test.definitions.sh

# Run test suite 
sh test.initialize.sh
sh test.fasta2brg.sh
sh test.index.sh
sh test.match.sh
sh test.localalign.sh
sh test.postprocess.sh
# Update md5sum
md5sum $OUTPUT_DIR/bfast* > $OUTPUT_DIR/tests.md5
# Archive
cd $DATA_DIR
tar -zxvf data.tar.gz
rm data.tar.gz
mv ../$OUTPUT_DIR/tests.md5 .
tar -cf data.tar *
gzip -9 --force data.tar
cd ..
#
sh test.cleanup.sh
