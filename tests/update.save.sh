#!/bin/sh

. test.definitions.sh

# Run test suite 
sh test.initialize.sh
sh test.bpreprocess.1.sh
sh test.bpreprocess.2.sh
sh test.bmatches.sh
sh test.balign.sh
sh test.bpostprocess.sh
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
