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
$MD5BIN $OUTPUT_DIR/bfast* > $OUTPUT_DIR/tests.md5
# Archive
cd $DATA_DIR
tar -jxvf data.tar.bz2
rm data.tar.bz2
mv ../$OUTPUT_DIR/tests.md5 .
tar -cf data.tar *
bzip2 -9 --force data.tar
cd ..
#
sh test.cleanup.sh
