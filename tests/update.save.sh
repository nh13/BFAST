#!/bin/sh

SAVE_DIR="save/";
OUTPUT_ID="phi.x.174";

# We update the bfast reference genome and bfast inded files for each new package version
sh test.initialize.sh
sh test.bpreprocess.1.$OUTPUT_ID.sh
sh test.bpreprocess.2.$OUTPUT_ID.sh
# Only copy over the reference genome and index files
cp output/bfast.rg.file*brg $SAVE_DIR.
cp output/bfast.index.file*bif $SAVE_DIR.
# Archive
cd $SAVE_DIR
tar -cf save.tar bfast*
gzip -9 --force save.tar
cd ..
# package
#sh test.cleanup.sh
