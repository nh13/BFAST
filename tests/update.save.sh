#!/bin/sh

# We update the bfast reference genome and bfast inded files for each new package version
sh test.initialize.sh
sh test.bpreprocess.1.phi.x.174.sh
sh test.bpreprocess.2.phi.x.174.sh
# Only copy over the reference genome and index files
cp output/bfast.rg.file*brg save/.
cp output/bfast.index.file*bif save/.
# Archive
cd save
tar -cf save.tar bfast*
gzip -9 --force save.tar
cd ..
# package
sh test.cleanup.sh
