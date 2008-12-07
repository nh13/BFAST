#!/bin/sh

OUTPUT_ID="phi.x.174";
OUTPUT_DIR="output/";
SAVE_DIR="save/";

echo "      Comparing output files.";
for PREFIX in rg index reads.filtered matches aligned not.aligned reported not.reported
do
	for SPACE in 0 1
	do
		NAME="bfast.$PREFIX.file.$OUTPUT_ID.$SPACE";
		echo "        Comparing $NAME*";

		diff -q $OUTPUT_DIR/$NAME* $SAVE_DIR/$NAME*;

		# Get return code
		if [ "$?" -ne "0" ]; then
			echo "        $NAME* did not match.";
			exit 1
		fi
	done
done
echo "      Output files are the same.";
