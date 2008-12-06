#!/bin/sh

OUTPUT_ID="phi.x.174";
OUTPUT_DIR="output/"
TMP_DIR="tmp/";
OFFSETS=$OUTPUT_DIR"offsets.0-100.txt";

echo "      Finding matches.";

for SPACE in 0 1
do
	echo "        Testing -A "$SPACE;

	RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.$SPACE.brg";
	MAIN=$OUTPUT_DIR"main.indexes.$OUTPUT_ID.$SPACE.txt";
	SECONDARY=$OUTPUT_DIR"secondary.indexes.$OUTPUT_ID.$SPACE.txt";
	READS=$OUTPUT_DIR"reads.$OUTPUT_ID.$SPACE.fa";

	# Make files holding paths to index files
	ls -1 $OUTPUT_DIR/*$OUTPUT_ID.$SPACE*bif > $MAIN;
	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi
	touch $SECONDARY;
	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi

	# Find matches
	../bmatches/bmatches -r $RG -i $MAIN -I $SECONDARY -O $OFFSETS -R $READS -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR 2> /dev/null > /dev/null; 
	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi
done

# Test passed!
echo "      Matches found.";
exit 0
