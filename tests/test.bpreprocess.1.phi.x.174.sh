#!/bin/sh

OUTPUT_ID="phi.x.174";
OUTPUT_DIR="output/";
RG_FASTA=$OUTPUT_DIR$OUTPUT_ID".fa";
TMP_DIR="tmp/";

echo "      Building a reference genome.";

for SPACE in 0 1
do  
	echo "        Testing -A "$SPACE;
	# Make reference genome
	../bpreprocess/bpreprocess -r $RG_FASTA -a 0 -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR 2> /dev/null > /dev/null;

	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi
done

# Test passed!
echo "      Reference genome successfully built.";
exit 0
