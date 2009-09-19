#!/bin/sh

. test.definitions.sh
OFFSETS=$OUTPUT_DIR"offsets.0-100.txt";

echo "      Finding matches.";

for SPACE in 0 1
do
	echo "        Testing -A "$SPACE;
	if [ "$SPACE" -eq "0" ]; then
		OUTPUT_ID=$OUTPUT_ID_NT;
	else
		OUTPUT_ID=$OUTPUT_ID_CS;
	fi
	RG_FASTA=$OUTPUT_DIR$OUTPUT_ID".fa";
	MAIN=$OUTPUT_DIR"main.indexes.$OUTPUT_ID.txt";
	READS=$OUTPUT_DIR"reads.$OUTPUT_ID.fastq";

	# Make files holding paths to index files
	ls -1 $OUTPUT_DIR/*$OUTPUT_ID.*bif > $MAIN;
	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi
	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi

	# Find matches
	CMD="${CMD_PREFIX}../bfast match -f $RG_FASTA -r $READS -A $SPACE -T $TMP_DIR > ${OUTPUT_DIR}bfast.matches.file.$OUTPUT_ID.bmf";
	eval $CMD 2> /dev/null;

	# Get return code
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		echo $CMD;
		eval $CMD;
		exit 1
	fi
done

# Test passed!
echo "      Matches found.";
exit 0
