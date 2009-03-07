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

	RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.$SPACE.brg";
	MAIN=$OUTPUT_DIR"main.indexes.$OUTPUT_ID..txt";
	SECONDARY=$OUTPUT_DIR"secondary.indexes.$OUTPUT_ID..txt";
	READS=$OUTPUT_DIR"reads.$OUTPUT_ID.fastq";

	# Make files holding paths to index files
	ls -1 $OUTPUT_DIR/*$OUTPUT_ID.*bif > $MAIN;
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
	CMD="../bmatches/bmatches -r $RG -i $MAIN -I $SECONDARY -O $OFFSETS -R $READS -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
	$CMD 2> /dev/null > /dev/null; 
	# Get return code
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		$CMD;
		exit 1
	fi
done

# Test passed!
echo "      Matches found.";
exit 0
