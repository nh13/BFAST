#!/bin/sh

OUTPUT_ID="phi.x.174";
OUTPUT_DIR="output/"
TMP_DIR="tmp/";
OFFSETS=$OUTPUT_DIR"offsets.0-100.txt";

echo "      Finding matches.";

for PAIRED_END in 0 1
do
	for SPACE in 0 1
	do
		if [ "$PAIRED_END" -eq "0" ]; then
			echo "        Testing -A "$SPACE;
		else
			echo "        Testing -A "$SPACE" -2";
		fi

		RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.$SPACE.brg";
		MAIN=$OUTPUT_DIR"main.indexes.$OUTPUT_ID.$SPACE.txt";
		SECONDARY=$OUTPUT_DIR"secondary.indexes.$OUTPUT_ID.$SPACE.txt";
		READS=$OUTPUT_DIR"reads.$OUTPUT_ID.$SPACE.$PAIRED_END.fa";

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
		if [ "$PAIRED_END" -eq "0" ]; then
			CMD="../bmatches/bmatches -r $RG -i $MAIN -I $SECONDARY -O $OFFSETS -R $READS -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
		else
			CMD="../bmatches/bmatches -r $RG -i $MAIN -I $SECONDARY -O $OFFSETS -R $READS -A $SPACE -2 -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
		fi
		$CMD 2> /dev/null > /dev/null; 
		# Get return code
		if [ "$?" -ne "0" ]; then
			# Run again without piping anything
			$CMD;
			exit 1
		fi
	done
done

# Test passed!
echo "      Matches found.";
exit 0
