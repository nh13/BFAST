#!/bin/sh 
. test.definitions.sh
LAYOUT_FILE=$OUTPUT_DIR"layouts.txt";

echo "      Building an index.";

for SPACE in 0 1
do
	echo "        Testing -A "$SPACE;
	if [ "$SPACE" -eq "0" ]; then
		OUTPUT_ID=$OUTPUT_ID_NT;
	else
		OUTPUT_ID=$OUTPUT_ID_CS;
	fi
	RG_FASTA=$OUTPUT_DIR$OUTPUT_ID".fa";

	# Make an index
	CMD=$CMD_PREFIX"../bfast index -f $RG_FASTA -A $SPACE -m 1111111111 -w 6 -d 1 -i 1 -T $TMP_DIR";
	$CMD 2> /dev/null > /dev/null; 

	# Get return code
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		echo $CMD;
		$CMD;
		exit 1
	fi
done

# Test passed!
echo "      Index successfully built.";
exit 0
