#!/bin/sh 
. test.definitions.sh
LAYOUT_FILE=$OUTPUT_DIR"layouts.txt";

echo "      Building an index.";

for SPACE in 0 1
do
	echo "        Testing -A "$SPACE;
	if [ "$SPACE" -eq "0" ]; then
		OUTPUT_ID=$OUTPUT_ID_NT;
		RG="$OUTPUT_DIR$OUTPUT_ID.fa.nt.brg";
	else
		OUTPUT_ID=$OUTPUT_ID_CS;
		RG="$OUTPUT_DIR$OUTPUT_ID.fa.cs.brg";
	fi


	# Make an index
	CMD=$CMD_PREFIX"../bfast index -r $RG -A $SPACE -i $LAYOUT_FILE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
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
