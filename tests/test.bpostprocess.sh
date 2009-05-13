#!/bin/sh

. test.definitions.sh

echo "      Running postprocessing.";

for SPACE in 0 1
do 
	echo "        Testing -A "$SPACE;
	if [ "$SPACE" -eq "0" ]; then
		OUTPUT_ID=$OUTPUT_ID_NT;
	else
		OUTPUT_ID=$OUTPUT_ID_CS;
	fi
	RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";
	ALIGN=$OUTPUT_DIR"bfast.aligned.file.$OUTPUT_ID.baf";

	# Run local alignment
	CMD=$CMD_PREFIX"../bpostprocess/bpostprocess -r $RG -i $ALIGN -a 3 -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR";
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
echo "      Postprocessing complete.";
exit 0
