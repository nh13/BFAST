#!/bin/sh

. test.definitions.sh

echo "      Running local alignment.";

for SPACE in 0 1
do
	echo "        Testing -A "$SPACE;
	if [ "$SPACE" -eq "0" ]; then
		OUTPUT_ID=$OUTPUT_ID_NT;
	else
		OUTPUT_ID=$OUTPUT_ID_CS;
	fi

	RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";
	PROFILE=$OUTPUT_DIR"bfast.indexes.profile.$OUTPUT_ID.txt";
	SCORING=$OUTPUT_DIR"scoring.$SPACE.txt";
	MATCHES=$OUTPUT_DIR"bfast.matches.file.$OUTPUT_ID.bmf";

	# Run local alignment
	CMD="../balign/balign -r $RG -m $MATCHES -x $SCORING -i $PROFILE -A $SPACE -O 15 -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
	$CMD 2> /dev/null > /dev/null; 
	# Get return code
	if [ "$?" -ne "0" ]; then
		echo $CMD;
		# Run again without piping anything
		$CMD;
		exit 1
	fi
done

# Test passed!
echo "      Local alignment complete.";
exit 0
