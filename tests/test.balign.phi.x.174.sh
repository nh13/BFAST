#!/bin/sh

. test.definitions.sh
RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";

echo "      Running local alignment.";

for PAIRED_END in 0 1
do
	for SPACE in 0 1
	do
		if [ "$PAIRED_END" -eq "0" ]; then
			echo "        Testing -A "$SPACE;
		else
			echo "        Testing -A "$SPACE" -2";
		fi

		SCORING=$OUTPUT_DIR"scoring.$SPACE.txt";
		MATCHES=$OUTPUT_DIR"bfast.matches.file.$OUTPUT_ID.$SPACE.$PAIRED_END.bmf";

		# Run local alignment
		if [ "$PAIRED_END" -eq "0" ]; then
			CMD="../balign/balign -r $RG -m $MATCHES -x $SCORING -A $SPACE -X $SPACE -O 15 -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
		else
			CMD="../balign/balign -r $RG -m $MATCHES -x $SCORING -A $SPACE -X $SPACE -2 -O 15 -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
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
echo "      Local alignment complete.";
exit 0
