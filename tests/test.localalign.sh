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
	RG_FASTA=$OUTPUT_DIR$OUTPUT_ID".fa";
	MATCHES=$OUTPUT_DIR"bfast.matches.file.$OUTPUT_ID.bmf";

	# Run local alignment
	CMD=$CMD_PREFIX"../bfast localalign -f $RG_FASTA -m $MATCHES -A $SPACE -o 15 -T $TMP_DIR > ${OUTPUT_DIR}bfast.aligned.file.$OUTPUT_ID.baf";
	eval $CMD 2> /dev/null;
	# Get return code
	if [ "$?" -ne "0" ]; then
		echo $CMD;
		# Run again without piping anything
		eval $CMD;
		exit 1
	fi
done

# Test passed!
echo "      Local alignment complete.";
exit 0
