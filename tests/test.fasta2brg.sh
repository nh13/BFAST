#!/bin/sh

. test.definitions.sh
TMP_DIR="tmp/";

echo "      Building a reference genome.";

for SPACE in 0 1
do  
	echo "        Testing -A "$SPACE;
	if [ "$SPACE" -eq "0" ]; then
		OUTPUT_ID=$OUTPUT_ID_NT;
	else
		OUTPUT_ID=$OUTPUT_ID_CS;
	fi
	RG_FASTA=$OUTPUT_DIR$OUTPUT_ID".fa";

	# Make reference genome in nt space always 
	CMD=$CMD_PREFIX"../bfast fasta2brg -f $RG_FASTA -A 0";
	eval $CMD 2> /dev/null;
	# Make reference genome in color space if necessary
	if [ "$SPACE" -eq "1" ]; then
		CMD=$CMD_PREFIX"../bfast fasta2brg -f $RG_FASTA -A $SPACE";
		eval $CMD 2> /dev/null;
	fi

	# Get return code
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		echo $CMD;
		eval $CMD;
		exit 1
	fi
done

# Test passed!
echo "      Reference genome successfully built.";
exit 0
