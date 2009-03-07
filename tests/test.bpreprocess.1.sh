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
	CMD="../bpreprocess/bpreprocess -r $RG_FASTA -a 0 -A 0 -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
	$CMD 2> /dev/null > /dev/null;
	# Make reference genome in color space if necessary
	if [ "$SPACE" -eq "1" ]; then
		CMD="../bpreprocess/bpreprocess -r $RG_FASTA -a 0 -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
		$CMD 2> /dev/null > /dev/null;
	fi

	# Get return code
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		$CMD;
		exit 1
	fi
done

# Test passed!
echo "      Reference genome successfully built.";
exit 0
