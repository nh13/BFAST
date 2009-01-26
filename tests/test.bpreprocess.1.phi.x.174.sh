#!/bin/sh

. test.definitions.sh
RG_FASTA=$OUTPUT_DIR$OUTPUT_ID".fa";
TMP_DIR="tmp/";

echo "      Building a reference genome.";

for SPACE in 0 1
do  
	echo "        Testing -A "$SPACE;
	# Make reference genome
	CMD="../bpreprocess/bpreprocess -r $RG_FASTA -a 0 -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR";
	$CMD 2> /dev/null > /dev/null;

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
