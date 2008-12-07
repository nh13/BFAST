#!/bin/sh

OUTPUT_ID="phi.x.174";
OUTPUT_DIR="output/"
RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";

echo "      Running postprocessing.";

for SPACE in 0 1
do 
	echo "        Testing -A "$SPACE;

	ALIGN=$OUTPUT_DIR"bfast.aligned.file.$OUTPUT_ID.$SPACE.baf";
	# Run local alignment
	CMD="../bpostprocess/bpostprocess -r $RG -i $ALIGN -o $OUTPUT_ID.$SPACE -d $OUTPUT_DIR";
	$CMD 2> /dev/null > /dev/null; 
	# Get return code
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		$CMD;
		exit 1
	fi
done

# Test passed!
echo "      Postprocessing complete.";
exit 0
