#!/bin/sh

OUTPUT_ID="phi.x.174";
OUTPUT_DIR="output/"
TMP_DIR="tmp/";
RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";

echo "      Running local alignment.";

for SPACE in 0 1
do
	echo "        Testing -A "$SPACE;

	MATCHES=$OUTPUT_DIR"bfast.matches.file.$OUTPUT_ID.$SPACE.-1.-1.0.0.0.0.0.0.bmf";
	SCORING=$OUTPUT_DIR"scoring.$SPACE.txt";

	# Run local alignment
	../balign/balign -r $RG -m $MATCHES -x $SCORING -A $SPACE -o $OUTPUT_ID -d $OUTPUT_DIR -T $TMP_DIR 2> /dev/null > /dev/null; 
	# Get return code
	if [ "$?" -ne "0" ]; then
		exit 1
	fi
done

# Test passed!
echo "      Local alignment complete.";
exit 0
