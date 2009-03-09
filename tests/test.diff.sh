#!/bin/sh

. test.definitions.sh

echo "      Double-checking output files.";

CMD="md5sum -c $OUTPUT_DIR/tests.md5";
$CMD 2> /dev/null > /dev/null;
# Get return code
if [ "$?" -ne "0" ]; then
	# Run again without piping anything
	echo $CMD;
	$CMD
	exit 1
fi

echo "      Output files are the same.";
