#!/bin/sh

SAVE_DIR="save/";
OUTPUT_DIR="output/";
TMP_DIR="tmp/";

echo "      Cleaning up files.";

rm -r $OUTPUT_DIR $TMP_DIR $SAVE_DIR/bfast*

# Get return code
if [ "$?" -ne "0" ]; then
	exit 1
fi

# Test passed!
echo "      Files cleaned up.";
exit 0
