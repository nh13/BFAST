#!/bin/sh

. test.definitions.sh

echo "      Cleaning up files.";

rm -r $OUTPUT_DIR $TMP_DIR 

# Get return code
if [ "$?" -ne "0" ]; then
	exit 1
fi

# Test passed!
echo "      Files cleaned up.";
exit 0
