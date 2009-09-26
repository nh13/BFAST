#!/bin/sh

. test.definitions.sh

echo "      Cleaning up files.";

rm -r $OUTPUT_DIR $TMP_DIR;
ls -1 $DATA_DIR/* | grep -v bz2 | xargs rm;

# Get return code
if [ "$?" -ne "0" ]; then
	exit 1
fi

# Test passed!
echo "      Files cleaned up.";
exit 0
