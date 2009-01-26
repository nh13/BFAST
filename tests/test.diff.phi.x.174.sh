#!/bin/sh

. test.definitions.sh

echo "      Double-checking output files.";
for PREFIX in rg index 
do
	for SPACE in 0 1
	do
		NAME="bfast.$PREFIX.file.$OUTPUT_ID.$SPACE";
		echo "        Comparing $NAME*";

		diff -q $OUTPUT_DIR/$NAME* $SAVE_DIR/$NAME*;

		# Get return code
		if [ "$?" -ne "0" ]; then
			echo "        $NAME* did not match.";
			exit 1
		fi
	done
done

for PREFIX in reads.filtered matches aligned not.aligned reported not.reported
do
	for PAIRED_END in 0 1
	do
		for SPACE in 0 1
		do
			NAME="bfast.$PREFIX.file.$OUTPUT_ID.$SPACE.$PAIRED_END";
			echo "        Comparing $NAME*";

			diff -q $OUTPUT_DIR/$NAME* $SAVE_DIR/$NAME*;

			# Get return code
			if [ "$?" -ne "0" ]; then
				echo "        $NAME* did not match.";
				exit 1
			fi
		done
	done
done
echo "      Output files are the same.";
