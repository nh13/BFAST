#!/bin/sh

. test.definitions.sh
RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";

echo "      Running postprocessing.";

for PAIRED_END in 0 1
do
	for SPACE in 0 1
	do 
		if [ "$PAIRED_END" -eq "0" ]; then
			echo "        Testing -A "$SPACE;
		else
			echo "        Testing -A "$SPACE" -2";
		fi

		ALIGN=$OUTPUT_DIR"bfast.aligned.file.$OUTPUT_ID.$SPACE.$PAIRED_END.baf";

		# Run local alignment
		if [ "$PAIRED_END" -eq "0" ]; then
			CMD="../bpostprocess/bpostprocess -r $RG -i $ALIGN -a 3 -o $OUTPUT_ID.$SPACE.$PAIRED_END -d $OUTPUT_DIR";
		else
			CMD="../bpostprocess/bpostprocess -r $RG -i $ALIGN -a 3 -2 -o $OUTPUT_ID.$SPACE.$PAIRED_END -d $OUTPUT_DIR";
		fi
		$CMD 2> /dev/null > /dev/null; 
		# Get return code
		if [ "$?" -ne "0" ]; then
			# Run again without piping anything
			$CMD;
			exit 1
		fi
	done
done

# Test passed!
echo "      Postprocessing complete.";
exit 0
