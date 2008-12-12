#!/bin/sh

# This is the commands to make the simulated reads for the test

SAVE_DIR="save/";
OUTPUT_DIR="output/";
OUTPUT_ID="phi.x.174";
DATA_DIR="data/";
READ_LENGTH=50;
NUM_READS=10;
RG=$OUTPUT_DIR"bfast.rg.file.$OUTPUT_ID.0.brg";

echo "Creating reads";
# Initialize
echo "Initializing.";
for CMD in "sh test.initialize.sh" "sh test.bpreprocess.1.$OUTPUT_ID.sh"
do
	$CMD 2> /dev/null > /dev/null;
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		$CMD;
		exit 1
	fi
done

# Extract data since we wish to overwrite
cd $DATA_DIR
tar -zxvf data.tar.gz 2> /dev/null > /dev/null;
cd ..

# Create reads
echo "Creating reads."
for SPACE in 0 1
do
	for PAIRED_END in 0 1
	do
		for((NUM_ERRORS=0;NUM_ERRORS<=10;NUM_ERRORS++))
		do
			# No indels
			CMD="../butil/bgeneratereads $RG $SPACE 0 0 0 0 $NUM_ERRORS $READ_LENGTH $PAIRED_END 50 $NUM_READS";
			$CMD 2> /dev/null > /dev/null;
			if [ "$?" -ne "0" ]; then
				# Run again without piping anything
				$CMD;
				exit 1
			fi

			# Indels
			for((INDEL_LENGTH=1;INDEL_LENGTH<=10;INDEL_LENGTH++))
			do
				# Deletions
				CMD="../butil/bgeneratereads $RG $SPACE 1 $INDEL_LENGTH 0 0 $NUM_ERRORS $READ_LENGTH $PAIRED_END 50 $NUM_READS";
				$CMD 2> /dev/null > /dev/null;
				if [ "$?" -ne "0" ]; then
					# Run again without piping anything
					$CMD;
					exit 1
				fi
				# Insertions
				CMD="../butil/bgeneratereads $RG $SPACE 2 $INDEL_LENGTH 1 0 $NUM_ERRORS $READ_LENGTH $PAIRED_END 50 $NUM_READS";
				$CMD 2> /dev/null > /dev/null;
				if [ "$?" -ne "0" ]; then
					# Run again without piping anything
					$CMD;
					exit 1
				fi
			done
		done
		cat reads.*fa 2> /dev/null > $DATA_DIR"reads.$OUTPUT_ID.$SPACE.$PAIRED_END.fa";
		rm reads.*fa 2> /dev/null > /dev/null;
	done
done

# Archive data 
echo "Archiving data";
cd $DATA_DIR
cp *fa ../$OUTPUT_DIR. 2> /dev/null > /dev/null;
tar -cf data.tar *txt *fa 2> /dev/null > /dev/null;
gzip -9 --force data.tar 2> /dev/null > /dev/null;
rm *fa *txt 2> /dev/null > /dev/null;
cd ..

# Re-run all
echo "Re-running bpreprocess.2, bmatches, balign, and bpostprocess";
for CMD in "sh test.bpreprocess.2.$OUTPUT_ID.sh" "sh test.bmatches.$OUTPUT_ID.sh" "sh test.balign.$OUTPUT_ID.sh" "sh test.bpostprocess.$OUTPUT_ID.sh"
do
	$CMD 2> /dev/null > /dev/null;
	if [ "$?" -ne "0" ]; then
		# Run again without piping anything
		$CMD;
		exit 1
	fi
done

# Archive output
echo "Archiving output";
cd $OUTPUT_DIR
tar -cf save.tar bfast* 2> /dev/null > /dev/null;
gzip -9 --force save.tar 2> /dev/null > /dev/null;
mv save.tar.gz ../$SAVE_DIR. 2> /dev/null > /dev/null;
cd ..

# Clean up
echo "Cleaning up";
sh test.cleanup.sh
echo "Reads created.";
