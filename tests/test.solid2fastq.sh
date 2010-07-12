#!/bin/sh
# 
. test.definitions.sh

run()
{
  # Get return code
	eval $CMD 2> /dev/null
  if [ "$?" -ne "0" ]; then
	  # Run again without piping anything
	  echo $CMD
	  eval $CMD
	  exit 1
  fi
}

echo "      Running s2f"
#echo "Current dir: `pwd`"
# BF mode
echo "        Testing F3 5"
CMD="${S2F} -n 2 -o ${OUTPUT_DIR}bfast.s2f.se.5 ${OUTPUT_DIR}s2f.5.F3.csfasta ${OUTPUT_DIR}s2f.5.F3_QV.qual"
run 
echo "        Testing F3 6"
CMD="${S2F} -n 2 -o ${OUTPUT_DIR}s2f.se.6 ${OUTPUT_DIR}s2f.6.F3.csfasta ${OUTPUT_DIR}s2f.6.F3_QV.qual"
run 
echo "        Testing V4"
CMD="$S2F -n 2 -o ${OUTPUT_DIR}s2f.se.v4 ${OUTPUT_DIR}s2f.v4.F3.csfasta ${OUTPUT_DIR}s2f.v4.R3.csfasta \
${OUTPUT_DIR}s2f.v4.F3_QV.qual ${OUTPUT_DIR}s2f.v4.R3_QV.qual"
run 
# bwa mode
echo "        Testing BWA single end"
CMD="$S2F -b -o ${OUTPUT_DIR}s2f.bwa.one_end  ${OUTPUT_DIR}s2f.5.F3.csfasta ${OUTPUT_DIR}s2f.5.F3_QV.qual"
run 
echo "        Testing BWA two ends"
CMD="$S2F -b -o ${OUTPUT_DIR}s2f.bwa ${OUTPUT_DIR}s2f.5.F3.csfasta ${OUTPUT_DIR}s2f.5.R3.csfasta \
${OUTPUT_DIR}s2f.5.F3_QV.qual ${OUTPUT_DIR}s2f.5.R3_QV.qual"
run 

cd ..

# Data generated
echo "      done"
exit 0
