#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>
#include <string.h>
#include "../blib/BError.h"

FILE *open_output_file(char*, int32_t, int32_t);
void print_fastq(FILE*, FILE*, FILE*);
char *get_read_name(FILE*, FILE*);
int32_t cmp_read_names(char*, char*);

int main(int argc, char *argv[])
{
	char *output_prefix=NULL;
	int32_t num_reads_per_file=-1;
	int32_t number_of_ends;
	char **csfasta_filenames=NULL;
	char **qual_filenames=NULL;
	FILE **fps_csfasta=NULL;
	FILE **fps_qual=NULL;
	FILE *fp_output=NULL;
	int32_t output_suffix_number;
	char c;
	int32_t i, j;

	// Get Parameters
	while((c = getopt(argc, argv, "n:o:")) >= 0) {
		switch(c) {
			case 'n':
				num_reads_per_file=atoi(optarg); break;
				break;
			case 'o':
				output_prefix=strdup(optarg); break;
				break;
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	// Print Usage
	if(argc == optind ||
			0 != ((argc - optind) % 2) ||
			(argc - optind) < 2) {
		fprintf(stderr, "solid2fastq %s\nCopyright 2009.\n",
				PACKAGE_VERSION);
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: solid2fastq [options] <list of .csfasta files> <lsit of .qual files>\n");
		fprintf(stderr, "\t-n\t\tnumber of reads per file.\n");
		fprintf(stderr, "\t-o\t\toutput prefix.\n");
		fprintf(stderr, "\n send bugs to %s\n", PACKAGE_BUGREPORT);
	}

	// Copy over the filenames
	number_of_ends = 0;
	for(i=optind;i<argc;i+=2) {
		number_of_ends++;
		// Reallocate memory
		csfasta_filenames = realloc(csfasta_filenames, sizeof(char*)*number_of_ends);
		if(NULL == csfasta_filenames) {
			PrintError("solid2fastq",
					"csfasta_filenames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		qual_filenames = realloc(qual_filenames, sizeof(char*)*number_of_ends);
		if(NULL == qual_filenames) {
			PrintError("solid2fastq",
					"qual_filenames",
					"Could not reallocate memory",
					Exit,
					ReallocMemory);
		}
		// Copy over filenames
		csfasta_filenames[number_of_ends-1] = strdup(argv[i]);
		qual_filenames[number_of_ends-1] = strdup(argv[i+1]);
	}

	// Allocate memory for input file pointers
	fps_csfasta = malloc(sizeof(FILE*)*number_of_ends);
	if(NULL == fps_csfasta) {
		PrintError("solid2fastq",
				"fps_csfasta",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	fps_qual = malloc(sizeof(FILE*)*number_of_ends);
	if(NULL == fps_qual) {
		PrintError("solid2fastq",
				"fps_qual",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	// Open input files
	for(i=0;i<number_of_ends;i++) {
		if(!(fps_csfasta[i] = fopen(csfasta_filenames[i], "r"))) {
			PrintError("solid2fastq",
					csfasta_filenames[i],
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
		if(!(fps_qual[i] = fopen(qual_filenames[i], "r"))) {
			PrintError("solid2fastq",
					qual_filenames[i],
					"Could not open file for reading",
					Exit,
					OpenFileError);
		}
	}

	// Skip over comments
	for(i=0;i<number_of_ends;i++) {
		fpos_t cur_pos;
		char cur_line[1024]="\0";

		assert(0 == fgetpos(fps_csfasta[i], &cur_pos));
		while(NULL != fgets(cur_line, 1024, fps_csfasta[i]) &&
				'#' == cur_line[0]) {
			assert(0 == fgetpos(fps_csfasta[i], &cur_pos));
		}
		assert(0 == fsetpos(fps_csfasta[i], &cur_pos));

		assert(0 == fgetpos(fps_qual[i], &cur_pos));
		while(NULL != fgets(cur_line, 1024, fps_qual[i]) &&
				'#' == cur_line[0]) {
			assert(0 == fgetpos(fps_qual[i], &cur_pos));
		}
		assert(0 == fsetpos(fps_qual[i], &cur_pos));
	}

	int32_t *to_print=NULL;
	int32_t more_fps_left=1;
	int32_t output_count=0;

	to_print = malloc(sizeof(int32_t)*number_of_ends);
	if(NULL == to_print) {
		PrintError("solid2fastq",
				"to_print",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	output_suffix_number = 1;
	more_fps_left=number_of_ends;
	output_count = 0;
	// Open output file
	fp_output = open_output_file(output_prefix, output_suffix_number, num_reads_per_file); 
	// TODO
	while(0 < more_fps_left) { // while an input file is still open
		// Get all with min read name
		char *min_read_name=NULL;
		for(i=0;i<number_of_ends;i++) {
			if(NULL != fps_csfasta[i] && NULL != fps_qual[i]) {
				fpos_t pos_csfasta, pos_qual;
				char *read_name=NULL;
				assert(0 == fgetpos(fps_csfasta[i], &pos_csfasta));
				assert(0 == fgetpos(fps_qual[i], &pos_qual));
				read_name = get_read_name(fps_csfasta[i], fps_qual[i]); // Get read name
				if(NULL == read_name) { // eof
					fclose(fps_csfasta[i]);
					fclose(fps_qual[i]);
					fps_csfasta[i] = fps_qual[i] = NULL;
				}
				else {
					// check if the read name is the min
					if(NULL == min_read_name || 0 == cmp_read_names(read_name, min_read_name)) {
						if(NULL == min_read_name) {
							min_read_name = strdup(read_name);
						}
						to_print[i] = 1;
					}
					else if(cmp_read_names(read_name, min_read_name) < 0) {
						free(min_read_name);
						min_read_name = strdup(read_name);
						// Re-initialize other fps
						for(j=0;j<i;j++) {
							to_print[j] = 0;
						}
						to_print[i] = 1;
					}
					else {
						to_print[i] = 0;
					}
					// Free
					free(read_name);
					// Seek back
					assert(0 == fsetpos(fps_csfasta[i], &pos_csfasta));
					assert(0 == fsetpos(fps_qual[i], &pos_qual));
				}
			}
			else {
				to_print[i] = 0;
			}
		}
		// Print all with min read name
		more_fps_left = 0;
		for(i=0;i<number_of_ends;i++) {
			if(1 == to_print[i]) {
				more_fps_left++;
				print_fastq(fp_output, fps_csfasta[i], fps_qual[i]);
			}
		}
		// Update counts
		output_count++;
		// Open a new output file if necessary
		if(0 < num_reads_per_file &&
				num_reads_per_file <= output_count) {
			output_suffix_number++;
			fp_output = open_output_file(output_prefix, output_suffix_number, num_reads_per_file); 
			output_count=0;
		}
	}

	// Free
	free(output_prefix);
	for(i=0;i<number_of_ends;i++) {
		free(csfasta_filenames[i]);
		free(qual_filenames[i]);
	}
	free(csfasta_filenames);
	free(qual_filenames);
	free(fps_csfasta);
	free(fps_qual);

	return 0;
}

FILE *open_output_file(char *output_prefix, int32_t output_suffix_number, int32_t num_reads_per_file)
{
	char *FnName="open_output_file";
	char output_filename[1024]="\0";
	FILE *fp_out=NULL;

	// Create output file name
	if(0 < num_reads_per_file) {
		assert(0 < sprintf(output_filename, "%s.%d.fastq", output_prefix, output_suffix_number));
	}
	else {
		assert(0 < sprintf(output_filename, "%s.fastq", output_prefix));
	}

	// Open an output file
	if(!(fp_out = fopen(output_filename, "w"))) {
		PrintError(FnName,
				output_filename,
				"Could not open file for writing",
				Exit,
				OpenFileError);
	}

	return fp_out;
}

void print_fastq(FILE *output_fp, FILE *csfasta_fp, FILE *qual_fp)
{
	char *FnName="print_fastq";
	char csfasta_name[SEQUENCE_NAME_LENGTH]="\0";
	char read[SEQUENCE_LENGTH]="\0";
	char qual_name[SEQUENCE_NAME_LENGTH]="\0";
	char qual[SEQUENCE_LENGTH]="\0";

	// Read in
	if(NULL == fgets(csfasta_name, SEQUENCE_NAME_LENGTH, csfasta_fp) ||
			NULL == fgets(reads, SEQUENCE_LENGTH, csfasta_fp)) {
		PrintError(FnName,
				"csfasta",
				"Could not read in from file",
				Exit,
				ReadFileError);
	}
	if(NULL == fgets(qual_name, SEQUENCE_NAME_LENGTH, qual_fp) ||
			NULL == fgets(qual, SEQUENCE_LENGTH, qual_fp)) {
		PrintError(FnName,
				"qual",
				"Could not read in from file",
				Exit,
				ReadFileError);
	}

	if(0 != strcmp(csfasta_name, qual_name)) {
		PrintError(FnName,
				"csfasta_name != qual_name",
				"Read names did not match",
				Exit,
				OutOfRange);
	}

	// Convert SOLiD qualities
	// TODO
	
	// Print out
	// TODO
}

/*
   char *get_read_name(FILE*, FILE*);
   int32_t cmp_read_names(char*, char*);
   */
