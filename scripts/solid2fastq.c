#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>
#include <string.h>
#include "../blib/BLibDefinitions.h"
#include "../blib/BError.h"
#include "../blib/BLib.h"

FILE *open_output_file(char*, int32_t, int32_t);
void print_fastq(FILE*, FILE*, FILE*);
char *get_read_name(FILE*, FILE*);
int32_t cmp_read_names(char*, char*);
void read_name_trim(char*);
char *strtok_mod(char*, char*, int32_t*);

int main(int argc, char *argv[])
{
	char *output_prefix=NULL;
	int32_t num_reads_per_file=-1;
	int32_t number_of_ends;
	int32_t num_ends_printed = 0;
	int64_t *end_counts=NULL;
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
		return 0;
	}

	// Copy over the filenames
	assert(0 == (argc - optind) % 2);
	number_of_ends = (argc - optind) / 2;
	// Allocate memory
	csfasta_filenames = malloc(sizeof(char*)*number_of_ends);
	if(NULL == csfasta_filenames) {
		PrintError("solid2fastq",
				"csfasta_filenames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	qual_filenames = malloc(sizeof(char*)*number_of_ends);
	if(NULL == qual_filenames) {
		PrintError("solid2fastq",
				"qual_filenames",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	end_counts = malloc(sizeof(int64_t)*number_of_ends);
	if(NULL == end_counts) {
		PrintError("solid2fastq",
				"end_counts",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}
	for(i=0;i<number_of_ends;i++) {
		csfasta_filenames[i] = strdup(argv[optind+i]);
		qual_filenames[i] = strdup(argv[optind+i+number_of_ends]);
		end_counts[i] = 0;
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
	int64_t output_count=0;

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
	fprintf(stderr, "Outputting, currently on:\n0");
	while(0 < more_fps_left) { // while an input file is still open

		if(0 == (output_count % 100000)) {
			fprintf(stderr, "\r%lld",
					(long long int)output_count);
		}
		/*
		   fprintf(stderr, "more_fps_left=%d\n",
		   more_fps_left);
		   */
		// Get all with min read name
		char *min_read_name=NULL;
		for(i=0;i<number_of_ends;i++) {
			to_print[i] = 0;
			/*
			   fprintf(stderr, "%d\t%d\t%d\n",
			   i,
			   NULL == fps_csfasta[i],
			   NULL == fps_qual[i]);
			   */
			if(NULL != fps_csfasta[i] && NULL != fps_qual[i]) {
				/*
				   fprintf(stderr, "HERE 1\n");
				   */
				fpos_t pos_csfasta, pos_qual;
				char *read_name=NULL;
				assert(0 == fgetpos(fps_csfasta[i], &pos_csfasta));
				assert(0 == fgetpos(fps_qual[i], &pos_qual));
				read_name = get_read_name(fps_csfasta[i], fps_qual[i]); // Get read name
				if(NULL == read_name) { // eof
					/*
					   fprintf(stderr, "EOF\n");
					   */
					fclose(fps_csfasta[i]);
					fclose(fps_qual[i]);
					fps_csfasta[i] = fps_qual[i] = NULL;
					to_print[i] = 0;
				}
				else {
					/*
					   fprintf(stderr, "%d\tread_name=%s\tmin_read_name=%s\n",
					   i,
					   read_name,
					   min_read_name);
					   */
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
		num_ends_printed = 0;
		for(i=0;i<number_of_ends;i++) {
			if(1 == to_print[i]) {
				more_fps_left++;
				print_fastq(fp_output, fps_csfasta[i], fps_qual[i]);
				num_ends_printed++;
			}
		}
		if(0 < num_ends_printed) {
			end_counts[num_ends_printed-1]++;
			// Update counts
			output_count++;
		}
		// Open a new output file if necessary
		if(0 < num_reads_per_file &&
				num_reads_per_file <= output_count) {
			output_suffix_number++;
			fp_output = open_output_file(output_prefix, output_suffix_number, num_reads_per_file); 
			output_count=0;
		}
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)output_count);

	fprintf(stderr, "Found\n%16s\t%16s\n", 
			"number_of_ends",
			"number_of_reads");

	for(i=0;i<number_of_ends;i++) {
		fprintf(stderr, "%16d\t%16lld\n",
				i+1,
				(long long int)end_counts[i]);
	}

	// Free
	free(output_prefix);
	for(i=0;i<number_of_ends;i++) {
		free(csfasta_filenames[i]);
		free(qual_filenames[i]);
	}
	free(csfasta_filenames);
	free(qual_filenames);
	free(end_counts);
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
	int32_t qual[SEQUENCE_LENGTH];
	int32_t i;

	if(0 != feof(csfasta_fp) || 0 != feof(qual_fp)) {
		return;
	}

	// Read in
	if(fscanf(csfasta_fp, "%s", csfasta_name) < 0 ||
			fscanf(csfasta_fp, "%s", read) < 0) {
		PrintError(FnName,
				"csfasta",
				"Could not read in from file",
				Exit,
				ReadFileError);
	}
	if(fscanf(qual_fp, "%s", qual_name) < 0) {
		PrintError(FnName,
				"qual_name",
				"Could not read in from file",
				Exit,
				ReadFileError);
	}
	StringTrimWhiteSpace(csfasta_name);
	StringTrimWhiteSpace(read);
	StringTrimWhiteSpace(qual_name);
	// Read in qual
	for(i=0;i<strlen(read)-1;i++) {
		if(fscanf(qual_fp, "%d", &qual[i]) < 0) {
			PrintError(FnName,
					"qual[i]",
					"Could not read in from file",
					Exit,
					ReadFileError);
		}
	}

	// Check that the read name and csfasta name match
	if(0 != strcmp(csfasta_name, qual_name)) {
		   fprintf(stderr, "csfasta_name=[%s]\nqual_name=[%s]\n", csfasta_name, qual_name);
		PrintError(FnName,
				"csfasta_name != qual_name",
				"Read names did not match",
				Exit,
				OutOfRange);
	}

	// Convert SOLiD qualities
	for(i=0;i<strlen(read)-1;i++) {
		/*
		   fprintf(stderr, "%c -> %c\n%d -> %d\n",
		   qual[i],
		   (qual[i] <= 93 ? qual[i] : 93) + 33,
		   qual[i],
		   (qual[i] <= 93 ? qual[i] : 93) + 33);
		   exit(1);
		   */
		qual[i] = (qual[i] <= 93 ? qual[i] : 93) + 33;
	}

	// Print out
	assert('>' == csfasta_name[0]);
	csfasta_name[0]='@'; // assumes '>' is the first character
	fprintf(output_fp, "%s\n%s\n+\n",
			csfasta_name,
			read);
	for(i=0;i<strlen(read)-1;i++) {
		fprintf(output_fp, "%c", (char)qual[i]);
	}
	fprintf(output_fp, "\n");
}

char *get_read_name(FILE *fp_csfasta, FILE *fp_qual)
{
	char *FnName="get_read_name";
	fpos_t pos_csfasta, pos_qual;
	char name_csfasta[SEQUENCE_NAME_LENGTH]="\0";
	char name_qual[SEQUENCE_NAME_LENGTH]="\0";

	/*
	   fprintf(stderr, "In %s\n", FnName);
	   */

	// Get position
	if(0 != fgetpos(fp_csfasta, &pos_csfasta) ||
			0 != fgetpos(fp_qual, &pos_qual)) {
		return NULL;
	}

	// Read
	if(EOF == fscanf(fp_csfasta, "%s", name_csfasta) ||
			EOF == fscanf(fp_qual, "%s", name_qual)) {
		return NULL;
	}

	// Reset position
	if(0 != fsetpos(fp_csfasta, &pos_csfasta) ||
			0 != fsetpos(fp_qual, &pos_qual)) {
		return NULL;
	}

	/*
	   fprintf(stderr, "Checking %s with %s\n",
	   name_csfasta,
	   name_qual);
	   */
	// Trim last _R3 or _F3 or _whatever
	read_name_trim(name_csfasta);
	read_name_trim(name_qual);

	// Check the names match
	if(0 != cmp_read_names(name_csfasta, name_qual)) {
		/*
		   fprintf(stderr, "\n%s\n%s\n",
		   name_csfasta,
		   name_qual);
		   */
		PrintError(FnName,
				NULL,
				"The read and qual names did not match",
				Exit,
				OutOfRange);
	}

	// This could be dangerous
	return strdup(name_csfasta);
}

int32_t cmp_read_names(char *name_one, char *name_two)
{
	char *name_one_cur = NULL;
	char *name_two_cur = NULL;
	int32_t name_one_num_state = 0, name_two_num_state = 0;
	int32_t return_value = 0;
	int32_t name_one_index = 0, name_two_index = 0;
	/*
	   fprintf(stderr, "comparing %s with %s\n",
	   name_one, name_two);
	   */

	name_one_cur = strtok_mod(name_one, "_", &name_one_index);
	name_two_cur = strtok_mod(name_two, "_", &name_two_index);

	while(NULL != name_one_cur && NULL != name_two_cur) {
		/*
		   fprintf(stderr, "name_one_cur=%s\nname_two_cur=%s\n",
		   name_one_cur,
		   name_two_cur);
		   */
		// assumes positive
		name_one_num_state = ( name_one_cur[0] < '0' || name_one_cur[0] > '9') ? 0 : 1;
		name_two_num_state = ( name_two_cur[0] < '0' || name_two_cur[0] > '9') ? 0 : 1;
		if(1 == name_one_num_state && 1 == name_two_num_state) {
			name_one_num_state = atoi(name_one_cur);
			name_two_num_state = atoi(name_two_cur);
			return_value = (name_one_num_state < name_two_num_state) ? -1 : ((name_one_num_state == name_two_num_state) ? 0 : 1); 
		}
		else {
			return_value = strcmp(name_one_cur, name_two_cur);
		}
		/*
		   fprintf(stderr, "return_value=%d\n", return_value);
		   */
		if(0 != return_value) {
			return return_value;
		}

		// Get next tokens
		free(name_one_cur);
		free(name_two_cur);
		name_one_cur = strtok_mod(name_one, "_", &name_one_index);
		name_two_cur = strtok_mod(name_two, "_", &name_two_index);
	}

	free(name_one_cur);
	free(name_two_cur);

	if(NULL != name_one_cur && NULL == name_two_cur) {
		return 1;
	}
	else if(NULL == name_one_cur && NULL != name_two_cur) {
		return -1;
	}
	else {
		return 0;
	}
}

void read_name_trim(char *name)
{
	int32_t i;

	// Trim last _R3 or _F3 or _whatever

	if(NULL == name) {
		return;
	}

	for(i=strlen(name)-1;0<i;i--) {
		switch(name[i]) {
			case '_':
				// truncate
				name[i]='\0';
				//fprintf(stderr, "return i=%d\n%s\n", i, name);
				return;
				break;
			default:
				break;
		}
	}

	assert('_' != name[0]);
}

char *strtok_mod(char *str, char *delim, int32_t *index)
{
	int32_t i, prev_index=(*index);
	char *r=strdup(str);

	assert(NULL != str);

	while((*index) < strlen(str)) {
		for(i=0;i<strlen(delim);i++) {
			if(delim[i] == str[(*index)]) {
				r[(*index)]='\0';
				(*index)++;
				return r;
			}
		}
		(*index)++;
	}

	if(prev_index == (*index)) {
		free(r);
		return NULL;
	}
	else {
		return r;
	}
}
