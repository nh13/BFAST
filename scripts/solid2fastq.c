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

typedef struct {
	int32_t to_print; // whether this entry should be printed (must be populated)
	int32_t is_pop; // whether this entry has been populated
	char name[SEQUENCE_NAME_LENGTH];
	char read[SEQUENCE_LENGTH];
	char qual[SEQUENCE_LENGTH];
} fastq_t;

FILE *open_output_file(char*, int32_t, int32_t);
void fastq_print(fastq_t*, FILE*);
void fastq_read(fastq_t*, FILE*, FILE*);
int32_t cmp_read_names(char*, char*);
void read_name_trim(char*);
char *strtok_mod(char*, char*, int32_t*);

int main(int argc, char *argv[])
{
	char *output_prefix=NULL;
	int32_t num_reads_per_file=-1;
	int32_t no_output=0;
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
	fastq_t *reads=NULL;
	int32_t more_fps_left=1;
	int64_t output_count=0;
	int64_t output_count_total=0;
	char *min_read_name=NULL;

	// Get Parameters
	while((c = getopt(argc, argv, "n:o:c")) >= 0) {
		switch(c) {
			case 'c':
				no_output=1; break;
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
		fprintf(stderr, "\t-c\t\tproduce no output.\n");
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

	reads = malloc(sizeof(fastq_t)*number_of_ends);
	if(NULL == reads) {
		PrintError("solid2fastq",
				"reads",
				"Could not allocate memory",
				Exit,
				MallocMemory);
	}

	for(i=0;i<number_of_ends;i++) {
		reads[i].to_print = 0;
		reads[i].is_pop = 0;
	}

	output_suffix_number = 1;
	more_fps_left=number_of_ends;
	output_count = output_count_total = 0;
	// Open output file
	if(0 == no_output) {
	fp_output = open_output_file(output_prefix, output_suffix_number, num_reads_per_file); 
	}
	fprintf(stderr, "Outputting, currently on:\n0");
	while(0 < more_fps_left) { // while an input file is still open

		if(0 == (output_count_total % 100000)) {
			fprintf(stderr, "\r%lld",
					(long long int)output_count_total);
		}
		/*
		   fprintf(stderr, "more_fps_left=%d\n",
		   more_fps_left);
		   */
		// Get all with min read name
		for(i=0;i<number_of_ends;i++) {
			// populate read if necessary
			if(0 == reads[i].is_pop &&
					NULL != fps_csfasta[i] &&
					NULL != fps_qual[i]) {
				fastq_read(&reads[i], fps_csfasta[i], fps_qual[i]); // Get read name
				if(0 == reads[i].is_pop) { // was not populated
					//fprintf(stderr, "EOF\n");
					fclose(fps_csfasta[i]);
					fclose(fps_qual[i]);
					fps_csfasta[i] = fps_qual[i] = NULL;
					reads[i].to_print = 0;
				}
			}
			else {
				reads[i].to_print = 0;
			}
			// check if the read name is the min
			if(1 == reads[i].is_pop) {
				/*
				fprintf(stdout, "i=%d\tmin_read_name=%s\treads[i].name=%s\n",
						i,
						min_read_name,
						reads[i].name);
						*/
				if(NULL == min_read_name || 
						0 == cmp_read_names(reads[i].name, min_read_name)) {
					if(NULL == min_read_name) {
						min_read_name = strdup(reads[i].name);
					}
					reads[i].to_print = 1;
				}
				else if(cmp_read_names(reads[i].name, min_read_name) < 0) {
					free(min_read_name);
					min_read_name = strdup(reads[i].name);
					// Re-initialize other fps
					for(j=0;j<i;j++) {
						reads[j].to_print = 0;
					}
					reads[i].to_print = 1;
				}
				else {
					reads[i].to_print = 0;
				}
			}
			else {
				assert(0 == reads[i].to_print);
			}
		}
		/*
		fprintf(stdout, "min_read_name was %s\n",
				min_read_name);
				*/
		free(min_read_name);
		min_read_name=NULL;

		// Print all with min read name
		more_fps_left = 0;
		num_ends_printed = 0;
		for(i=0;i<number_of_ends;i++) {
			if(1 == reads[i].to_print) {
				more_fps_left++;
				num_ends_printed++;
				if(0 == no_output) {
					fastq_print(&reads[i], fp_output);
				}
				reads[i].is_pop = reads[i].to_print = 0;
			}
		}
		if(0 < num_ends_printed) {
			end_counts[num_ends_printed-1]++;
			// Update counts
			output_count++;
			output_count_total++;
		}
		// Open a new output file if necessary
		if(0 < num_reads_per_file &&
				num_reads_per_file <= output_count) {
			output_suffix_number++;
			if(0 == no_output) {
			fclose(fp_output);
			fp_output = open_output_file(output_prefix, output_suffix_number, num_reads_per_file); 
			}
			output_count=0;
		}
	}
	if(0 < output_count && 0 == no_output) {
		fclose(fp_output);
	}
	fprintf(stderr, "\r%lld\n",
			(long long int)output_count_total);

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
	free(reads);

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

void fastq_print(fastq_t *read, FILE *output_fp)
{
	//char *FnName="fastq_print";
	int32_t i;

	// Print out
	fprintf(output_fp, "@%s\n%s\n+\n",
			read->name,
			read->read);
	for(i=0;i<strlen(read->read)-1;i++) {
		fprintf(output_fp, "%c", read->qual[i]);
	}
	fprintf(output_fp, "\n");
}

void fastq_read(fastq_t *read, FILE *fp_csfasta, FILE *fp_qual) 
{
	char *FnName="fastq_read";
	char qual_name[SEQUENCE_NAME_LENGTH]="\0";
	int32_t qual[SEQUENCE_LENGTH];
	int32_t i;

	assert(0 == read->is_pop);
	assert(0 == read->to_print);

	if(NULL == fp_csfasta ||
			NULL == fp_qual ||
			0 != feof(fp_csfasta) || 
			0 != feof(fp_qual)) {
		return;
	}

	// Read in
	if(fscanf(fp_csfasta, "%s", read->name) < 0 ||
			fscanf(fp_csfasta, "%s", read->read) < 0) {
		return;
	}
	if(fscanf(fp_qual, "%s", qual_name) < 0) {
		return;
	}
	StringTrimWhiteSpace(read->name);
	StringTrimWhiteSpace(read->read);
	StringTrimWhiteSpace(qual_name);
	// Read in qual
	for(i=0;i<strlen(read->read)-1;i++) {
		if(fscanf(fp_qual, "%d", &qual[i]) < 0) {
			PrintError(FnName,
					"qual[i]",
					"Could not read in from file",
					Exit,
					ReadFileError);
		}
	}

	// Check that the read name and csfasta name match
	if(0 != strcmp(read->name, qual_name)) {
		fprintf(stderr, "csfasta_name=[%s]\nqual_name=[%s]\n", read->name, qual_name);
		PrintError(FnName,
				"read->name != qual_name",
				"Read names did not match",
				Exit,
				OutOfRange);
	}
	// Remove leading '@' from the read name
	for(i=1;i<strlen(read->name);i++) {
		read->name[i-1] = read->name[i];
	}
	read->name[i]='\0';

	// Convert SOLiD qualities
	for(i=0;i<strlen(read->read)-1;i++) {
		/*
		   fprintf(stderr, "%c -> %c\n%d -> %d\n",
		   qual[i],
		   (qual[i] <= 93 ? qual[i] : 93) + 33,
		   qual[i],
		   (qual[i] <= 93 ? qual[i] : 93) + 33);
		   exit(1);
		   */
		qual[i] = 33 + (qual[i] <= 0 ? 1 : (qual[i] <= 93 ? qual[i] : 93));
		read->qual[i] = (char)qual[i];
	}

	// Trim last _R3 or _F3 or _whatever
	read_name_trim(read->name);

	read->is_pop = 1;
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
			free(name_one_cur);
			free(name_two_cur);
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
	char *r=strdup(str + (*index));

	assert(NULL != str);

	while((*index) < strlen(str)) {
		for(i=0;i<strlen(delim);i++) {
			if(delim[i] == str[(*index)]) {
				r[(*index)-prev_index]='\0';
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
