#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <config.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "../bfast/BLibDefinitions.h"
#include "../bfast/BError.h"
#include "../bfast/BLib.h"
#include "../bfast/aflib.h"

#define Name "solid2fastq"

#define MALLOC_AND_TEST(VAR, SIZE) \
	VAR = malloc(SIZE); \
	if(NULL == VAR) { \
		PrintError(Name, #VAR, "Could not allocate memory", Exit, MallocMemory); \
	}

#define CHECK_RM_EMPTY_FILE(OUTPUT_COUNTS, FN, FH, CREATE_FN, OUT_COMP)	\
	if (0 == OUTPUT_COUNTS) {					\
		CREATE_FN;						\
		switch(OUT_COMP) {					\
 			case AFILE_GZ_COMPRESSION:			\
				strcat(FN, ".gz"); break;		\
			case AFILE_BZ2_COMPRESSION:			\
				strcat(FN, ".bz2"); break;		\
			 default:					\
				break;					\
		}							\
		if ((FH = fopen(FN, "r")) != NULL) {			\
			if(remove(FN)) PrintError(Name, #FN, "Cannot remove file", Exit, DeleteFileError); \
			fclose(FH);					\
		}							\
	}


enum fastq_read_type { read1, read2, single, combined, undefined };

typedef struct {
	int32_t to_print; // whether this entry should be printed (must be populated)
	int32_t is_pop; // whether this entry has been populated
	char name[SEQUENCE_NAME_LENGTH];
	char read[SEQUENCE_LENGTH];
	char qual[SEQUENCE_LENGTH];
} fastq_t;

void open_output_file(char*, enum fastq_read_type, int32_t, int32_t, int32_t, AFILE **, int32_t, int32_t);
void open_output_files(char*, int32_t**, int32_t, int32_t, AFILE **, int32_t, int32_t, int32_t);
void fastq_print(fastq_t*, AFILE**, int32_t, char *, int32_t, fastq_t *, int32_t, int32_t, int32_t, int64_t **);
void dump_read(AFILE *afp_output, fastq_t *read, int64_t *output_count);
int32_t cmp_read_names(char*, char*, int32_t);
void read_name_trim(char*);
char *strtok_mod(char*, char*, int32_t*);
void add_read_name_prefix(fastq_t *, char *, int32_t);
void to_bwa(fastq_t *, int32_t);
void close_fds(AFILE **, int32_t, char *, int32_t, int32_t);
//int is_empty(char *path);
#ifdef HAVE_FSEEKO
void fastq_read(fastq_t*, AFILE*, AFILE*, off_t, off_t, int32_t, int32_t);
int32_t read_line(AFILE *afp, off_t end_pos, char *line);
#else
void fastq_read(fastq_t*, AFILE*, AFILE*, int32_t, int32_t);
int32_t read_line(AFILE *afp, char *line);
#endif
int print_usage ()
{
	fprintf(stderr, "solid2fastq %s\n", PACKAGE_VERSION);
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: solid2fastq [options] <list of .csfasta files> <list of .qual files>\n");
	fprintf(stderr, "\t-c\t\tproduce no output.\n");
	fprintf(stderr, "\t-n\tINT\tnumber of reads per file.\n");
	fprintf(stderr, "\t-o\tSTRING\toutput prefix.\n");
	fprintf(stderr, "\t-p\tSTRING\tread name prefix (added to front of all read names).\n");
	fprintf(stderr, "\t-j\t\tinput files are bzip2 compressed.\n");
	fprintf(stderr, "\t-z\t\tinput files are gzip compressed.\n");
	fprintf(stderr, "\t-J\t\toutput files are bzip2 compressed.\n");
	fprintf(stderr, "\t-Z\t\toutput files are gzip compressed.\n");
	fprintf(stderr, "\t-t\tINT\ttrim INT bases from the 3' end of the reads.\n");
	fprintf(stderr, "\t-b\t\tEnable bwa output (for 'bwa aln', not for 'bfast bwaaln').\n");
	fprintf(stderr, "\t-w\t\tCreate a single file to dump reads with only one end.\n");
#ifdef HAVE_FSEEKO
        fprintf(stderr, "\t-s\tINT,INT[,INT,INT]\tstart reading at the given byte location in each file\n");
        fprintf(stderr, "\t-e\tINT,INT[,INT,INT]\tstop reading at the given byte location in each file\n");
	fprintf(stderr, "\t-N\tINT\tFile number (appended to end of output prefix)\n");
	fprintf(stderr, "\t-h\t\tprint this help message.\n");
        fprintf(stderr, "\n\n-s and -e are used to split files in parallel on a cluster.");
#else
	fprintf(stderr, "\t-h\t\tprint this help message.\n");
#endif
        fprintf(stderr, "\n");
	fprintf(stderr, "\n send bugs to %s\n", PACKAGE_BUGREPORT);
	return 1;
}

#ifdef HAVE_FSEEKO
int parse_locations(char *options, int64_t *csfasta_loc, int64_t *qual_loc) {
	char *str_tmp1;
	char *str_tmp2;
	char *str_tmp3;
	char *str_tmp4;

	str_tmp1 = strtok(options, ",");
	str_tmp2 = strtok(NULL, ",");
	str_tmp3 = strtok(NULL, ",");
	str_tmp4 = strtok(NULL, ",");

	if (str_tmp1 == NULL || str_tmp2 == NULL || (str_tmp3 != NULL && str_tmp4 == NULL)) {
		PrintError(Name, "csfasta_loc, qual_loc", "not enough byte location arguments", Exit, OutOfRange);
	}

	if (str_tmp3 == NULL) {
		csfasta_loc[0] = atoll(str_tmp1);
		qual_loc[0] = atoll(str_tmp2);
		return 1; // total number of csfasta/qual pairs found
	}

	csfasta_loc[0] = atoll(str_tmp1);
	csfasta_loc[1] = atoll(str_tmp2);
	qual_loc[0] = atoll(str_tmp3);
	qual_loc[1] = atoll(str_tmp4);
	
	return 2;  // Number of csfasta/qual pairs found
}
#endif

int main(int argc, char *argv[])
{
	char *output_prefix=NULL;
	char *read_name_prefix=NULL;
	int32_t rnp_len=0; // length of read_name_prefix
	int32_t num_reads_per_file=-1;
	int32_t no_output=0;
	int32_t bwa_output=0;
	int32_t single_output=0;
	int32_t number_of_ends;
	int32_t num_ends_printed = 0;
	int64_t *end_counts=NULL;
	char **csfasta_filenames=NULL;
	char **qual_filenames=NULL;
	AFILE **afps_csfasta=NULL;
	AFILE **afps_qual=NULL;
	AFILE *afp_output[3]; // Necessary for BWA, single output (it uses three file descriptors (read1, read2, single-end reads)
	enum fastq_read_type output_read_type[3];
	int32_t in_comp=AFILE_NO_COMPRESSION;
	int32_t out_comp=AFILE_NO_COMPRESSION;
	int32_t *output_suffix_number;
	int c;
	int32_t i, j;
	int32_t trim_end = 0;
	fastq_t *reads=NULL;
	int32_t more_afps_left=1;
	int64_t *output_counts;
	int64_t output_count_total=0;
	char *min_read_name=NULL;
	int32_t prev=0;
	int32_t num_output_files = 1;
#ifdef HAVE_FSEEKO
	off_t csfasta_start_pos[2];
	off_t csfasta_end_pos[2];
	off_t qual_start_pos[2];
	off_t qual_end_pos[2];
	int start_pos_counts = 0;
	int end_pos_counts = 0;
	int file_number = -1;
#endif

	// Get Parameters
	while((c = getopt(argc, argv, "n:o:t:p:chjzJZbws:e:N:")) >= 0) {
		switch(c) {
			case 'b':
				bwa_output=1; break;
			case 'c':
				no_output=1; break;
			case 'h':
				return print_usage(); break;
			case 'j':
				in_comp=AFILE_BZ2_COMPRESSION; break;
			case 'n':
				num_reads_per_file=atoi(optarg); break;
			case 'o':
				output_prefix=strdup(optarg); break;
			case 'p':
				MALLOC_AND_TEST(read_name_prefix, strlen(optarg) + 3);
				strcpy(read_name_prefix, optarg); 
				strcat(read_name_prefix, ":");
				rnp_len = strlen(read_name_prefix);
				break;
			case 't':
				trim_end=atoi(optarg); break;
			case 'w':
				single_output=1; break;
			case 'z':
				in_comp=AFILE_GZ_COMPRESSION; break;
			case 'J':
				out_comp=AFILE_BZ2_COMPRESSION; break;
			case 'Z':
				out_comp=AFILE_GZ_COMPRESSION; break;
#ifdef HAVE_FSEEKO
			case 's':
				start_pos_counts = parse_locations(optarg, csfasta_start_pos, qual_start_pos); break;
			case 'e':
				end_pos_counts = parse_locations(optarg, csfasta_end_pos, qual_end_pos); break;
			case 'N':
				file_number = atoi(optarg); break;
#endif
			default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}

	// Print Usage
	if(argc == optind ||
			0 != ((argc - optind) % 2) ||
			(argc - optind) < 2) {
		return print_usage();
	}

	// Validate argument counts
	assert(0 == (argc - optind) % 2);
	number_of_ends = (argc - optind) / 2;

#ifdef HAVE_FSEEKO
	if (start_pos_counts > 0 || end_pos_counts > 0) {
		if (start_pos_counts != end_pos_counts ||
		    start_pos_counts != number_of_ends) {
			PrintError(Name, "start_pos_counts", 
				   "Number of start/end positions not consistent with number of files",
				   Exit, OutOfRange);
		}
	}
#endif

	// Single output mode is only useful for paired-end data
	if (1 == number_of_ends) {
		single_output = 0;
	}

	if (2 == number_of_ends && (1 == bwa_output || 1 == single_output)) {
		num_output_files = 3;
	} 
	else {
		num_output_files = 1;
	}

	// bwa expects a read_name_prefix
	if (1 == bwa_output && read_name_prefix == NULL && output_prefix != NULL) {
		char *prefix_start;

		// Remove path from beginning of output prefix
		if (NULL == (prefix_start = strrchr(output_prefix, '/'))) {
			prefix_start = output_prefix;
		}
		else {
			prefix_start += 1; // need to move past '/'
		}

		// Add a colon to the read_name_prefix
		MALLOC_AND_TEST(read_name_prefix, strlen(prefix_start) + 3);
		strcpy(read_name_prefix, prefix_start); 
		strcat(read_name_prefix, ":");
		rnp_len = strlen(read_name_prefix);
	}

#ifdef HAVE_FSEEKO
	// Add number to output_prefix, if necessary
	// Do this AFTER copying output_prefix to read_name_prefix (block above)

	if (file_number > 0 && NULL != output_prefix) {
		char *tmp_output_prefix = strdup(output_prefix);
		free(output_prefix);
		MALLOC_AND_TEST(output_prefix, strlen(tmp_output_prefix) + (int)(log10((double)file_number)) + 3);
		sprintf(output_prefix, "%s.%d", tmp_output_prefix, file_number);
	}
#endif

	// Copy over the filenames

	// Allocate memory
	MALLOC_AND_TEST(csfasta_filenames, sizeof(char*)*number_of_ends);
	MALLOC_AND_TEST(qual_filenames, sizeof(char*)*number_of_ends);
	MALLOC_AND_TEST(end_counts, sizeof(int64_t)*number_of_ends);

	for(i=0;i<number_of_ends;i++) {
		csfasta_filenames[i] = strdup(argv[optind+i]);
		qual_filenames[i] = strdup(argv[optind+i+number_of_ends]);
		end_counts[i] = 0;
	}

	// Allocate memory for input file pointers
	MALLOC_AND_TEST(afps_csfasta, sizeof(AFILE*)*number_of_ends);
	MALLOC_AND_TEST(afps_qual, sizeof(AFILE*)*number_of_ends);

	// Open input files
	for(i=0;i<number_of_ends;i++) {
		if(!(afps_csfasta[i] = AFILE_afopen(csfasta_filenames[i], "rb", in_comp))) {
			PrintError(Name, csfasta_filenames[i], "Could not open file for reading", Exit, OpenFileError);
		}
		if(!(afps_qual[i] = AFILE_afopen(qual_filenames[i], "rb", in_comp))) {
			PrintError(Name, qual_filenames[i], "Could not open file for reading", Exit, OpenFileError);
		}
#ifdef HAVE_FSEEKO
		if (start_pos_counts > 0) {
			AFILE_afseek(afps_csfasta[i], csfasta_start_pos[i], SEEK_SET);
			AFILE_afseek(afps_qual[i], qual_start_pos[i], SEEK_SET);
		}
#endif
	}

	MALLOC_AND_TEST(reads, sizeof(fastq_t)*number_of_ends);

	for(i=0;i<number_of_ends;i++) {
		reads[i].to_print = 0;
		reads[i].is_pop = 0;
	}

	MALLOC_AND_TEST(output_suffix_number, sizeof(int32_t)*3);
	MALLOC_AND_TEST(output_counts, sizeof(int64_t)*3);

	for(i=0;i<num_output_files;i++) {
		output_suffix_number[i] = 1;
		output_counts[i] = 0;
	}

	more_afps_left=number_of_ends;
	output_count_total = 0;

	// Set up output file read types
	if (bwa_output) {
		output_read_type[0] = single;
		output_read_type[1] = read2;
		output_read_type[2] = read1;
	}
	else if (single_output) {
		output_read_type[0] = single;
		output_read_type[1] = read1;
		output_read_type[2] = read2;
	}
	else {
		output_read_type[0] = combined;
		output_read_type[1] = undefined;
		output_read_type[2] = undefined;
	}


	// Open output file

	afp_output[0] = afp_output[1] = afp_output[2] = NULL;

	if(NULL == output_prefix) {
		if(!(afp_output[0] = AFILE_afdopen(fileno(stdout), "wb", out_comp))) {
			PrintError(Name, "stdout", "Could not open for writing", Exit, WriteFileError);
		}
	}
	else if(0 == no_output) {
		open_output_files(output_prefix, &output_suffix_number, num_reads_per_file, out_comp, afp_output, bwa_output, number_of_ends, single_output);
	}

	fprintf(stderr, "Outputting, currently on:\n0");
	while(0 < more_afps_left) { // while an input file is still open

		if(0 == (output_count_total % 100000)) {
			fprintf(stderr, "\r%lld", (long long int)output_count_total);
		}
		/*
		   fprintf(stderr, "more_afps_left=%d\n", more_afps_left);
		   */
		// Get reads (one at a time) and set the reads data structure (0: read1, 1: read2)
		for(i=0;i<number_of_ends;i++) {
			// populate read if necessary
			if(0 == reads[i].is_pop &&
					NULL != afps_csfasta[i] &&
					NULL != afps_qual[i]) {
#ifdef HAVE_FSEEKO
				fastq_read(&reads[i], afps_csfasta[i], afps_qual[i], csfasta_end_pos[i], qual_end_pos[i], trim_end, bwa_output); // Get read name
#else
				fastq_read(&reads[i], afps_csfasta[i], afps_qual[i], trim_end, bwa_output); // Get read name
#endif
				if(0 == reads[i].is_pop) { // was not populated
					//fprintf(stderr, "EOF\n");
					AFILE_afclose(afps_csfasta[i]);
					AFILE_afclose(afps_qual[i]);
					afps_csfasta[i] = afps_qual[i] = NULL;
					reads[i].to_print = 0;
				}
				else {
					// add read name prefix
					if (NULL != read_name_prefix) {
						add_read_name_prefix(&reads[i], read_name_prefix, rnp_len);
					}

					// BWA output enabled, transform data (read/qual) to bwa expected format
					if (1 == bwa_output) { 
						to_bwa(&reads[i], i);
					}
				}

			}
			else {
				reads[i].to_print = 0;
			}

			if(1 == reads[i].is_pop) {
				/* fprintf(stdout, "i=%d\tmin_read_name=%s\treads[i].name=%s\n", i, min_read_name, reads[i].name); */
				if(NULL == min_read_name ||
				   0 == cmp_read_names(reads[i].name, min_read_name, rnp_len)) {
					if(NULL == min_read_name) {
						min_read_name = strdup(reads[i].name);
					}
					reads[i].to_print = 1;
				}
				else if(cmp_read_names(reads[i].name, min_read_name, rnp_len) < 0) {
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
		} // end for
		/*
		   fprintf(stdout, "min_read_name was %s\n", min_read_name);
		   */
		free(min_read_name);
		min_read_name=NULL;

		// Print all with min read name
		more_afps_left = 0;
		num_ends_printed = 0;
		for(i=0;i<number_of_ends;i++) {
			// In bwa mode, I need to know if the first end was available or not to drop the 
			// read in single or in its proper end file (read1 / read 2)
			if (i==0 && number_of_ends == 2) { prev = reads[i].to_print; }
			if(1 == reads[i].to_print) {
				more_afps_left++;
				num_ends_printed++;
				if(0 == no_output) {
					fastq_print(&reads[i], afp_output, bwa_output, output_prefix, number_of_ends, reads, single_output, i, prev, &output_counts);
				}
				reads[i].is_pop = reads[i].to_print = 0;
			}
		}
		if(0 < num_ends_printed) {
			end_counts[num_ends_printed-1]++;
			// Update counts
			output_count_total++;

			// Open a new output file if necessary
			if(0 < num_reads_per_file) {
				for(i=0;i<num_output_files;i++) {
					if (0 == i && 1 == single_output) {
						// in single_output node, all single reads are dumped into
						// one file, so that file is not rotated.
						continue;
					}
					if (num_reads_per_file <= output_counts[i]) {
						if (NULL != output_prefix && 0 == no_output) {
							assert(NULL != afp_output[i]);
							AFILE_afclose(afp_output[i]);

							output_suffix_number[i]++;
							open_output_file(output_prefix, output_read_type[i], output_suffix_number[i], num_reads_per_file, out_comp, &afp_output[i], bwa_output, single_output);
							output_counts[i]=0;
						}
					}
				}
			}
		}
	} // end while

	//	if(0 < output_count && 0 == no_output) {
	if(0 == no_output) {
		close_fds(afp_output, bwa_output, output_prefix, number_of_ends, single_output);
	}

	// Remove last fastq file when total input reads % num_reads_per_file == 0
	// We don't want an empty file
	if(0 == no_output && NULL != output_prefix && num_reads_per_file > 0) {
		if (0 == bwa_output && 0 == single_output) {
			char empty_fn[4096]="\0";
			FILE *f;
			CHECK_RM_EMPTY_FILE( output_counts[0], empty_fn, f, 
					     sprintf(empty_fn, "%s.%d.fastq", output_prefix, output_suffix_number[0]),
					     out_comp);
		}

		if(1 == bwa_output) {
			char empty_fn[4096]="\0";
			FILE *f;
			CHECK_RM_EMPTY_FILE( output_counts[2], empty_fn, f, 
					     sprintf(empty_fn, "%s.read1.%d.fastq", output_prefix, output_suffix_number[2]),
					     out_comp );
			CHECK_RM_EMPTY_FILE( output_counts[1], empty_fn, f, 
					     sprintf(empty_fn, "%s.read2.%d.fastq", output_prefix, output_suffix_number[1]),
					     out_comp );
			CHECK_RM_EMPTY_FILE( output_counts[0], empty_fn, f, 
					     sprintf(empty_fn, "%s.single.%d.fastq", output_prefix, output_suffix_number[0]),
					     out_comp );

		}

		if(1 == single_output) {
			char empty_fn[4096]="\0";
			FILE *f;
			CHECK_RM_EMPTY_FILE( output_counts[1], empty_fn, f, 
					     sprintf(empty_fn, "%s.rd1.%d.fastq", output_prefix, output_suffix_number[1]),
					     out_comp );
			CHECK_RM_EMPTY_FILE( output_counts[2], empty_fn, f, 
					     sprintf(empty_fn, "%s.rd2.%d.fastq", output_prefix, output_suffix_number[2]),
					     out_comp );
			CHECK_RM_EMPTY_FILE( output_counts[0], empty_fn, f, 
					     sprintf(empty_fn, "%s.single.%d.fastq", output_prefix, output_suffix_number[0]),
					     out_comp );
		}
	}

	fprintf(stderr, "\r%lld\n", (long long int)output_count_total);
	fprintf(stderr, "Found\n%16s\t%16s\n", "number_of_ends", "number_of_reads");
	for(i=0;i<number_of_ends;i++) {
		fprintf(stderr, "%16d\t%16lld\n", i+1, (long long int)end_counts[i]);
	}

	// Free
	if(NULL != output_prefix) free(output_prefix);
	for(i=0;i<number_of_ends;i++) {
		free(csfasta_filenames[i]);
		free(qual_filenames[i]);
	}
	free(csfasta_filenames);
	free(qual_filenames);
	free(end_counts);
	free(afps_csfasta);
	free(afps_qual);
	free(reads);
	free(output_counts);
	free(output_suffix_number);

	return 0;
}

void open_output_file(char *output_prefix, enum fastq_read_type file_read_type, int32_t output_suffix_number, int32_t num_reads_per_file, 
		      int32_t out_comp, AFILE **fps, int32_t bwa_output, int32_t single_output) {

	char *FnName="open_output_file";
	char output_filename[4096]="\0";

	assert(file_read_type != undefined);

	if (file_read_type == combined) {
		// Create output file name
		if(0 < num_reads_per_file) {
			assert(0 < sprintf(output_filename, "%s.%d.fastq", output_prefix, output_suffix_number));
		}
		else {
			assert(0 < sprintf(output_filename, "%s.fastq", output_prefix));
		}
	}
	else {
		char *read_type = "";

		switch (file_read_type) {
		case read1: 
			if (1 == bwa_output) read_type = "read1";
			else if (1 == single_output) read_type = "rd1";
			break;
		case read2:
			if (1 == bwa_output) read_type = "read2";
			else if (1 == single_output) read_type = "rd2";
			break;
		case single:
			read_type = "single";
			break;
		default:
			break;
		}

		// For single mode, when splitting files by number,
		// only read1 and read2 have a suffix number.
		// All single reads are dumped into the same file (why?)

		if(0 < num_reads_per_file && !(1 == single_output && single == file_read_type)) {
			assert(0 < sprintf(output_filename, "%s.%s.%d.fastq", output_prefix, read_type, output_suffix_number));
		}
		else {
			assert(0 < sprintf(output_filename, "%s.%s.fastq", output_prefix, read_type));
		}
	}

	switch(out_comp) {
	case AFILE_GZ_COMPRESSION:
		strcat(output_filename, ".gz"); break;
	case AFILE_BZ2_COMPRESSION:
		strcat(output_filename, ".bz2"); break;
	default: 
		break;
	}

	// Open output files
	if(!(fps[0] = AFILE_afopen(output_filename, "wb", out_comp))) {
		PrintError(FnName, output_filename, "Could not open file for writing", Exit, OpenFileError);
	}


}

void open_output_files(char *output_prefix, int32_t **output_suffix_number, int32_t num_reads_per_file, 
		int32_t out_comp, AFILE **fps, int32_t bwa_output, int32_t number_of_ends, int32_t single_output)
{
	if (0 == bwa_output && 0 == single_output) { // Regular BF mode
		open_output_file(output_prefix, combined, (*output_suffix_number)[0], num_reads_per_file,
				 out_comp, &fps[0], bwa_output, single_output);
	}
	else if (1 == single_output) { // BF single mode (read1, read2, single)
		open_output_file(output_prefix, single, (*output_suffix_number)[0], num_reads_per_file,
				 out_comp, &fps[0], bwa_output, single_output);
		if (2 == number_of_ends) {
			open_output_file(output_prefix, read1, (*output_suffix_number)[1], num_reads_per_file,
					 out_comp, &fps[1], bwa_output, single_output);
			open_output_file(output_prefix, read2, (*output_suffix_number)[2], num_reads_per_file,
					 out_comp, &fps[2], bwa_output, single_output);
		}
	}
	else if (1 == bwa_output) { // BWA output detected
		open_output_file(output_prefix, single, (*output_suffix_number)[0], num_reads_per_file,
				 out_comp, &fps[0], bwa_output, single_output);
		if (2 == number_of_ends) {
			open_output_file(output_prefix, read2, (*output_suffix_number)[1], num_reads_per_file,
					 out_comp, &fps[1], bwa_output, single_output);
			open_output_file(output_prefix, read1, (*output_suffix_number)[2], num_reads_per_file,
					 out_comp, &fps[2], bwa_output, single_output);
		}
	}
}

void fastq_print(fastq_t *read, AFILE **fps, int32_t bwa_output, char *output_prefix, 
		 int32_t number_of_ends, fastq_t *reads, int32_t single_output, int32_t rend, int32_t prev, int64_t **output_counts)
{
	assert(rend == 1 || rend == 0);

	if (0 == bwa_output || NULL==output_prefix) {
		if (1 == single_output && 2 == number_of_ends) { 
			if (0 != strcmp(reads[0].name, reads[1].name)) {  // single output, read has only 1 end
				dump_read(fps[0], read, &(*output_counts)[0]);
			}
			else { // single output enabled, the read has two ends
				dump_read(fps[rend+1], read, &(*output_counts)[rend+1]);
			}
		}
		else { // No single_output enabled
			dump_read(fps[0], read, &(*output_counts)[0]);
		}
	}
	else if (1 == bwa_output && 2 == number_of_ends) { // bwa_output (two ends)

		//fprintf(stdout, "DRD: %d %d %d\n", rend, reads[0].to_print, reads[1].to_print);
		if (rend == 0 && reads[0].to_print == 1 && reads[1].to_print == 0) { // end 0, one read
			dump_read(fps[0], read, &(*output_counts)[0]); 
		}
		else if (rend == 0 && reads[0].to_print == 1 && reads[1].to_print == 1) { // end 0, two reads
			dump_read(fps[1], read, &(*output_counts)[1]);
		}
		else if (rend == 1 && reads[0].to_print == 0 && reads[1].to_print == 1 && prev == 1) { // end 1, two reads
			dump_read(fps[2], read, &(*output_counts)[2]);
		}
		else { // end 1, one read
			dump_read(fps[0], read, &(*output_counts)[0]); 
		}

		/*
		int i;
		int read_type = 0;

		// Find read type (read1 or read2)
		for (i=0; read->name[i]!='\0'; ++i);
		read_type = read->name[i-1];
		assert(read_type == 49 || read_type == 50);

		if (read_type == 50) {
			if (0 == cmp_read_names(reads[0].name, reads[1].name)) {
				dump_read(fps[2], read); // dump read1 -- will dump read2 in the next func call
			}
			else {
				dump_read(fps[0], read); // dump read1 in single
			}
		}

		if (read_type == 49) {
			if (0 == cmp_read_names(reads[0].name, reads[1].name)) {
				dump_read(fps[1], read);      // dump read2
			}
		}
		*/
	}
	else if (1 == bwa_output && 1 == number_of_ends) {
		dump_read(fps[0], read, &(*output_counts)[0]); // single
	}
}

void dump_read(AFILE *afp_output, fastq_t *read, int64_t *output_count)
{
	char at = '@';
	char plus = '+';
	char new_line = '\n';
	int32_t i;

	// Dump the read
	// Name
	AFILE_afwrite(&at, sizeof(char), 1, afp_output);
	for(i=0;i<strlen(read->name);i++) {
		AFILE_afwrite(&read->name[i], sizeof(char), 1, afp_output);
	}
	AFILE_afwrite(&new_line, sizeof(char), 1, afp_output);

	// Sequence
	for(i=0;i<strlen(read->read);i++) {
		AFILE_afwrite(&read->read[i], sizeof(char), 1, afp_output);
	}
	AFILE_afwrite(&new_line, sizeof(char), 1, afp_output);

	// Comment
	AFILE_afwrite(&plus, sizeof(char), 1, afp_output);
	AFILE_afwrite(&new_line, sizeof(char), 1, afp_output);

	// Quality
	for(i=0;i<strlen(read->qual);i++) {
		AFILE_afwrite(&read->qual[i], sizeof(char), 1, afp_output);
	}

	AFILE_afwrite(&new_line, sizeof(char), 1, afp_output);

	(*output_count)++;
}

void add_read_name_prefix(fastq_t *read, char *read_name_prefix, int rnp_len) {
	char tmp_name[SEQUENCE_NAME_LENGTH]="\0";

	strcpy(tmp_name, read->name);
	strcpy(read->name, read_name_prefix);  // read name prefix already has a ":" at the end
	strncpy(read->name + rnp_len, tmp_name, SEQUENCE_NAME_LENGTH - rnp_len); // strncpy fills the remaining space with zeros
}

/*
 * @427_67_118                --> @bwa:427_67_118/2
 * T3233100000020000000000000 --> GTTCAAAAAAGAAAAAAAAAAAAA
 * B@.7>/+-8:0.<8:/%@;280>=>  --> @.7>/+-8:0.<8:/%@;280>=> 
 */
void to_bwa(fastq_t *read, int32_t n_end)
{
	int32_t i;
	char tmp_read[SEQUENCE_LENGTH]="\0";
	char de[52]="\0";

	// Prepare the double encode convertion (TODO: make this global for speed)
	de[46]='N'; de[48]='A'; de[49]='C'; de[50]='G'; de[51]='T';

	if (n_end == 1) {
		strcat(read->name, "/2"); // In bwa's script, they encode R3 as 1 and F3 as 2.
	}
	else {
		strcat(read->name, "/1");
	}

	// Remove last base of the primer and first color call
	// and double encode the rest of the color space sequence
	strcpy(tmp_read, read->read);
	for(i=2; i<strlen(read->read); i++) {
		tmp_read[i-2] = de[(int)read->read[i]];
	}
	tmp_read[i-2] = '\0';
	strcpy(read->read, tmp_read);

	// Remove first quality value
	for(i=0; i<strlen(read->qual); i++) {
		read->qual[i]= read->qual[i+1];
	}
}

#ifdef HAVE_FSEEKO
void fastq_read(fastq_t *read, AFILE *afp_csfasta, AFILE *afp_qual, off_t csfasta_end_pos, off_t qual_end_pos, int32_t trim_end, int32_t bwa_output)
#else
void fastq_read(fastq_t *read, AFILE *afp_csfasta, AFILE *afp_qual, int32_t trim_end, int32_t bwa_output)
#endif
{
	char *FnName="fastq_read";
	char qual_name[SEQUENCE_NAME_LENGTH]="\0";
	char qual_line[SEQUENCE_NAME_LENGTH]="\0";
	int32_t qual[SEQUENCE_LENGTH];
	int32_t i;
	char *pch=NULL, *saveptr=NULL;

	assert(0 == read->is_pop);
	assert(0 == read->to_print);

	if(NULL == afp_csfasta &&
			NULL == afp_qual) {
		//fprintf(stderr, "return 1\n"); // HERE
		return;
	}

	// Read in
#ifdef HAVE_FSEEKO
	if(read_line(afp_csfasta, csfasta_end_pos, read->name) < 0 || 
			read_line(afp_csfasta, csfasta_end_pos, read->read) < 0 ||
			read_line(afp_qual, qual_end_pos, qual_name) < 0 ||
			read_line(afp_qual, qual_end_pos, qual_line) < 0) {
		//fprintf(stderr, "return 2\n"); // HERE
		return;
	}
#else
	if(read_line(afp_csfasta, read->name) < 0 || 
			read_line(afp_csfasta, read->read) < 0 ||
			read_line(afp_qual, qual_name) < 0 ||
			read_line(afp_qual, qual_line) < 0) {
		//fprintf(stderr, "return 2\n"); // HERE
		return;
	}
#endif

	StringTrimWhiteSpace(read->name);
	StringTrimWhiteSpace(read->read);
	StringTrimWhiteSpace(qual_name);
	StringTrimWhiteSpace(qual_line);

	// Parse qual line
	pch = strtok_r(qual_line, " ", &saveptr);
	for(i=0;pch!=NULL;i++) {
		qual[i] = atoi(pch);
		pch = strtok_r(NULL, " ", &saveptr);
	}

	// Check that the read name and csfasta name match
	if(0 != strcmp(read->name, qual_name)) {
		fprintf(stderr, "read->name=[%s]\nqual_name=[%s]\n", read->name, qual_name);
		PrintError(FnName, "read->name != qual_name", "Read names did not match", Exit, OutOfRange);
	}
	// Remove leading '@' from the read name
	for(i=1;i<strlen(read->name);i++) {
		read->name[i-1] = read->name[i];
	}
	read->name[i-1]='\0';

	// Convert SOLiD qualities
	for(i=0;i<strlen(read->read)-1;i++) {
		/*
		   fprintf(stderr, "%c -> %c\n%d -> %d\n", qual[i], (qual[i] <= 93 ? qual[i] : 93) + 33, qual[i], (qual[i] <= 93 ? qual[i] : 93) + 33);
		   exit(1);
		   */
		if (bwa_output == 0) {
			qual[i] = 33 + (qual[i] <= 0 ? 0 : (qual[i] <= 93 ? qual[i] : 93));
		}
		else { // bwa's solid2fastq does not force that max 93
			qual[i] = 33 + (qual[i] <= 0 ? 0 : qual[i]);
		}
		read->qual[i] = (char)qual[i];
	}

	// Trim last _R3 or _F3 or _whatever
	read_name_trim(read->name);

	// Trim last few colors
	if(0 < trim_end) {
		if(strlen(read->read) + 1 <= trim_end) {
			PrintError(FnName, "-t", "Trimming all the colors", Exit, OutOfRange);
		}
		read->read[strlen(read->read)-trim_end]='\0';
		read->qual[strlen(read->qual)-trim_end]='\0';
	}

	read->is_pop = 1;
}

int32_t cmp_read_names(char *name_one, char *name_two, int32_t rnp_len)
{
	char *name_one_cur = NULL;
	char *name_two_cur = NULL;
	int32_t name_one_num_state = 0, name_two_num_state = 0;
	int32_t return_value = 0;
	int32_t name_one_index = 0, name_two_index = 0;
	/*
	   fprintf(stderr, "comparing %s with %s\n", name_one, name_two);
	   */

	// for safety, skip the read name prefix
	if (rnp_len > 0) { // @read_name:5_123_456;
		name_one_cur = strtok_mod(name_one, ":", &name_one_index);
		name_two_cur = strtok_mod(name_two, ":", &name_two_index);
		free(name_one_cur);
		free(name_two_cur);
	}

	name_one_cur = strtok_mod(name_one, "_", &name_one_index);
	name_two_cur = strtok_mod(name_two, "_", &name_two_index);

	while(NULL != name_one_cur && NULL != name_two_cur) {
		/*
		   fprintf(stderr, "name_one_cur=%s\nname_two_cur=%s\n", name_one_cur, name_two_cur);
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
	int32_t l;

	// Trim last _R3 or _F3 : >427_67_118_R3
        // For V4 (PE), read2 file uses: F5-P2: >427_67_118_F5-P2
	if(NULL == name) {
		return;
	}

	l=strlen(name);
	if(3 < l &&
			name[l-3]=='_' &&
			(name[l-2]=='F' || name[l-2]=='R') &&
			name[l-1]=='3') { 
		name[l-3]='\0';
	}
	else if(3 < l &&
			name[l-3]=='-' &&
			name[l-2]=='P' &&
			name[l-4]=='5') {
		name[l-6]='\0';
        }
	else if(3 < l &&
			name[l-3]=='-' &&
			name[l-2]=='B' &&
			name[l-4]=='C') {
		name[l-6]='\0';
        }
        else {
            // try to trim off until the last underscore
            l--; // 0-based
            while(0 <= l) {
                if('_' == name[l]) {
                    name[l]='\0';
                    break;
                }
                name[l]='\0';
                l--;
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

#ifdef HAVE_FSEEKO
int32_t read_line(AFILE *afp, off_t end_pos, char *line)
#else
int32_t read_line(AFILE *afp, char *line)
#endif
{
	//char *FnName="read_line";
	char c=0;
	int32_t i=0, p=0;
	int32_t state=0;

	if(NULL == afp) return -1;

#ifdef HAVE_FSEEKO
	if(end_pos > 0 && AFILE_aftell(afp) >= end_pos) return -1;
#endif

	// States:
	// 0 - no characteres in line
	// 1 - reading valid line 
	// 2 - reading comment line

	// Read a character at a time
	while(0 != AFILE_afread(&c, sizeof(char), 1, afp)) {
		if(EOF == c) {
			if(0 < i) { // characters have been read
				line[i]='\0';
				//fprintf(stderr, "return read_line 1\n");
				return  1;
			}
			else { // nothing to report
				//fprintf(stderr, "return read_line 2\n");
				return -1;
			}
		}
		else if('\n' == c) { // endline
			if(1 == state) { // end of the valid line
				line[i]='\0';
				//fprintf(stderr, "return read_line 3\n");
				return 1;
			}
			i=state=0;
		}
		else {
			if(0 == state) { // first character in the line
				if('#' == c) { // comment
					state = 2;
				}
				else { // valid line
					state = 1;
					assert(i==0);
					assert(p==0);
				}
			}
			if(1 == state) { // valid line
				// if previous was not whitespace or 
				// current is not a whitespace
				if(0 == p || 0 == IsWhiteSpace(c)) { 
					line[i]=c;
					i++;
					p=IsWhiteSpace(c);
				}
			}
		}
	}
	// must have hit eof
	if(0 < i) { // characters have been read
		line[i]='\0';
		//fprintf(stderr, "return read_line 4\n");
		return  1;
	}
	return -1;
}

/*
 * Close open File descriptors
 * Notice that in bwa mode the single file is in afp_output[0] and the ends are in afp_output[1]/[2]
 * In single_output mode (BF) afp_output[0]/[1] have the ends and [2] has the single.
 * There is only 1 single file but many ends files.
 */
void close_fds(AFILE **afp_output, int32_t bwa_output, char *prefix_output, int32_t number_of_ends, int32_t single_output)
{
	assert(NULL != afp_output[0]);
	AFILE_afclose(afp_output[0]);
	
	if (prefix_output != NULL && number_of_ends == 2 && (1 == bwa_output || 1 == single_output)) {
		assert(NULL != afp_output[1]);
		assert(NULL != afp_output[2]);
		AFILE_afclose(afp_output[1]);
		AFILE_afclose(afp_output[2]);
	}
}

