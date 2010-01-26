#ifndef AFLIB_H_
#define AFLIB_H_

enum {AFILE_NO_COMPRESSION=0, AFILE_BZ2_COMPRESSION, AFILE_GZ_COMPRESSION};
enum {AFILE_BZ2_READ=0, AFILE_BZ2_WRITE};

typedef struct {
	FILE *fp;
	BZFILE *bz2;
	gzFile gz;
	int32_t c;

	// bzip2
	char unused[BZ_MAX_UNUSED];
	int32_t n_unused, bzerror, open_type;
} AFILE;

AFILE *AFILE_afopen(const char* path, const char *mode, int32_t compression);
AFILE *AFILE_afdopen(int fildes, const char *mode, int32_t compression);
void AFILE_afclose(AFILE *afp); 
size_t AFILE_afread(void *ptr, size_t size, size_t count, AFILE *afp);
size_t AFILE_afwrite(void *ptr, size_t size, size_t count, AFILE *afp); 

#endif
