#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <zlib.h>
#include <stdint.h>
#include <limits.h>
#include "BError.h"
#include "util.h"

void *my_realloc(void *ptr, size_t size, char *fn_name)
{
	ptr = realloc(ptr, size);
	if(NULL == ptr) {
    fprintf(stderr, "my_realloc: requesting: %ld\n", size);
		PrintError(fn_name, NULL, "Could not reallocate memory", Exit, ReallocMemory);
	}
	return realloc(ptr, size);
}

void *my_malloc(size_t size, char *fn_name)
{
	void *ptr = malloc(size);
	if(NULL == ptr) {
		PrintError(fn_name, NULL, "Could not allocate (malloc) memory", Exit, MallocMemory);
	}
	return ptr;
}

void *my_calloc(size_t num, size_t size, char *fn_name)
{
	void *ptr = calloc(num, size);
	if(NULL == ptr) {
		PrintError(fn_name, NULL, "Could not allocate (calloc) memory", Exit, MallocMemory);
	}
	return ptr;
}
