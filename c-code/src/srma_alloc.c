#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "srma_error.h"
#include "srma_alloc.h"

inline void *srma_malloc(size_t size, const char *fn_name, const char *variable_name)
{
	void *ptr = malloc(size);
	if(NULL == ptr) {
		srma_error(fn_name, variable_name, Exit, MallocMemory);
	}
	return ptr;
}

inline void *srma_realloc(void *ptr, size_t size, const char *fn_name, const char *variable_name)
{
	ptr = realloc(ptr, size);
	if(NULL == ptr) {
		srma_error(fn_name, variable_name, Exit, ReallocMemory);
	}
	return ptr;
}

inline void *srma_calloc(size_t num, size_t size, const char *fn_name, const char *variable_name)
{
	void *ptr = calloc(num, size);
	if(NULL == ptr) {
		srma_error(fn_name, variable_name, Exit, MallocMemory);
	}
	return ptr;
}

inline char *srma_strdup(const char *str, const char *fn_name)
{
	char *s = strdup(str);
	if(NULL == s) {
		srma_error(fn_name, str, Exit, MallocMemory);
	}
	return s;
}
