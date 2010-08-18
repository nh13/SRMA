#ifndef DAT_ALLOC_H_
#define DAT_ALLOC_H_

/*! @function
  @abstract              wrapper function for malloc
  @param  size           the size of the memory block, in bytes
  @param  fn_name        the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return                upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
*/
inline void *srma_malloc(size_t size, const char *fn_name, const char *variable_name);

/*! @function
  @abstract              wrapper function for realloc
  @param  ptr            the pointer to a memory block previously allocated
  @param  size           the size of the memory block, in bytes
  @param  fn_name        the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return 				 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  @discussion			 the ptr must be a memory block previously allocated with malloc, calloc, or realloc to be reallocated; if the ptr is NULL, a new block of memory will be allocated. 
*/
inline void *srma_realloc(void *ptr, size_t size, const char *fn_name, const char *variable_name);

/*! @function
  @abstract  wrapper function for calloc
  @param  num            the number of elements to be allocated
  @param  size           the size of the memory block, in bytes
  @param  fn_name        the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return                upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
*/
inline void *srma_calloc(size_t num, size_t size, const char *fn_name, const char *variable_name);

/*! @function
  @abstract        wrapper for 'strdup' that checks memory allocation
  @param  str      string to be copied
  @param  fn_name  the calling function name 
  @return          a pointer to the copied string
*/
inline char *srma_strdup(const char *str, const char *fn_name);

#endif
