#ifndef SRMA_FAIDX_H_
#define SRMA_FAIDX_H_

#include "samtools/faidx.h"


/*! @typedef
  @abtract
  @param  ref  pointer to the reference sequence
  @param  tid  the contig (0-based)
  @param  beg  the start position (0-based)
  @param  end  the end position (0-based)
*/
typedef struct {
	    char *ref;
		int32_t tid, beg, end;
} ref_t;

/*! @typedef
  @abtract    a structure for holding a range to be used with a FASTA index
  @field  tid  the contig (0-based)
  @field  beg  the start position (0-based)
  @field  end  the end position (0-based)
*/
typedef struct {
	int32_t tid;
	int32_t beg, end;
} srma_range_t;

/*! @typedef
  @abtract       a structure for holding ranges to be used with a FASTA index
  @field  i       index of the current range
  @field  length  the number of ranges
  @field  ranges  an array of ranges
*/
typedef struct {
	int32_t i;
	int32_t length;
	srma_range_t *ranges;
} srma_ranges_t;

/*! @function
  @abstract          initialize the ranges structure
  @param  fai        pointer to the FASTA index structure
  @param  fn_ranges  a file name storing a list of ranges
  @param  range      the string holding the range
  @param  offset     the number of bases to expand this range 
  @param             pointer to the ranges structure
  @discussion        at most one of fn_ranges and ranges should be specified
*/
srma_ranges_t *srma_ranges_init(const faidx_t *fai, char *fn_ranges, char *range, int32_t offset);

/*! @function
  @abstract    add to the the current ranges structure 
  @param  r    pointer to the ranges structure
  @param  tid  the contig (0-based)
  @param  beg  the start position (0-based)
  @param  end  the end position (0-based)
*/
void srma_ranges_add(srma_ranges_t *r, int32_t tid, int32_t beg, int32_t end);

/*! @function
  @abstract    initialize the range structure from the fai structure
  @param  fai  pointer to the the FASTA index structure
  @param       pointer to the ranges structure
  @discussion  all bases will be included
*/
srma_ranges_t *srma_ranges_init_from_fai(faidx_t *fai);

/*! @function
  @abstract  get the current range structure 
  @param  r  pointer to the ranges structure
  @return    pointer to the next range structure, NULL if all have already been returned
*/
inline srma_range_t *srma_ranges_peek(srma_ranges_t *r);

/*! @function
  @abstract  get the next range structure 
  @param  r  pointer to the ranges structure
  @return    pointer to the next range structure, NULL if all have already been returned
*/
inline srma_range_t *srma_ranges_poll(srma_ranges_t *r);

/*! @function
  @abstract  free memory associated with this ranges structure
  @param  r  pointer to the ranges structure
*/
void srma_ranges_free(srma_ranges_t *r);

/*! @function
  @abstract    retrieve a reference sequence
  @param  fai  pointer to the the FASTA index structure
  @param  tid  the contig (0-based)
  @param  beg  the start position (0-based)
  @param  end  the end position (0-based)
*/
void srma_fai_fetch(const faidx_t *fai, ref_t *ref, int32_t tid, int32_t beg, int32_t end);

/*! @function
  @abstract    retrieve the contig name
  @param  fai  pointer to the the FASTA index structure
  @param  tid  the contig (0-based)
*/
inline char *srma_fai_name(faidx_t *fai, int tid);

#endif
