#ifndef SRMA_UTIL_H_
#define SRMA_UTIL_H_

#include <stdint.h>
#include "samtools/bam.h"
#include "srma_faidx.h"

#define ALN_GAP '-'
#define SRMA_SPACE_NT 0
#define SRMA_SPACE_CS 1
#define SRMA_SW_HEAP_MIN 0
#define SRMA_SW_HEAP_MAX 1
#define SRMA_CORRECT_BASE_QUALITY_PENALTY 20

extern char *bam_nt16_rev_table;

/*! @typedef
  @abstract                a user-friendly alignment structure
  @field  read              the read, with gaps
  @field  ref				the reference, with gaps
  @field  positions         for each base in the read (without gaps), the index of the corresponding base in the reference (without gaps)
  @field  positions_index   for each base in the read (without gaps), the index of the corresponding bases in the alignment
  @field  length            the length of the read and ref arrays
  @field  positions_length  the length of the positions array
  */
typedef struct {
	uint8_t *read;
	uint8_t *ref;
	uint32_t *positions;
	uint32_t *positions_index;
	uint32_t length, positions_length;
} bam_aln_t;

/*! @function
  @abstract    initialize the alignment structure
  @param  bam  the SAM/BAM structure to convert
  @param  ref  a pointer to the reference sequence structure
  @return      a pointer to the initialized alignment structure
  */
bam_aln_t *bam_aln_init(bam1_t *bam, ref_t *ref);

/*! @function
  @abstract    free memory associated with this structure
  @param  aln  pointer to the alignment structure
  */
void bam_aln_free(bam_aln_t *aln);

/*! @function
  @abstract     convert ascii quality to integer quality
  @param  qual  the ascii quality to convert
  @return       the integer quality value, bounded by 0 and 255
  */
inline int32_t srma_char2qual(char qual);

/*! @function
  @abstract     convert integer quality to ascii quality
  @param  qual  the integer quality to convert
  @return       the ascii quality value
  */
inline char srma_qual2char(int32_t qual);

/*! @typedef
  @abstract                   structure to hold allele cut-offs
  @field  min_allele_coverage  array of coverage cut-offs
  @field  max_coverage         the maximum coverage for which a value is stored
  @field  length               the length of a cut-off array
*/
typedef struct {
	int32_t *min_allele_coverage;
	int32_t max_coverage, length;
} cov_cutoffs_t;

/*! @abstract
  @param min_allele_coverage  the absolute minimum allele coverage regardless of total coverage
  @param min_allele_prob      the minimum binomial probability
  @return                     a pointer to the initialzed cut-off structure
*/
cov_cutoffs_t *cov_cutoffs_init(int32_t min_allele_coverage, double min_allele_prob);

/*! @abstract       for a given total coverage, the minumum allele coverage 
  @param  c         a pointer to the cut-off structure
  @param  coverage  the total coverage
  @return           the minumum allele coverage 
*/
int32_t cov_cutoffs_get(cov_cutoffs_t *c, int32_t coverage);

/*! @abstract
  @param  c  a pointer to the cut-off structure
*/
void cov_cutoffs_free(cov_cutoffs_t *c);

#ifndef roundup32
/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.  This was pulled from "samtools/bam.h"
*/
#define roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#endif
