#ifndef SRMA_SAM_IO_H_
#define SRMA_SAM_IO_H_

#include "node.h"
#include "samtools/sam.h"

/*! @typedef
  @abstract 
  @field  fp_i  the file index from which the bam structure was read
  @field  b     pointer to the bam structure
  @field  prev  pointer to the previous bam record in the list, if any 
  @field  next  pointer to the next bam record in the list, if any 
  @field  node  pointer to the start node in this alignment
*/
typedef struct __bam_record_t {
	int32_t fp_i;
	bam1_t *b;
	struct __bam_record_t *prev, *next;
	node_t *node;
} bam_record_t;

/*! @function
  @abstract           free memory associated with this record
  @param  bam_record  pointer to the bam record
*/
void bam_record_free(bam_record_t *bam_record);

/*! @typedef
  @abstract     structure to help implement a linked list of bam records
  @param  head  pointer to the head of the linked list
  @param  tail  pointer to the tail of the linked list
  @param  iter  pointer used to iterate through the list
  @param  size  the size of hte linked list
*/
typedef struct __bam_record_ll_t {
	bam_record_t *head, *tail, *iter;
	int32_t size;
} bam_record_ll_t;

/*! @function
  @abstract  initialize an empty linked list structure
  @return  a pointer to the initialized linked list structure
*/
bam_record_ll_t *bam_record_ll_init();

/*! @function
  @abstract     add the bam structure to the linked list structure
  @param  ll    pointer to the linked list structure
  @param  b     pointer to the bam record structure
*/
void bam_record_ll_add(bam_record_ll_t *ll, bam_record_t *b);

/*! @function
  @abstract
  @param  ll    pointer to the linked list structure
  @param  b     pointer to the bam structure
  @param  fp_i  the file pointer index from which the bam structure was read
  @param  node  pointer to the start node in this alignment
*/
void bam_record_ll_add1(bam_record_ll_t *ll, bam1_t *b, int32_t fp_i, node_t *node);

/*! @function
  @abstract
  @param  ll    pointer to the linked list structure
  @return       pointer to the bam record structure, NULL if the linked list is empty
*/
bam_record_t *bam_record_ll_remove(bam_record_ll_t *ll);

/*! @function
  @abstract     remove all records from this linked list structure
  @param  ll    pointer to the linked list structure
*/
void bam_record_ll_clear(bam_record_ll_t *ll);

/*! @function
  @abstract     remove all records from this linked list structure and destroy the linked list structure
  @param  ll    pointer to the linked list structure
  @param  b     pointer to the bam structure
  @param  fp_i  the file pointer index from which the bam structure was read
*/
void bam_record_ll_destroy(bam_record_ll_t *ll);

#define SRMA_SAM_IO_TYPE_SAM 0
#define SRMA_SAM_IO_TYPE_BAM 1
#define SRMA_SAM_IO_TYPE_BAM_ITER 2

/*! @typedef
  @abstract              the sam io structure for multiple inputs/outputs and range querying
  @field  fps_in          array of sam files from which to read
  @field  fps_in_type     type of read function to perform
  @field  fps_in_num      number of input sam files
  @field  fps_out         array of sam files to which to write
  @field  fps_out_num     number of output sam files
  @field  bam_indexes     bam indexes for eah BAM input; the value is NULL for a SAM input 
  @field  bam_iters       bam iterators for eah BAM input; the value is NULL for a SAM input 
  @field  buffer_head     pointer to the head of the bam record buffer
  @field  buffer_tail     pointer to the tail of the bam record buffer
*/
typedef struct {
	samfile_t **fps_in;
	int32_t *fps_in_type;
	int32_t fps_in_num;
	samfile_t **fps_out;
	int32_t fps_out_num;
	bam_index_t **bam_indexes;
	bam_iter_t *bam_iters;
	bam_record_ll_t *buffer;
} srma_sam_io_t;

/*! @function
  @abstract
  @param  fn_inputs         a list of input file names 
  @param  fn_inputs_num     the number of input file names
  @param  fn_outputs        a list of output file names
  @param  fn_outputs_num    the number of output file names, which must equal 1 or the number of input file names
  @param  fn_output_header  file name for the header to use when there are multiple inputs and a single output
  @param  use_ranges        0 if we are not to use ranges, 1 otherwise
  @discussion               if use_ranges is 1, then sam indexes will be loaded
  @return 
*/
srma_sam_io_t *srma_sam_io_init(char **fn_inputs, int32_t fn_inputs_num, char **fn_outputs, int32_t fn_outputs_num, char *fn_output_header, int32_t use_ranges);

/*! @function
  @abstract   
  @param  s  the pointer to the io stucture
  @return    1 if there are more bam records, 0 otherwise 
*/
inline int32_t srma_sam_io_has_next(srma_sam_io_t *s);

/*! @function
  @abstract
  @param  s  the pointer to the io stucture
  @return    a pointer to the next bam record stucture, NULL if none exists
*/
bam_record_t *srma_sam_io_get_next(srma_sam_io_t *s);

/*! @function
  @abstract  initialize the buffer
  @param  s  the pointer to the io stucture
*/
void srma_sam_io_init_buffer(srma_sam_io_t *s);

/*! @function
  @abstract               move to a new range to process within the bam files
  @param  s               the pointer to the io stucture
  @param  sequence_index  0-based contig index (tid)
  @param  position_start  0-based start position
  @param  position_end    0-based end position
  @discussion             only works for BAM files.
*/
void srma_sam_io_change_range(srma_sam_io_t *s, int32_t sequence_index, int32_t position_start, int32_t position_end);

/*! @function
  @abstract  destroys the input iterators and buffer
  @param  s  the pointer to the io stucture
*/
void srma_sam_io_close_inputs(srma_sam_io_t *s);

/*! @function
  @abstract  writes the bam associated with the pointer record to file
  @param  s  the pointer to the io stucture
  @param  b  the pointer to the program record
  @discussion  the output file index is stored within the bam record
*/
void srma_sam_io_write(srma_sam_io_t *s, bam_record_t *b);

/*! @function
  @abstract  destroy the io structure
  @param  s  the pointer to the io stucture
  @discussion  this also closes the input/output/index files
*/
void srma_sam_io_free(srma_sam_io_t *s);

#endif
