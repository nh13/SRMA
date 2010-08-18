#ifndef SW_NODE_H_
#define SW_NODE_H_

#include "node.h"

extern uint8_t nt2int_table[256];

/*! @typedef
  @abstract               structure for a pseudo-cell in the dynamic programming model
  @param  node            the node of the last base in this alignment
  @param  read_offset     number of read bases used in this alignment
  @param  score           score of alignment
  @param  coverage_sum    sum of coverage of alignment
  @param  start_position  start position of the alignment (1-based)
  @param  prev_i          index previous pseudo-cell, if any, otherwise -1
*/
typedef struct {
	node_t *node;
	uint32_t read_offset;
	int32_t score;
	uint32_t coverage_sum;
	uint32_t start_position;
	int32_t cur_i; // IMPORTANT: must be set in sw_heap_t
	int32_t prev_i;
} sw_node_t;

/*! @function
  @abstract               initialize and extend the alignment
  @param  sw_node_cur     pointer to the pseudo-cell to initialize
  @param  sw_node_prev    previous pseudo-cell in the alignment, -1 if none
  @param  cur_node        the node of the last base in the alignment
  @param  coverage        coverage to the current base
  @param  base            read base (integer)
  @param  qual            read base quality
  @param  use_qualities   1 if we are to score with qualities, 0 otherwise 
  @param  space           encoding space
  @return                 pointer to an initialized pseudo-cell
*/
void sw_node_init(sw_node_t *sw_node, sw_node_t *sw_node_prev, node_t *cur_node, int32_t coverage, uint8_t base, char qual, uint8_t use_qualities, uint8_t space);

/*! @function
  @abstract     compare to pseudo-cells based on heap type
  @param  n1    the first node to compare
  @param  n2    the second node to compare
  @param  type  the type of heap (max/min)
  @return       will return -1, 0, 1 based on genomic coordinate, read offset, node type, node base, and score, in that order, depending on the heap type
*/
int sw_node_compare(sw_node_t *n1, sw_node_t *n2, int32_t type);

/*! @function
  @abstract     free a pseudo-cell
  @param  node  the pseudo-cell to free
*/
void sw_node_free(sw_node_t *node);


#endif
