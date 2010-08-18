#ifndef SW_ALIGN_H_
#define SW_ALIGN_H_

#include "graph.h"
#include "sw_heap.h"

/*! @function
  @abstract                    perform re-alignment of a bam structure to the graph structure
  @param  g                    pointer to the graph structure
  @param  b                    pointer to the original bam structure
  @param  n                    pointer to the starting node of the original alignment
  @param  heap                 pointer to a alignment heap
  @param  rg_id                the program record id
  @param  offset               number of bases around the original alignment to consider
  @param  cutoffs              for filtering nodes
  @param  correct_bases        update the bases to match the alignment (no effect for color space data)
  @param  use_qualities        use quality scores to weight the alignment
  @param  max_total_coverage   the maximum total coverage over a given position to consider
  @param  max_heap_size        the maximum number of nodes allowed in the heap
  @return                      a pointer to the bam structure
  @discussion                  tries to re-align the given BAM.  If unsuccessful, the pointer to the original bam structure will be returned, otherwise a pointer to a new bam structure is returned.  In the latter case, the original bam structure is destroyed.
*/
bam1_t *sw_align(graph_t *g, bam1_t *b, node_t *n, sw_heap_t *heap, char *rg_id, int32_t offset, cov_cutoffs_t *cutoffs, uint8_t correct_bases, uint8_t use_qualities, int32_t max_total_coverage, int32_t max_heap_size);

/*! @function
  @abstract
  @param  bam_old          pointer to the original bam structure
  @param  rg_id            the program record id
  @param  heap             the heap from which the best node comes 
  @param  sw_node_best_i   the index of the best psuedo-cell from which to backtrace
  @param  space            the encoding space
  @param  colors           a string of colors including adaptor, NULL if not color space
  @param  color_qualities  a string of qualities for the colors, NULL if not color space
  @param  strand           0 for forward, 1 for reverse
  @param  correct_bases    update the bases to match the alignment (no effect for color space data)
  @return                  a pointer to the bam structure
  @discussion              a pointer to the original bam structure if no best node was found, otherwise a pointer to a new bam structure based on the best alignment.  In the lastter case, the original bam is destroyed. 
*/
bam1_t *sw_align_update_bam(bam1_t *bam_old, char *rg_id, sw_heap_t *heap, int32_t sw_node_best_i, uint8_t space, char *colors, char *color_qualities, uint8_t strand, uint8_t correct_bases);

#endif
