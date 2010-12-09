#ifndef GRAPH_H_
#define GRAPH_H_

#include "samtools/sam.h"
#include "node.h"

/*! @typedef
  @abstract              a graph of nodes that span some part of a contig
  @field  contig          the current contig (1-based)
  @field  position_start  the start position of the graph (1-based)
  @field  position_end    the end position of the graph (1-based)
  @field  nodes_mem       the memory allocated end position of the graph 
  @field  nodes           a list of list of nodes
  @field  is_empty        are there any nodes in this graph?
*/
typedef struct {
	uint32_t contig;
	uint32_t position_start, position_end;
	uint32_t nodes_mem;
	node_list_t **nodes;
	uint8_t is_empty;
} graph_t;

/*! @function
  @abstract       initialize a new graph
  @return         a pointer to the initialized graph structure
*/
graph_t *graph_init();

/*! @function
  @abstract            add the alignment to the graph 
  @param  g            pointer to the graph
  @param  b            pointer to the bam structure
  @param  ref          a pointer the the reference sequence structure
  @param  use_threads  1 if this function should be thread-safe, 0 otherwise
  @return              pointer the first node in the alignment
  @discussion          this function is thread-safe
*/
node_t *graph_add_sam(graph_t *g, bam1_t *b, ref_t *ref, int32_t use_threads);

/*! @function
  @abstract
  @param  g            pointer to the graph
  @param  node         pointer to the node to be added
  @param  prev         pointer to the previous node
  @param  use_threads  1 if this function should be thread-safe, 0 otherwise
  @return              pointer to the node in the graph
  @discussion          if node is already in the graph, it will be freed
*/
node_t *graph_add_node(graph_t *g, node_t *node, node_t *prev, int32_t use_threads);

/*! @function
  @param  g     pointer to the graph
  @param  node  pointer to the node to be tested
  @return       the node in the graph, NULL otherwise
*/
node_t *graph_contains(graph_t *g, node_t *node);

/*! @function
  @param  g         pointer to the graph
  @param  position  minimum position to be returned
  @return           return the 0-based index of the first node with position greater than or equal to 'position' 
*/
uint32_t graph_get_node_list_index_at_or_after(graph_t *g, uint32_t position);

/*! @function
  @param  g         pointer to the graph
  @param  position  maximum position to be returned
  @return           return the 0-based index of the first node with position less than or equal to 'position' 
*/
uint32_t graph_get_node_list_index_at_or_before(graph_t *g, uint32_t position);

/*! @function
  @param  g         pointer to the graph
  @param  position  the position of the node list
  @return           the node list structure associated with this position
*/
node_list_t *graph_get_node_list(graph_t *g, uint32_t position);

/*! @function
  @param  g         pointer to the graph
  @param  position  the position of the coverage value
  @return           the coverage value associated with this position
*/
uint16_t graph_get_coverage(graph_t *g, uint32_t position);

/*! @function
  @param  g                pointer to the graph
  @param  contig_index     the alignment contig index (0-based)
  @param  alignment_start  the start position of the alignment (1-based)
  @param  offset           the offset of the alignment algorithm
*/
void graph_prune(graph_t *g, uint32_t contig_index, uint32_t alignment_start, int32_t offset);

/*! @function
  @param  g                pointer to the graph
*/
void graph_free(graph_t *g);

#endif
