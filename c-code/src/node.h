#ifndef NODE_H_
#define NODE_H_

#include <stdint.h>

/*! @typedef
  @abstract        the node structure
  @field  base      the type (upper 4 bits) and base (lower 4 bits)
  @field  contig    the one-based contig index
  @field  position  the one-based position
  @field  offset    the offset of this node if an insertion, zero otherwise
  @field  coverage  the number of alignments that pass through this node
  @field  next      the list of downstream nodes
  @field  prev      the list of upstream nodes
*/
typedef struct {
	uint8_t base; // lower-4 base, upper-4 type
	uint32_t contig, position; // one-based
	uint8_t offset; // for insertions
	uint16_t coverage; 
	struct __edge_list_t *next, *prev; // NULL terminated list
} node_t;

/*! @macro
  @abstract      get the node's type
  @param  _node  pointer to the node
  @return        the nodes type
*/
#define __node_type(_node) (((_node->base) >> 4))

/*! @macro
  @abstract      get the node's base
  @param  _node  pointer to the node
  @return        the nodes base 
*/
#define __node_base(_node) (((_node->base) & 0x0f))

/*! @enum
  @abstract node types
*/
enum {
	NODE_MATCH     = 0,
	NODE_MISMATCH  = 1,
	NODE_INSERTION = 2,
	NODE_DELETION  = 3
};

/*! @function
  @abstract         initialize a node structure
  @param  base      the base to be represented by the node
  @param  type      the type of node
  @param  contig    the one-based contig index
  @param  position  the one-based position
  @param  offset    the offset of this node if an insertion, zero otherwise
  @param  prev      the connecting upstream node
  @return
*/
node_t *node_init(char base, uint8_t type, uint32_t contig, uint32_t position, node_t *prev);

/*! @function
  @abstract   compares two nodes based on contig, position, offset, type, and base, in that order
  @param  n1  the first node to compare
  @param  n2  the second node to compare
  @return     -1 if the first node is less than the second node, 0 if equal, 1 otherwise
*/
inline int node_compare(node_t *n1, node_t *n2);

/*! @function
  @abstract   compares two nodes based on offset, type, and base, in that order
  @param  n1  the first node to compare
  @param  n2  the second node to compare
  @return     -1 if the first node is less than the second node, 0 if equal, 1 otherwise
  @discussion assumes contig and position are qual for n1 and n2
*/
inline int node_compare2(node_t *n1, node_t *n2);

/*! @function
  @abstract     add the node to the current node's downstream list
  @param  cur   the node whos list should be updated
  @param  node  the node to be added
*/
inline void node_add_to_next(node_t *cur, node_t *node);

/*! @function
  @abstract     add the node to the current node's upstream list
  @param  cur   the node whos list should be updated
  @param  node  the node to be added
*/
inline void node_add_to_prev(node_t *cur, node_t *node);

/*! @function
  @abstract     remove the node to the current node's downstream list
  @param  cur   the node whos list should be updated
  @param  node  the node to be removed
*/
inline void node_remove_from_next(node_t *cur, node_t *node);

/*! @function
  @abstract     remove the node to the current node's upstream list
  @param  cur   the node whos list should be updated
  @param  node  the node to be removed
*/
inline void node_remove_from_prev(node_t *cur, node_t *node);


/*! @function
  @abstract     free memory associated with this node
  @param  node  the node to be freed
*/
void node_free(node_t *node);

/*! @typedef
  @abstract         a structure for representing edges in the graph 
  @field  coverages  the coverages to the nodes
  @field  nodes      pointers to the nodes
  @field  length     the length of the list
*/
typedef struct __edge_list_t {
	uint16_t *coverages;
	node_t **nodes;
	int32_t length;
} edge_list_t;

/*! @function
  @abstract  initialize an empty list
  @return    a pointer to initialized memory
*/
inline edge_list_t *edge_list_init();

/*! @function
  @abstract     free memory associated with this list
  @param  list  a pointer to the list to be freed
*/
void edge_list_free(edge_list_t *list);

/*! @function
  @abstract     insert a node into the list
  @param  list  a pointer to the list 
  @param  node  the node to be inserted
*/
void edge_list_add(edge_list_t *list, node_t *node);

/*! @function
  @abstract     remove a node into the list
  @param  list  a pointer to the list 
  @param  node  the node to be inserted
  @return       1 if successful, 0 otherwise
*/
int32_t edge_list_remove(edge_list_t *list, node_t *node);

/*! @typedef
  @abstract         a structure for a priority queue of nodes
  @field  nodes      pointers to the nodes
  @field  length     the length of the list
*/
typedef struct __node_list_t {
	node_t **nodes;
	int32_t length;
} node_list_t;

/*! @function
  @abstract  initialize an empty list
  @return    a pointer to initialized memory
*/
inline node_list_t *node_list_init();

/*! @function
  @abstract     free memory associated with this list
  @param  list  a pointer to the list to be freed
*/
void node_list_free(node_list_t *list);

/*! @function
  @abstract     free memory associated with this list, but not the structure itself
  @param  list  a pointer to the list to be freed
*/
void node_list_clear(node_list_t *list);

/*! @function
  @abstract     insert a node into the list
  @param  list  a pointer to the list 
  @param  node  the node to be inserted
*/
void node_list_add(node_list_t *list, node_t *node);

#endif
