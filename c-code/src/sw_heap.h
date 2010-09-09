#ifndef SW_HEAP_H_
#define SW_HEAP_H_

#include "sw_node.h"

/*! @typedef
  @abstract            queue structure for sw_node_t structures
  @field  type          type of queue (min/max)
  @field  queue         array of pointers to nodes in the queue
  @field  queue_start   the index of the first node in the queue
  @fiel   queue_end     the index of the last node in the queue
  @field  queue_mem     number of node pointers allocated 
  @field  nodes         memory pool of nodes
  @field  nodes_mem     number of nodes allocated
  @field  nodes_next    pointer the next unused node
*/ 
typedef struct {
	int32_t type;
	int32_t *queue;
	int32_t queue_start, queue_end;
	int32_t queue_mem;

	sw_node_t *nodes;
	int32_t nodes_mem, nodes_next;
} sw_heap_t;

/*! @typedef
  @abstract     initialize the queue structure
  @param  type  the type of queue (min/max)
  @param  size  size of the intial heap
  @return       a pointer to the initialized memory strucutre 
*/
sw_heap_t *sw_heap_init(int32_t type, int32_t size);

/*! @function
  @abstract   get a node structure
  @param   h  a pointer to the heap
  @return     the index of an unused node
*/
int32_t sw_heap_get_node_i(sw_heap_t *h);

/*! @function
  @abstract      add the node to the heap
  @param  h      a pointer to the heap
  @param  cur_i  the index of the node to add
*/
void sw_heap_add_i(sw_heap_t *h, int32_t cur_i);

/*! @function
  @abstract  poll a node in the queue
  @param  h  a pointer to the heap
  @return    the index of the next node in the queue
  @discussion  the node is removed from the queue
*/
int32_t sw_heap_poll_i(sw_heap_t *h);

/*! @function
  @abstract  peek at a node in the queue
  @param  h  a pointer to the heap
  @return    the index of the next node in the queue
  @discussion  the node is not removed from the queue
*/
int32_t sw_heap_peek_i(sw_heap_t *h);

/*! @function
  @abstract  reset the queue structure
  @param  h  a pointer to the heap
  @discussion  this removes all nodes in the queue, but does not clear past nodes
*/
void sw_heap_reset(sw_heap_t *h);

/*! @function
  @abstract  clear the queue structure
  @param  h  a pointer to the heap
  @discussion  this removes all nodes in the queue and clears past nodes
*/
void sw_heap_clear(sw_heap_t *h);

/*! @function
  @abstract  clear all memory associated with this heap
  @param  h  a pointer to the heap
*/
void sw_heap_free(sw_heap_t *h);

#endif
