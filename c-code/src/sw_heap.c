#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>

#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_util.h"
#include "sw_node.h"
#include "sw_heap.h" 

sw_heap_t *sw_heap_init(int32_t type, int32_t size)
{
	sw_heap_t *h = NULL;

	h = srma_malloc(sizeof(sw_heap_t), __func__, "h");

	h->type = type;
	h->queue_start = 0;
	h->queue_end = -1;
	h->queue_mem = size;
	h->queue = srma_calloc(h->queue_mem, sizeof(int32_t), __func__, "h->queue");

	h->nodes_mem = size;
	h->nodes_next = 0;
	h->nodes = srma_calloc(h->nodes_mem, sizeof(sw_node_t), __func__, "h->nodes");

	return h;
}

int32_t sw_heap_get_node_i(sw_heap_t *h)
{
	int32_t prev_mem, i;
	if(h->nodes_mem < h->nodes_next+1) {
		prev_mem = h->nodes_mem;
		h->nodes_mem = h->nodes_next+1;
		roundup32(h->nodes_mem); // from bam.h->nodes_mem;
		h->nodes = srma_realloc(h->nodes, sizeof(sw_node_t)*h->nodes_mem, __func__, "h->nodes");
		for(i=prev_mem;i<h->nodes_mem;i++) {
			// just to be safe
			h->nodes[i].node = NULL; 
			h->nodes[i].read_offset = 0;
			h->nodes[i].score = 0;
			h->nodes[i].coverage_sum = 0;
			h->nodes[i].start_position = 0;
			h->nodes[i].cur_i = -1;
			h->nodes[i].prev_i = -1;
		}
	}
	assert(h->nodes_next < h->nodes_mem); // DEBUG
	h->nodes[h->nodes_next].cur_i = h->nodes_next;
	h->nodes[h->nodes_next].node = NULL;
	h->nodes[h->nodes_next].prev_i = -1;
	h->nodes_next++;

	return (h->nodes_next-1);
}

static inline int sw_heap_get_i(sw_heap_t *h, int32_t cur_i)
{
	int32_t low, mid, high;

	if(h->queue_end < h->queue_start) {
		return h->queue_start;
	}

	low = h->queue_start;
	high = h->queue_end;
	mid = (low + high) / 2;
	while(low <= high) {
		mid = (low + high) / 2;
		int c = sw_node_compare(&h->nodes[mid], &h->nodes[cur_i], h->type);
		//fprintf(stderr, "low=%d mid=%d high=%d c=%d\n", low, mid, high, c);
		if(0 == c) {
			//fprintf(stderr, "returning 1 mid=%d\n", mid);
			return mid;
		}
		else if(c < 0) {
			high = mid-1;
		}
		else {
			low = mid+1;
		}
	}
	//fprintf(stderr, "returning 2 mid=%d\n", mid);
	return mid;
}

void sw_heap_add_i(sw_heap_t *h, int32_t cur_i)
{
	int32_t i, j;

	// binary search 
	i = sw_heap_get_i(h, cur_i);

	// check if there is enough memory 
	h->queue_end++;
	if(h->queue_mem <= h->queue_end) {
		h->queue_mem = h->queue_end;
		roundup32(h->queue_mem);
		h->queue = srma_realloc(h->queue, sizeof(sw_node_t*)*h->queue_mem, __func__, "h->queue");
	}
	// shift 
	for(j=h->queue_end;i<j;j--) {
		h->queue[j] = h->queue[j-1];
	}
	h->queue[i] = cur_i;
}

int32_t sw_heap_poll_i(sw_heap_t *h)
{
	int32_t ret = -1;
	//fprintf(stderr, "POLL h->queue_start=%d h->queue_end=%d h->queue_mem=%d\n", h->queue_start, h->queue_end, h->queue_mem); // DEBUG
	if(h->queue_end < h->queue_start) {
		return -1;
	}
	ret = h->queue[h->queue_start];
	h->queue_start++;
	return ret;
}

int32_t sw_heap_peek_i(sw_heap_t *h)
{
	if(h->queue_end < h->queue_start) {
		return -1;
	}
	return h->queue[h->queue_start];
}

void sw_heap_reset(sw_heap_t *h)
{
	// remove all from the top of the heap
	// do not destroy previous nodes in the heap
	while(-1 != sw_heap_poll_i(h)) {
	}	
}

void sw_heap_clear(sw_heap_t *h)
{
	int32_t i;
	for(i=0;i<h->queue_mem;i++) {
		h->queue[i] = -1;
	}
	h->queue_start = 0;
	h->queue_end = -1;
	for(i=0;i<h->nodes_mem;i++) {
		// just to be safe
		h->nodes[i].node = NULL; 
		h->nodes[i].read_offset = 0;
		h->nodes[i].score = 0;
		h->nodes[i].coverage_sum = 0;
		h->nodes[i].start_position = 0;
		h->nodes[i].prev_i = -1;
		h->nodes[i].cur_i = -1;
	}
	h->nodes_next = 0;
}

void sw_heap_free(sw_heap_t *h)
{
	free(h->queue);
	free(h->nodes);
	free(h);
}
