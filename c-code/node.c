#include <stdlib.h>
#include <stdio.h>
#include "srma_error.h"
#include "srma_alloc.h"
#include "node.h"

// acgtnACGTN
// 01234
uint8_t nt2int_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

node_t *node_init(char base, uint8_t type, uint32_t contig, uint32_t position, node_t *prev)
{
	node_t *node = NULL;

	node = srma_malloc(sizeof(node_t), __func__, "node");

	node->base = (type << 4) | nt2int_table[(int)base];
	node->contig = contig;
	node->position = position;
	node->offset = 0;
	node->coverage = 1;
	node->next = edge_list_init();
	node->prev = edge_list_init();

	if(NULL != prev) {
		if(NODE_INSERTION == __node_type(prev) && NODE_INSERTION == __node_type(node)) {
			node->offset = prev->offset + 1;
		}
	}

	return node;
}

inline int node_compare(node_t *n1, node_t *n2)
{
	if(n1->contig == n2->contig
			&& n1->position == n2->position
			&& n1->offset == n2->offset
			&& __node_type(n1) == __node_type(n2)
			&& __node_base(n1) == __node_base(n2))
	{
		return 0;
	}
	else if(n1->contig < n2->contig
			|| (n1->contig == n2->contig
				&& n1->position < n2->position)
			|| (n1->contig == n2->contig
				&& n1->position == n2->position
				&& n1->offset < n2->offset)
			|| (n1->contig == n2->contig
				&& n1->position == n2->position
				&& n1->offset == n2->offset
				&& __node_type(n1) < __node_type(n2))
			|| (n1->contig == n2->contig
				&& n1->position == n2->position
				&& n1->offset == n2->offset
				&& __node_type(n1) == __node_type(n2)
				&& __node_base(n1) < __node_base(n2)))
	{
		return -1;
	}
	else {
		return 1;
	}

}	

inline int node_compare2(node_t *n1, node_t *n2)
{
	if(n1->contig != n2->contig
			|| n1->position != n2->position) {
		srma_error(__func__, "Control reach unexpected point", Exit, OutOfRange);
	}

	if(n1->offset == n2->offset
			&& __node_type(n1) == __node_type(n2)
			&& __node_base(n1) == __node_base(n2))
	{
		return 0;
	}
	else if(n1->offset < n2->offset
			|| (n1->offset == n2->offset
				&& __node_type(n1) < __node_type(n2))
			|| (n1->offset == n2->offset
				&& __node_type(n1) == __node_type(n2)
				&& __node_base(n1) < __node_base(n2))) 
	{
		return -1;
	}
	else {
		return 1;
	}

}	

inline void node_add_to_next(node_t *cur, node_t *node)
{
	edge_list_add(cur->next, node);
}

inline void node_add_to_prev(node_t *cur, node_t *node)
{
	edge_list_add(cur->prev, node);
}

inline void node_remove_from_next(node_t *cur, node_t *node)
{
	if(0 == edge_list_remove(cur->next, node)) {
		srma_error(__func__, "Error removing node from the edge list", Exit, OutOfRange);
	}
}

inline void node_remove_from_prev(node_t *cur, node_t *node)
{
	if(0 == edge_list_remove(cur->prev, node)) {
		srma_error(__func__, "Error removing node from the edge list", Exit, OutOfRange);
	}
}

void node_free(node_t *node)
{
	int32_t i;
	for(i=0;i<node->next->length;i++) {
		node_remove_from_prev(node->next->nodes[i], node);
	}
	for(i=0;i<node->prev->length;i++) {
		node_remove_from_next(node->prev->nodes[i], node);
	}
	edge_list_free(node->next);
	edge_list_free(node->prev);
	free(node);
}

inline edge_list_t *edge_list_init()
{
	return srma_calloc(1, sizeof(edge_list_t), __func__, "edge_list_t");
}

void edge_list_free(edge_list_t *list)
{
	if(NULL == list) return;
	free(list->coverages);
	free(list->nodes);
	free(list);
}

// binary search
static inline int edge_list_get_i(edge_list_t *list, node_t *n)
{
	int32_t low, mid, high;

	if(0 == list->length) {
		return 0;
	}

	low = 0;
	high = list->length-1;
	mid = (low + high) / 2;
	while(low <= high) {
		mid = (low + high) / 2;
		int c = node_compare(list->nodes[mid], n);
		if(0 == c) {
			return mid;
		}
		else if(c < 0) {
			high = mid-1;
		}
		else {
			low = mid+1;
		}
	}
	// mid will be before, so move forward 
	while(mid < list->length && 0 < node_compare(list->nodes[mid], n)) {
		mid++;
	}
	return mid;
}

void edge_list_add(edge_list_t *list, node_t *node)
{
	int32_t i, j;

	i = edge_list_get_i(list, node);
	if(0 < list->length && i < list->length) {
		if(0 == node_compare(list->nodes[i], node)) { // same node, increment coverage
			list->coverages[i]++;
			return;
		}
	}

	list->length++;
	list->nodes = srma_realloc(list->nodes, sizeof(node_t*)*list->length, __func__, "list->nodes");
	list->coverages = srma_realloc(list->coverages, sizeof(uint16_t)*list->length, __func__, "list->coverages");
	for(j=list->length-1;i<j;j--) { // shift over
		list->nodes[j] = list->nodes[j-1];
		list->coverages[j] = list->coverages[j-1];
	}
	list->nodes[j] = node;
	list->coverages[j] = 1;
}

int32_t edge_list_remove(edge_list_t *list, node_t *node)
{
	int32_t i, j;
	if(0 == list->length) {
		return 0;
	}
	i = edge_list_get_i(list, node);
	list->coverages[i]--; // decrement coverage
	for(j=i+1;j<list->length;j++) {
		list->nodes[j-1] = list->nodes[j];
	}
	list->length--;
	if(0 == list->length) {
		free(list->nodes);
		free(list->coverages);
		list->nodes = NULL;
		list->coverages = NULL;
	}
	else {
		list->nodes = srma_realloc(list->nodes, sizeof(node_t*)*list->length, __func__, "list->nodes");
		list->coverages = srma_realloc(list->coverages, sizeof(uint16_t)*list->length, __func__, "list->coverages");
	}
	return 1;
}

inline node_list_t *node_list_init()
{
	return srma_calloc(1, sizeof(node_list_t), __func__, "node_list_t");
}

void node_list_free(node_list_t *list)
{
	int32_t i;
	if(NULL == list) return;
	for(i=0;i<list->length;i++) {
		node_free(list->nodes[i]);
	}
	free(list->nodes);
	free(list);
}

void node_list_clear(node_list_t *list)
{
	int32_t i;
	if(NULL == list) return;
	for(i=0;i<list->length;i++) {
		node_free(list->nodes[i]);
	}
	free(list->nodes);
	list->nodes = NULL;
	list->length = 0;
}

// binary search
static inline int node_list_get_i(node_list_t *list, node_t *n)
{
	int32_t low, mid, high;

	if(0 == list->length) {
		return 0;
	}

	low = 0;
	high = list->length-1;
	mid = (low + high) / 2;
	while(low <= high) {
		mid = (low + high) / 2;
		int c = node_compare(list->nodes[mid], n);
		if(0 == c) {
			return mid;
		}
		else if(c < 0) {
			high = mid-1;
		}
		else {
			low = mid+1;
		}
	}
	// mid will be before, so move forward 
	while(mid < list->length && 0 < node_compare(list->nodes[mid], n)) {
		mid++;
	}
	return mid;
}

void node_list_add(node_list_t *list, node_t *node)
{
	uint32_t i, j;

	i = node_list_get_i(list, node);

	if(0 < list->length && i < list->length && 0 == node_compare(list->nodes[i], node)) { // already exists
		return;
	}

	list->length++;
	list->nodes = srma_realloc(list->nodes, sizeof(node_t*)*list->length, __func__, "list->nodes");
	for(j=list->length-1;i<j;j--) { // shift over
		list->nodes[j] = list->nodes[j-1];
	}
	list->nodes[j] = node;
}
