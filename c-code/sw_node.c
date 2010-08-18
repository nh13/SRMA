#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_util.h"
#include "node.h"
#include "sw_node.h"

void sw_node_init(sw_node_t *sw_node_cur, sw_node_t *sw_node_prev, node_t *cur_node, int32_t coverage, uint8_t base, char qual, uint8_t use_qualities, uint8_t space)
{
	uint8_t b;

	if(NULL == sw_node_cur) {
		srma_error(__func__, "sw_node_cur", Exit, OutOfRange);
	}
	if(NULL == cur_node) {
		srma_error(__func__, "cur_node", Exit, OutOfRange);
	}

	sw_node_cur->node = cur_node;

	if(NULL == sw_node_prev) { // first
		sw_node_cur->read_offset = 0;
		sw_node_cur->score = 0;
		sw_node_cur->start_position = cur_node->position;
		sw_node_cur->prev_i = -1;
		sw_node_cur->coverage_sum = coverage;
		b = base;
	}
	else {
		sw_node_cur->read_offset = sw_node_prev->read_offset + 1;
		sw_node_cur->score = sw_node_prev->score;
		sw_node_cur->start_position = sw_node_prev->start_position;
		sw_node_cur->prev_i = sw_node_prev->cur_i;
		sw_node_cur->coverage_sum = sw_node_prev->coverage_sum + coverage;
		if(space == SRMA_SPACE_CS) {
			b = __node_base(sw_node_prev->node) ^ base; // color space
		}
		else {
			b = base;
		}
	}
	if(use_qualities == 1) {
		sw_node_cur->score += (b == __node_base(cur_node)) ? 0 : -1*srma_char2qual(qual);
	}
	else {
		sw_node_cur->score += (b == __node_base(cur_node)) ? 0 : -1;
	}
}

int sw_node_compare(sw_node_t *n1, sw_node_t *n2, int32_t type)
{
	// sort by:
	// - MIN/MAX genomic coordinate
	// - MIN read offset
	// - min node type
	// - min base
	// - min score

	// contig
	if(n1->node->contig < n2->node->contig) {
		return (SRMA_SW_HEAP_MIN == type) ? -1 : 1;
	}
	else if(n1->node->contig > n2->node->contig) {
		return (SRMA_SW_HEAP_MIN == type) ? 1 : -1;
	}

	// position
	if(n1->node->position < n2->node->position) {
		return (SRMA_SW_HEAP_MIN == type) ? -1 : 1;
	}
	else if(n1->node->position > n2->node->position) {
		return (SRMA_SW_HEAP_MIN == type) ? 1 : -1;
	}
	// read_offset
	if(n1->read_offset < n2->read_offset) {
		return -1;
	}
	else if(n1->read_offset > n2->read_offset) {
		return 1;
	}
	// type
	if(__node_type(n1->node) < __node_type(n2->node)) {
		return -1;
	}
	else if(__node_type(n1->node) > __node_type(n2->node)) {
		return 1;
	}
	// base
	if(__node_base(n1->node) < __node_base(n2->node)) {
		return -1;
	}
	else if(__node_base(n1->node) > __node_base(n2->node)) {
		return 1;
	}
	// score
	if(n1->score < n2->score) {
		return -1;
	}
	else if(n1->score > n2->score) {
		return 1;
	}
	// same
	return 0;

}

void sw_node_free(sw_node_t *node)
{
	free(node);
}
