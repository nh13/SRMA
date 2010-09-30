#include <stdlib.h>
#include <pthread.h>
#include <assert.h>

#include "samtools/sam.h"
#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_util.h"
#include "node.h"
#include "graph.h"

// See the following comments for synhronization
// --- SYNC ON --- 
// --- SYNC OFF --- 
static pthread_mutex_t graph_mutex = PTHREAD_MUTEX_INITIALIZER;

// DEBUG
void graph_check_nodes(graph_t *g)
{
	int32_t i, j;
	for(i=0;i<g->position_end-g->position_start+1;i++) {
		assert(NULL != g->nodes[i]);
		for(j=0;j<g->nodes[i]->length;j++) {
			assert(NULL != g->nodes[i]->nodes[j]);
			assert(g->position_start + i == g->nodes[i]->nodes[j]->position);
			assert(g->contig == g->nodes[i]->nodes[j]->contig);
		}
	}
	for(i=g->position_end-g->position_start+1;i<g->nodes_mem;i++) {
		assert(NULL != g->nodes[i]);
		assert(0 == g->nodes[i]->length);
		assert(NULL == g->nodes[i]->nodes);
	}
}

static inline void graph_nodes_realloc(graph_t *g, int32_t size)
{
	int32_t i, prev_mem;
	if(size < g->nodes_mem) {
		return;
	}
	prev_mem = g->nodes_mem;
	g->nodes_mem = size;
	roundup32(g->nodes_mem); // round up
	assert(size <= g->nodes_mem);
	g->nodes = srma_realloc(g->nodes, sizeof(node_list_t*)*g->nodes_mem, __func__, "g->nodes");
	g->coverages = srma_realloc(g->coverages, sizeof(uint16_t)*g->nodes_mem, __func__, "g->coverages");
	for(i=prev_mem;i<g->nodes_mem;i++) {
		g->nodes[i] = node_list_init();
		g->coverages[i] = 0;
	}
}

graph_t *graph_init()
{
	graph_t *g = NULL;

	g = srma_malloc(sizeof(graph_t), __func__, "g");
	g->contig = 1;
	g->position_start = 1;
	g->position_end = 1;
	g->nodes = NULL;
	g->coverages = NULL;
	g->is_empty = 1;

	// Add a dummy node
	g->nodes = srma_malloc(sizeof(node_list_t*), __func__, "g->nodes");
	g->nodes[0] = node_list_init();
	g->coverages = srma_calloc(1, sizeof(uint16_t), __func__, "g->coverages");
	g->coverages[0] = 0;
	g->nodes_mem = 1;

	return g;
}

node_t *graph_add_sam(graph_t *g, bam1_t *b, ref_t *ref, int32_t use_threads)
{
	bam_aln_t *aln = NULL;
	int32_t aln_start, aln_index, ref_index, aln_ref_index;
	int32_t i;
	node_t *prev_node=NULL, *cur_node=NULL, *ret_node=NULL;
	uint8_t type, strand;

	aln_start = b->core.pos+1;
	aln_ref_index = b->core.tid;
	aln = bam_aln_init(b, ref);
	strand =  bam1_strand(b);

	// --- SYNC ON --- 
	if(1 == use_threads) pthread_mutex_lock(&graph_mutex); // synchronize start
	if(aln_start < g->position_start) {
		int32_t diff = g->position_start - aln_start;
		graph_nodes_realloc(g, g->position_end - aln_start + 1); // alloc more memory if needed
		// shift up
		for(i=g->position_end-g->position_start;0<=i;i--) {
			node_list_t *list = g->nodes[i+diff];
			uint16_t coverage = g->coverages[i+diff];
			// swap
			g->nodes[i+diff] = g->nodes[i];
			g->coverages[i+diff] = g->coverages[i];
			g->nodes[i] = list;
			g->coverages[i] = coverage;
		}
		g->position_start = aln_start;
	}

	if(1 == g->is_empty) {
		for(i=0;i<g->position_end - g->position_start + 1;i++) {
			node_list_clear(g->nodes[i]);
			assert(0 == g->nodes[i]->length); // DEBUG
			g->coverages[i] = 0;
		}
		g->position_start = aln_start;
		if(ALN_GAP == aln->ref[0]) {
			g->position_start--;
		}
		g->position_end = g->position_start;
		g->contig = aln_ref_index + 1;
		g->is_empty = 0;
	}
	if(1 == use_threads) pthread_mutex_unlock(&graph_mutex); // synchronize end 
	// --- SYNC OFF --- 
	for(aln_index=0,ref_index=-1;aln_index<aln->length;aln_index++,prev_node=cur_node) {

		// Skip over a deletion
		while(ALN_GAP == aln->read[aln_index]) {
			aln_index++;
			ref_index++;
		}

		if(aln->read[aln_index] == aln->ref[aln_index]) { // match
			type = NODE_MATCH;
		}
		else if(aln->ref[aln_index] == ALN_GAP) { // insertion
			type = NODE_INSERTION;
		}
		else { // mismatch
			type = NODE_MISMATCH;
		}
		if(NULL == prev_node || NODE_INSERTION != __node_type(prev_node)) { // previous was an insertion, already on the position
			ref_index++;
		}

		cur_node = graph_add_node(g, 
				node_init(aln->read[aln_index], type, g->contig, aln_start + ref_index, prev_node),
				prev_node,
				use_threads);
	
		if(NULL == prev_node && 0 == strand) { // first node and forward strand
			ret_node = cur_node;
		}
	}

	if(1 == strand) {
		ret_node = cur_node;
	}

	bam_aln_free(aln);
	return ret_node;
}

node_t *graph_add_node(graph_t *g, node_t *node, node_t *prev, int32_t use_threads)
{
	node_t *cur = NULL;

	// --- SYNC ON --- 
	if(1 == use_threads) pthread_mutex_lock(&graph_mutex); // synchronize start

	cur = graph_contains(g, node);
	if(NULL == cur) {
		if(node->contig != g->contig) {
			srma_error(__func__, "Control reached unexpected point", Exit, OutOfRange);
		}
		graph_nodes_realloc(g, node->position - g->position_start + 1); // alloc more memory if needed
		assert(node->position - g->position_start < g->nodes_mem); // DEBUG
		assert(NULL != g->nodes[node->position - g->position_start]); // DEBUG
		node_list_add(g->nodes[node->position - g->position_start], node);
		if(NODE_INSERTION != __node_type(node) || 0 == node->offset) {
			g->coverages[node->position - g->position_start] += node->coverage;
		}
		if(g->position_end < node->position) {
			g->position_end = node->position;
		}
		cur = node;
		g->is_empty = 0;
	}
	else {
		node_free(node);
		node=NULL;
		cur->coverage++;
		if(NODE_INSERTION != __node_type(cur) || 0 == cur->offset) {
			g->coverages[cur->position - g->position_start]++;
		}
	}

	if(NULL != prev) {
		node_add_to_prev(cur, prev);
		node_add_to_next(prev, cur);
		assert(0 < cur->prev->length);
		assert(0 < prev->next->length);
	}

	if(1 == use_threads) pthread_mutex_unlock(&graph_mutex); // synchronize end 
	// --- SYNC OFF --- 

	return cur;
}

node_t *graph_contains(graph_t *g, node_t *node)
{
	uint32_t i;
	node_list_t *list = NULL;

	list = graph_get_node_list(g, node->position);
	if(NULL == list) return NULL;

	for(i=0;i<list->length;i++) {
		assert(node->contig == list->nodes[i]->contig); // DEBUG
		assert(node->position == list->nodes[i]->position); // DEBUG
		if(0 == node_compare2(node, list->nodes[i])) {
			return list->nodes[i];
		}
	}

	return NULL;
}

uint32_t graph_get_node_list_index_at_or_after(graph_t *g, uint32_t position)
{
	if(position < g->position_start) {
		return g->position_start;
	}

	while(position <= g->position_end) {
		if(g->nodes[position - g->position_start]->length > 0) {
			return position;
		}
		position++;
	}
	return 0;
}

uint32_t graph_get_node_list_index_at_or_before(graph_t *g, uint32_t position)
{
	if(g->position_end < position) {
		return g->position_end;
	}

	while(g->position_start <= position) {
		if(g->nodes[position - g->position_start]->length > 0) {
			return position;
		}
		position--;
	}
	return 0;
}

node_list_t *graph_get_node_list(graph_t *g, uint32_t position)
{
	if(position < g->position_start || g->position_end < position) {
		return NULL;
	}
	else {
		return g->nodes[position - g->position_start];
	}
}

uint16_t graph_get_coverage(graph_t *g, uint32_t position)
{
	if(position < g->position_start || g->position_end < position) {
		return 0;
	}
	else {
		return g->coverages[position - g->position_start];
	}
}

// Not thread-safe
void graph_prune(graph_t *g, uint32_t contig_index, uint32_t alignment_start, int32_t offset)
{
	uint32_t i;
	uint8_t clear = 0;

	if(g->contig != 1 + contig_index) {
		clear = 1;
	}
	else if(g->position_start < alignment_start - offset) { // unreachable start node
		if(g->position_end < alignment_start - offset) { // all nodes unreachable
			clear = 1;
		}
		else {
			// free and shift down
			uint32_t diff = alignment_start - offset - g->position_start;
			for(i=0;i<g->position_end-(alignment_start-offset)+1;i++) {
				node_list_t *ptr = g->nodes[i];
				node_list_clear(ptr);
				assert(0 == ptr->length); // DEBUG
				// swap
				assert(i+diff <= g->position_end);
				g->nodes[i] = g->nodes[i+diff];
				g->nodes[i+diff] = ptr;
				g->coverages[i] = g->coverages[i+diff];
				g->coverages[i+diff] = 0;
			}
			for(;i<g->position_end-g->position_start+1;i++) { // clear the rest
				node_list_clear(g->nodes[i]);
				assert(0 == g->nodes[i]->length); // DEBUG
			}
			g->position_start = alignment_start - offset;
		}
	}
	if(1 == clear) {
		for(i=0;i<g->position_end - g->position_start + 1;i++) {
			node_list_clear(g->nodes[i]);
			assert(0 == g->nodes[i]->length); // DEBUG
			g->coverages[i] = 0;
		}
		g->contig = contig_index + 1;
		g->position_start = alignment_start;
		g->position_end = alignment_start;

		g->is_empty = 1;
	}
}

void graph_free(graph_t *g)
{
	uint32_t i;
	for(i=0;i<g->nodes_mem;i++) {
		node_list_free(g->nodes[i]);
	}
	free(g->nodes);
	free(g->coverages);
	free(g);
}
