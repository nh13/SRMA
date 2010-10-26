#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "samtools/bam.h"
#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_util.h"
#include "node.h"
#include "graph.h"
#include "sw_node.h"
#include "sw_heap.h"
#include "sw_align.h"

// 1 for A, 2 for C, 4 for G,
// 8 for T and 15 for N
int nt4bit_to_int[16] = {
	-1, 0, 1, -1,
	2, -1, -1, -1, 
	3, -1, -1, -1,
	-1, -1, -1, 4
};

int int_to_nt4bit[5] = {1, 2, 4, 8, 15};

// DEBUG FUNCTION
/*
   static void print_flag(bam1_t *b)
   {
   static char *bam_flag2char_table = "pPuUrR12sfd\0\0\0\0\0";
   int32_t i;
   for (i = 0; i < 16; ++i)
   if ((b->core.flag & 1<<i) && bam_flag2char_table[i])
   fputc(bam_flag2char_table[i], stderr);
   fputc('\n', stderr);
   }
   */

static char *sw_align_get_cs(bam1_t *b)
{
	uint8_t *c = bam_aux_get(b, "CS");
	// return the base if the tag was not found
	if(0 == c) return NULL;
	return bam_aux2Z(c);
}

static char *sw_align_get_cq(bam1_t *b)
{
	uint8_t *c = bam_aux_get(b, "CQ");
	// return the base if the tag was not found
	if(0 == c) return NULL;
	return bam_aux2Z(c);
}

static int pass_filters(graph_t *graph, node_t *node, int32_t to_node_coverage, cov_cutoffs_t *cutoffs, int32_t max_total_coverage)
{
	int32_t total_coverage;

	total_coverage = graph_get_coverage(graph, node->position);
	// DEBUG
	/*
	   fprintf(stderr, "max_total_coverage=%d total_coverage=%d",
	   max_total_coverage, total_coverage);
	   fprintf(stderr, " cov_cutoffs_get(cutoffs, total_coverage)=%d to_node_coverage=%d\n",
	   cov_cutoffs_get(cutoffs, total_coverage), to_node_coverage);
	   */
	if(max_total_coverage < total_coverage) {
		return -1;
	}
	else if(cov_cutoffs_get(cutoffs, total_coverage) <= to_node_coverage) {
		return 0;
	}
	else {
		return 1;
	}
}

static int pass_filters1(graph_t *graph, node_t *node, cov_cutoffs_t *cutoffs, int32_t max_total_coverage)
{
	return pass_filters(graph, node, node->coverage, cutoffs, max_total_coverage);
}

static int32_t sw_align_bound(graph_t *g, bam1_t *b, node_t *n, sw_heap_t *heap, uint8_t strand, const char *colors, const char *color_qualities, uint8_t space, cov_cutoffs_t *cutoffs, uint8_t use_qualities, int32_t max_total_coverage, int32_t max_heap_size)
{
	int32_t sw_node_i=-1, sw_node_best_i=-1, sw_node_cur_i=-1, sw_node_next_i=-1;
	int32_t i;
	char base, qual;

	if(0 != pass_filters1(g, n, cutoffs, max_total_coverage)) {
		//fprintf(stderr, "NOT BOUNDED 1\n"); // DEBUG
		return -1;
	}

	{ // add the start node to the heap
		// Get first base
		if(SRMA_SPACE_CS == space) {
			base = nt2int_table[(int)colors[1]];
			qual = nt2int_table[(int)color_qualities[0]]; 
		}
		else {
			if(strand) {
				base = nt4bit_to_int[bam1_seqi(bam1_seq(b), b->core.l_qseq-1)];
				qual = bam1_qual(b)[b->core.l_qseq-1];
			}
			else {
				base = nt4bit_to_int[bam1_seqi(bam1_seq(b), 0)];
				qual = bam1_qual(b)[0];
			}
		}
		sw_node_i = sw_heap_get_node_i(heap);
		sw_node_init(&heap->nodes[sw_node_i], NULL, n, n->coverage, base, qual, use_qualities, space); 
		sw_heap_add_i(heap, sw_node_i);
	}

	sw_node_cur_i = sw_heap_poll_i(heap);
	assert(0 <= sw_node_cur_i); // DEBUG
	while(0 <= sw_node_cur_i) {
		//fprintf(stderr, "sw_node_cur_i=%d\n", sw_node_cur_i); // DEBUG
		if(max_heap_size <- heap->queue_end - heap->queue_start + 1) {
			// too many to consider
			//fprintf(stderr, "NOT BOUNDED 2\n"); // DEBUG
			return -1;
		}

		sw_node_next_i = sw_heap_peek_i(heap);
		assert(0 <= sw_node_cur_i); // DEBUG
		while(NODE_INSERTION != __node_type(heap->nodes[sw_node_cur_i].node)
				&& 0 <= sw_node_next_i
				&& 0 == sw_node_compare(&heap->nodes[sw_node_cur_i], &heap->nodes[sw_node_next_i], heap->type)) {
			if(heap->nodes[sw_node_cur_i].score < heap->nodes[sw_node_next_i].score ||
					(heap->nodes[sw_node_cur_i].score == heap->nodes[sw_node_next_i].score &&
					 heap->nodes[sw_node_cur_i].coverage_sum < heap->nodes[sw_node_next_i].coverage_sum)) { 
				sw_node_cur_i = sw_heap_poll_i(heap);
			}
			else {
				// ignore the next node
				sw_heap_poll_i(heap);
			}
			sw_node_next_i = sw_heap_peek_i(heap);
		}
		sw_node_next_i = -1;

		// DEBUG
		/*
		   fprintf(stderr, "read_offset=%d l_qseq-1=%d\n",
		   heap->nodes[sw_node_cur_i].read_offset,
		   b->core.l_qseq-1);
		   */
		if(heap->nodes[sw_node_cur_i].read_offset == b->core.l_qseq-1) { // found, keep best
			if(sw_node_best_i < 0 ||
					heap->nodes[sw_node_best_i].score < heap->nodes[sw_node_cur_i].score ||
					(heap->nodes[sw_node_best_i].score == heap->nodes[sw_node_cur_i].score && 
					 heap->nodes[sw_node_best_i].coverage_sum < heap->nodes[sw_node_cur_i].coverage_sum)) {
				sw_node_best_i = sw_node_cur_i;
			}
		}
		else {
			edge_list_t *list = NULL;
			if(1 == strand) { // reverse
				list = heap->nodes[sw_node_cur_i].node->prev;
			}
			else {
				list = heap->nodes[sw_node_cur_i].node->next;
			}
			{ // get the aligned base and quality
				// do not use color space data for bounding
				if(strand) {
					base = nt4bit_to_int[bam1_seqi(bam1_seq(b), b->core.l_qseq-1-heap->nodes[sw_node_cur_i].read_offset-1)];
					qual = bam1_qual(b)[b->core.l_qseq-1-heap->nodes[sw_node_cur_i].read_offset-1];
				}
				else {
					base = nt4bit_to_int[bam1_seqi(bam1_seq(b), (heap->nodes[sw_node_cur_i].read_offset+1))];
					qual = bam1_qual(b)[(heap->nodes[sw_node_cur_i].read_offset+1)];
				}
			}
			//fprintf(stderr, "list->length=%d\n", list->length); // DEBUG
			for(i=0;i<list->length;i++) {
				node_t *node_cur= list->nodes[i];
				// DEBUG
				/*
				   fprintf(stderr, "%d:%d __node_base(node_cur)=%d base=%d __node_type(node_cur)=%d coverages_cur=%d\n",
				   node_cur->contig, node_cur->position,
				   __node_base(node_cur), base, __node_type(node_cur), list->coverages[i]);
				   */
				// base should match unless filters don't pass
				if(__node_base(node_cur) == base) {
					uint16_t coverage_cur = list->coverages[i];
					int32_t pass = pass_filters(g, node_cur, coverage_cur, cutoffs, max_total_coverage);
					//fprintf(stderr, "pass=%d\n", pass); // DEBUG
					if(0 == pass) {
						if(SRMA_SPACE_CS == space) { // use color space data
							base = nt2int_table[(int)colors[1 + (heap->nodes[sw_node_cur_i].read_offset+1)]];
							qual = nt2int_table[(int)color_qualities[heap->nodes[sw_node_cur_i].read_offset+1]]; 
						}
						// add to the heap
						sw_node_i = sw_heap_get_node_i(heap);
						// DEBUG
						assert(0 <= sw_node_cur_i);
						assert(0 <= heap->nodes[sw_node_cur_i].read_offset);
						sw_node_init(&heap->nodes[sw_node_i], &heap->nodes[sw_node_cur_i], node_cur, coverage_cur, base, qual, use_qualities, space); 
						sw_heap_add_i(heap, sw_node_i);
					}
					else if(pass < 0) {
						//fprintf(stderr, "NOT BOUNDED 3\n"); // DEBUG
						return -1;
					}
				}
			}
		}
		// get the next node
		sw_node_cur_i = sw_heap_poll_i(heap);
	}

	//fprintf(stderr, "BOUNDED %d\n", sw_node_best_i); // DEBUG
	return sw_node_best_i;
}

static int32_t sw_align_get_soft_clip(bam1_t *b, int32_t is_end)
{
	int32_t n;
	n = bam1_cigar(b)[(0 == is_end) ? 0 : b->core.n_cigar-1];
	if(BAM_CSOFT_CLIP == (n & BAM_CIGAR_MASK)) { // soft-clipping
		return (n >> BAM_CIGAR_SHIFT);
	}
	return 0;
}

// TODO soft-clipping
bam1_t *sw_align(graph_t *g, bam1_t *b, node_t *n, sw_heap_t *heap, char *rg_id, int32_t offset, cov_cutoffs_t *cutoffs, uint8_t correct_bases, uint8_t use_qualities, int32_t max_total_coverage, int32_t max_heap_size)
{
	char *colors = NULL;
	char *color_qualities = NULL;
	char base, qual;
	uint8_t space = SRMA_SPACE_NT;
	uint8_t strand;
	int32_t i, j, aln_start;
	int32_t num_start_nodes_added=0;
	int32_t sw_node_i=-1, sw_node_best_i=-1, sw_node_cur_i=-1, sw_node_next_i=-1;
	int32_t soft_clip_start_l = 0, soft_clip_end_l = 0;

	strand = bam1_strand(b);

	// soft-clipping
	if(1 == strand) { //reverse
		// going from 3'->5'
		soft_clip_start_l = sw_align_get_soft_clip(b, 1); 
		soft_clip_end_l = sw_align_get_soft_clip(b, 0);
	}
	else {
		// going from 5'->3'
		soft_clip_start_l = sw_align_get_soft_clip(b, 0); 
		soft_clip_end_l = sw_align_get_soft_clip(b, 1);
	}
	// FOR NOW
	if(0 < soft_clip_start_l || 0 < soft_clip_end_l) {
		return b;
	}

	// Check color space
	colors = sw_align_get_cs(b);
	if(NULL == colors) {
		space = SRMA_SPACE_NT;
	}
	else {
		space = SRMA_SPACE_CS;
		color_qualities  = sw_align_get_cq(b);
		// Some aligners include a quality value for the adapter.  A quality value
		// IMHO should not be given for an unobserved (assumed) peice of data.  Trim
		// the first quality in this case
		if(strlen(colors) == strlen(color_qualities)) {  // ignore leading quality
			color_qualities++;
		}
		if(0 < soft_clip_start_l || 0 < soft_clip_end_l) {
			srma_error(__func__, "Soft clipping not supported for color space", Exit, OutOfRange);
		}
	}	

	// remove mate info 
	b->core.flag &= ~(BAM_FPROPER_PAIR | BAM_FMREVERSE | BAM_FMUNMAP);
	b->core.mtid = -1;
	b->core.mpos = -1;
	b->core.isize = 0;

	// re-type heap
	heap->type = (1 == strand) ? SRMA_SW_HEAP_MAX : SRMA_SW_HEAP_MIN;

	// bound with original alignment
	sw_node_best_i = sw_align_bound(g, b, n, heap, strand, colors, color_qualities, space, cutoffs, use_qualities, max_total_coverage, max_heap_size);
	if(0 <= sw_node_best_i) {
		sw_heap_reset(heap); // reset the heap, keep old nodes
		/*
		   fprintf(stderr, "BOUNDED score=%d coverage_sum=%hu\n", 
		   heap->nodes[sw_node_best_i].score,
		   heap->nodes[sw_node_best_i].coverage_sum); // DEBUG
		   */
	}
	else {
		//fprintf(stderr, "NOT BOUNDED\n"); // DEBUG
		// nodes do not need to be preserved
		sw_heap_clear(heap);
	}
	//return b; // HERE DEBUG HERE BUG

	// add start nodes
	if(strand) {
		if(SRMA_SPACE_CS == space) {
			base = nt2int_table[(int)colors[1]];
			qual = nt2int_table[(int)color_qualities[0]]; 
		}
		else {
			base = nt4bit_to_int[bam1_seqi(bam1_seq(b), b->core.l_qseq-1)];
			qual = bam1_qual(b)[b->core.l_qseq-1];
		}
		aln_start = bam_calend(&b->core, bam1_cigar(b));
		for(i=aln_start+offset;aln_start-offset<=i;i--) {
			int32_t pos = graph_get_node_list_index_at_or_before(g, i);
			node_list_t *list = graph_get_node_list(g, pos);
			if(0 != pos && NULL != list) {
				for(j=0;j<list->length;j++) {
					node_t *node = list->nodes[j];
					int32_t pass = pass_filters1(g, node, cutoffs, max_total_coverage);
					if(0 == pass) {
						sw_node_i = sw_heap_get_node_i(heap);
						sw_node_init(&heap->nodes[sw_node_i], NULL, node, node->coverage, base, qual, use_qualities, space); 
						sw_heap_add_i(heap, sw_node_i);
					}
					else if(pass < 0) {
						sw_heap_clear(heap); // clear heap
						return b;
					}
					if(node->position < i) {
						i = node->position;
					}
					num_start_nodes_added++;
				}
			}
		}
	}
	else {
		if(SRMA_SPACE_CS == space) {
			base = nt2int_table[(int)colors[1]];
			qual = nt2int_table[(int)color_qualities[0]]; 
		}
		else {
			base = nt4bit_to_int[bam1_seqi(bam1_seq(b), 0)];
			qual = bam1_qual(b)[0];
		}
		aln_start = b->core.pos;
		for(i=aln_start-offset;i<=aln_start+offset;i++) {
			int32_t pos = graph_get_node_list_index_at_or_after(g, i);
			node_list_t *list = graph_get_node_list(g, pos);
			if(0 != pos && NULL != list) {
				for(j=0;j<list->length;j++) {
					node_t *node = list->nodes[j];
					int32_t pass = pass_filters1(g, node, cutoffs, max_total_coverage);
					if(0 == pass) {
						sw_node_i = sw_heap_get_node_i(heap);
						sw_node_init(&heap->nodes[sw_node_i], NULL, node, node->coverage, base, qual, use_qualities, space); 
						sw_heap_add_i(heap, sw_node_i);
					}
					else if(pass < 0) {
						sw_heap_clear(heap); // clear heap
						return b;
					}
					if(node->position < i) {
						i = node->position;
					}
					num_start_nodes_added++;
				}
			}
		}
	}
	if(0 == num_start_nodes_added) {
		srma_error(__func__, "Did not add any start nodes", Exit, OutOfRange);
	}

	sw_node_cur_i = sw_heap_poll_i(heap);
	while(0 <= sw_node_cur_i) {
		if(max_heap_size <- heap->queue_end - heap->queue_start + 1) {
			// too many to consider
			sw_heap_clear(heap); // clear heap
			return b;
		}

		sw_node_next_i = sw_heap_peek_i(heap);
		assert(0 <= sw_node_cur_i); // DEBUG
		while(NODE_INSERTION != __node_type(heap->nodes[sw_node_cur_i].node)
				&& 0 <= sw_node_next_i
				&& 0 == sw_node_compare(&heap->nodes[sw_node_cur_i], &heap->nodes[sw_node_next_i], heap->type)) {
			if(heap->nodes[sw_node_cur_i].score < heap->nodes[sw_node_next_i].score ||
					(heap->nodes[sw_node_cur_i].score == heap->nodes[sw_node_next_i].score &&
					 heap->nodes[sw_node_cur_i].coverage_sum < heap->nodes[sw_node_next_i].coverage_sum)) { 
				sw_node_cur_i = sw_heap_poll_i(heap);
			}
			else {
				// ignore the next node
				sw_heap_poll_i(heap);
			}
			sw_node_next_i = sw_heap_peek_i(heap);
		}
		sw_node_next_i = -1;

		if(heap->nodes[sw_node_cur_i].read_offset == b->core.l_qseq-1) { // found, keep best
			if(sw_node_best_i < 0 ||
					heap->nodes[sw_node_best_i].score < heap->nodes[sw_node_cur_i].score ||
					(heap->nodes[sw_node_best_i].score == heap->nodes[sw_node_cur_i].score && 
					 heap->nodes[sw_node_best_i].coverage_sum < heap->nodes[sw_node_cur_i].coverage_sum)) {
				//fprintf(stderr, "FOUND BEST\n"); // DEBUG
				sw_node_best_i = sw_node_cur_i;
			}
		}
                else if(0 <= sw_node_best_i && 
                        heap->nodes[sw_node_cur_i].score < heap->nodes[sw_node_best_i].score) {
                        // ignore, under the assumption that scores can only
                        // become more negative.
                }
		else {
			edge_list_t *list = NULL;
			if(1 == strand) { // reverse
				list = heap->nodes[sw_node_cur_i].node->prev;
			}
			else {
				list = heap->nodes[sw_node_cur_i].node->next;
			}
			{ // get the base and quality
				if(SRMA_SPACE_CS == space) {
					base = nt2int_table[(int)colors[1 + (heap->nodes[sw_node_cur_i].read_offset+1)]];
					qual = nt2int_table[(int)color_qualities[heap->nodes[sw_node_cur_i].read_offset+1]]; 
				}
				else {
					if(strand) {
						base = nt4bit_to_int[bam1_seqi(bam1_seq(b), b->core.l_qseq-1-heap->nodes[sw_node_cur_i].read_offset-1)];
						qual = bam1_qual(b)[b->core.l_qseq-1-heap->nodes[sw_node_cur_i].read_offset-1];
					}
					else {
						base = nt4bit_to_int[bam1_seqi(bam1_seq(b), (heap->nodes[sw_node_cur_i].read_offset+1))];
						qual = bam1_qual(b)[(heap->nodes[sw_node_cur_i].read_offset+1)];
					}
				}
			}
			/*
			   node_t *node = heap->nodes[sw_node_cur_i].node;
			   fprintf(stderr, "NODE %d:%d offset=%d coverage=%d base=%d\n",
			   node->contig, node->position, node->offset, node->coverage, node->base);
			   fprintf(stderr, "SW_NODE read_offset=%d score=%d coverage_sum=%d start_position=%d space=%d\n",
			   heap->nodes[sw_node_cur_i].read_offset, heap->nodes[sw_node_cur_i].score, heap->nodes[sw_node_cur_i].coverage_sum, heap->nodes[sw_node_cur_i].start_position, space);
			   */
			for(i=0;i<list->length;i++) {
				node_t *node_cur = list->nodes[i];
				uint16_t coverage_cur = list->coverages[i];
				int32_t pass = pass_filters(g, node_cur, coverage_cur, cutoffs, max_total_coverage);
				if(0 == pass) {
					// add to the heap
					sw_node_i = sw_heap_get_node_i(heap);
					// DEBUG
					assert(0 <= sw_node_cur_i);
					assert(0 <= heap->nodes[sw_node_cur_i].read_offset);
					sw_node_init(&heap->nodes[sw_node_i], &heap->nodes[sw_node_cur_i], node_cur, coverage_cur, base, qual, use_qualities, space); 
					sw_heap_add_i(heap, sw_node_i);
				}
				else if(pass < 0) {
					sw_heap_clear(heap); // clear heap
					return b;
				}
			}
		}
		// get the next node
		sw_node_cur_i = sw_heap_poll_i(heap);
	}

	/*
	// fprintf(stderr, "sw_node_best_i=%d\n", sw_node_best_i); // DEBUG
	if(0 <= sw_node_best_i) {
	fprintf(stderr, "END score=%d coverage_sum=%hu\n", 
	heap->nodes[sw_node_best_i].score,
	heap->nodes[sw_node_best_i].coverage_sum); // DEBUG
	}
	*/
	// update SAM/BAM
	b = sw_align_update_bam(b, rg_id, heap, sw_node_best_i, space, colors, color_qualities, strand, correct_bases);
	sw_heap_clear(heap); // clear heap
	return b;
}	

static inline void sw_align_bam_alloc_data(bam1_t *bam, int size)
{
	if (bam->m_data < size) {
		bam->m_data = size;
		roundup32(bam->m_data); // from bam.h
		bam->data = (uint8_t*)realloc(bam->data, bam->m_data);
	}
}

// Bound quality score
// For: sw_align_update_bam
#define __bound_qual(_qual) ((_qual <= 0) ? 1 : ((93 < _qual) ? 93 : _qual))

// Copy old tag into the new bam
// For: sw_align_update_bam
#define __copy_old(_tag) do { \
	s = bam_aux_get(bam_old, _tag); \
	if(NULL != s) { \
		bam_aux_append(bam_new, _tag, *s, 1+strlen((char*)(s+1)), s+1); \
	} \
} while(0)

// TODO soft clipping
bam1_t *sw_align_update_bam(bam1_t *bam_old, char *rg_id, sw_heap_t *heap, int32_t sw_node_best_i, uint8_t space, char *colors, char *color_qualities, uint8_t strand, uint8_t correct_bases)
{
	bam1_t *bam_new=NULL;
	int32_t sw_node_cur_i=-1, sw_node_prev_i=-1;
	int32_t i;
	int32_t cigar_cur_op, cigar_prev_op;
	int32_t cigar_cur_length, cigar_prev_length;
	uint32_t read_index;

	if(sw_node_best_i < 0) { // none found, do not modify alignment
		return bam_old;
	}

	bam_new = srma_calloc(1, sizeof(bam1_t), __func__, "bam_new");

	if(1 == strand) {
		read_index = 0;
	}
	else {
		read_index = bam_old->core.l_qseq-1;
	}

	{ // query name
		bam_new->core.l_qname = bam_old->core.l_qname;
		bam_new->data_len += bam_new->core.l_qname;
		sw_align_bam_alloc_data(bam_new, bam_new->data_len);
		memcpy(bam1_qname(bam_new), bam1_qname(bam_old), bam_old->core.l_qname);
	}
	{ // flag
		bam_new->core.flag = bam_old->core.flag;
	}
	{ // tid, pos, qual
		bam_new->core.tid = heap->nodes[sw_node_best_i].node->contig-1; // it is one-based, we want zero-based
		if(1 == strand) { // reverse strand
			bam_new->core.pos = heap->nodes[sw_node_best_i].node->position-1;
		}
		else {
			bam_new->core.pos = heap->nodes[sw_node_best_i].start_position-1; // zero-based
		}
		bam_new->core.qual = bam_old->core.qual; // should we change the mapping quality?
		bam_new->core.mtid = -1;
		bam_new->core.mpos = -1;
		bam_new->core.isize = 0;
	}
	{ // cigar length
		bam_new->core.n_cigar = 0;
		cigar_cur_op = cigar_prev_op = -1;
		sw_node_cur_i = sw_node_best_i;
		while(0 <= sw_node_cur_i) {
			if(0 <= sw_node_prev_i  && BAM_CDEL == cigar_prev_op && 1 < fabs(heap->nodes[sw_node_cur_i].node->position - heap->nodes[sw_node_prev_i].node->position)) {
				cigar_cur_op = BAM_CDEL;
			}	
			else {
				switch(__node_type(heap->nodes[sw_node_cur_i].node)) {
					case NODE_MATCH:
					case NODE_MISMATCH:
						cigar_cur_op = BAM_CMATCH;
						break;
					case NODE_INSERTION:
						cigar_cur_op = BAM_CINS;
						break;
					default:
						srma_error(__func__, "unknown node type", Exit, OutOfRange);
				}
			}
			if(cigar_prev_op != cigar_cur_op) {
				// update the previous cigar operator
				cigar_prev_op = cigar_cur_op;
				bam_new->core.n_cigar++;
			}
			// Update
			if(BAM_CDEL != cigar_cur_op) {
				sw_node_prev_i = sw_node_cur_i;
				sw_node_cur_i = heap->nodes[sw_node_cur_i].prev_i;
			}
		}
	}

	{ // cigar and seq
		uint32_t *cigar_ptr=NULL;
		uint8_t *seq_ptr=NULL;
		uint32_t cigar_i = 0;
		// cigar
		bam_new->data_len += bam_new->core.n_cigar*sizeof(uint32_t);
		sw_align_bam_alloc_data(bam_new, bam_new->data_len);
		cigar_ptr = bam1_cigar(bam_new);
		// seq
		bam_new->core.l_qseq = bam_old->core.l_qseq;
		bam_new->data_len += (bam_new->core.l_qseq + 1)/2;
		sw_align_bam_alloc_data(bam_new, bam_new->data_len);
		seq_ptr = bam1_seq(bam_new);
		// fill in cigar and seq
		cigar_i = (1 == strand) ? bam_new->core.n_cigar-1 : 0;
		cigar_cur_op = cigar_prev_op = -1;
		cigar_cur_length = cigar_prev_length = -1;
		sw_node_cur_i = sw_node_best_i;
		while(0 <= sw_node_cur_i) {
			if(0 <= sw_node_prev_i && BAM_CDEL == cigar_prev_op && 1 < fabs(heap->nodes[sw_node_cur_i].node->position - heap->nodes[sw_node_prev_i].node->position)) {
				cigar_cur_op = BAM_CDEL;
			}	
			else {
				switch(__node_type(heap->nodes[sw_node_cur_i].node)) {
					case NODE_MATCH:
					case NODE_MISMATCH:
						cigar_cur_op = BAM_CMATCH;
						break;
					case NODE_INSERTION:
						cigar_cur_op = BAM_CINS;
						break;
					default:
						srma_error(__func__, "unknown node type", Exit, OutOfRange);
				}
				// pack sequence
				if(1 == strand && 0 == read_index%2) {
					seq_ptr[read_index/2] = 0;
				}
				else if(0 == strand && 1 == read_index%2) {
					seq_ptr[read_index/2] = 0;
				}
				// DEBUG
				/*
				   fprintf(stderr, "read_index=%d base=%d\n",
				   read_index, __node_base(heap->nodes[sw_node_cur_i].node));
				   */
				seq_ptr[read_index/2] |= int_to_nt4bit[__node_base(heap->nodes[sw_node_cur_i].node)] << 4*(1-(read_index%2));
				if(1 == strand) {
					read_index++;
				}
				else {
					read_index--;
				}
			}
			if(cigar_prev_op != cigar_cur_op) {
				// add the previous cigar operator
				if(-1 != cigar_prev_op) {
					bam1_cigar(bam_new)[cigar_i] = (cigar_prev_length << BAM_CIGAR_SHIFT) | cigar_prev_op;
					if(1 == strand) { // reverse strand
						cigar_i--;
					}
					else {
						cigar_i++;
					}
				}

				// update the previous cigar operator
				cigar_prev_op = cigar_cur_op;
				if(cigar_cur_op == BAM_CDEL) {
					// deletion length
					cigar_prev_length = (int)fabs(heap->nodes[sw_node_cur_i].node->position - heap->nodes[sw_node_cur_i].node->position) - 1;
				}
				else {
					cigar_prev_length = 1;
				}
			}
			else {
				cigar_prev_length++;
			}
			// Update
			if(BAM_CDEL != cigar_cur_op) {
				sw_node_prev_i = sw_node_cur_i;
				sw_node_cur_i = heap->nodes[sw_node_cur_i].prev_i;
			}
		}
		if(0 < cigar_prev_length) {
			if(-1 == cigar_prev_op || BAM_CDEL == cigar_prev_op) {
				srma_error(__func__, "Alignment ended with a null cigar or a deletion", Exit, OutOfRange);
			}	
			bam1_cigar(bam_new)[cigar_i] = (cigar_prev_length << BAM_CIGAR_SHIFT) | cigar_prev_op;
			// DEBUG
			if(1 == strand) { // reverse strand
				assert(cigar_i == 0);
			}
			else {
				assert(cigar_i == bam_new->core.n_cigar-1);
			}
		}
	}

	{ // qualities
		uint8_t *qual_ptr = NULL;
		char qual;

		bam_new->data_len += bam_new->core.l_qseq;
		sw_align_bam_alloc_data(bam_new, bam_new->data_len);
		qual_ptr = bam1_qual(bam_new);

		if(space == SRMA_SPACE_CS) {
			// Get new base qualities based on color qualities
			for(i=0;i<bam_new->core.l_qseq;i++) {
				// use MAQ 0.7.1 conversion
				if(i == bam_new->core.l_qseq - 1) { // at the end of the alignment
					qual = srma_qual2char(color_qualities[i]);
				}
				else {
					int m1, m2;
					m1 = ((nt2int_table[(int)colors[i]] ^ nt2int_table[(int)colors[i+1]]) ==  nt4bit_to_int[bam1_seqi(bam1_seq(bam_new), i)]) ? 1 : 0;
					m2 = ((nt2int_table[(int)colors[i]] ^ nt2int_table[(int)colors[i+1]]) ==  nt4bit_to_int[bam1_seqi(bam1_seq(bam_new), i)]) ? 1 : 0;
					if(1 == m1 && 1 == m2) {
						qual = srma_char2qual(color_qualities[i]) + srma_char2qual(color_qualities[i+1]) + 10; 
					}
					else if(1 == m1) {
						qual = srma_char2qual(color_qualities[i]) - srma_char2qual(color_qualities[i+1]);
					}
					else if(1 == m2) {
						qual = srma_char2qual(color_qualities[i+1]) - srma_char2qual(color_qualities[i]);
					}
					else {
						qual = 1;
					}
				}
				bam1_qual(bam_new)[i] = __bound_qual(qual);
			}
		}
		else if(1 == correct_bases) {
			// Get new base qualities
			for(i=0;i<bam_new->core.l_qseq;i++) {
				if(bam1_seqi(bam1_seq(bam_new), i) == bam1_seqi(bam1_seq(bam_old), i)) {
					bam1_qual(bam_new)[i] = __bound_qual(bam1_qual(bam_old)[i]);
				}
				else {
					bam1_qual(bam_new)[i] = __bound_qual(bam1_qual(bam_old)[i] - SRMA_CORRECT_BASE_QUALITY_PENALTY);
				}
			}
		}
		else {
			// Copy old quality
			memcpy(bam1_qual(bam_new), bam1_qual(bam_old), bam_new->core.l_qseq);
		}
	}

	// TODO soft-clipping

	{ // Add in any auxiliary data as necessary
		uint8_t *s;
		bam_new->l_aux = 0;

		if(space == SRMA_SPACE_CS) {
			__copy_old("CS");
			__copy_old("CQ");
			__copy_old("RG");
		}
		// TODO 
		// XO/XQ
		// AS
		// XC
		// PG
	}


	// destroy the old bam structure
	bam_destroy1(bam_old);

	return bam_new;
}	
