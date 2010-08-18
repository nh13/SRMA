#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include "samtools/faidx.h"
#include "samtools/sam.h"
#include "samtools/bam.h"
#include "srma_alloc.h"
#include "srma_error.h"
#include "srma_sam_io.h"

#define __is_sam(_fn) ((0 == strcmp(".sam", _fn + strlen(_fn) - 4)) ? 1 : 0) 

void bam_record_free(bam_record_t *bam_record)
{
	if(NULL != bam_record->b) { // it could have been removed
		bam_destroy1(bam_record->b);
	}
	free(bam_record);
}

bam_record_ll_t *bam_record_ll_init()
{
	bam_record_ll_t *ll = NULL;

	ll = srma_malloc(sizeof(bam_record_ll_t), __func__, "ll->");
	ll->head = ll->tail = NULL;
	ll->size = 0;

	return ll;
}

void bam_record_ll_add(bam_record_ll_t *ll, bam_record_t *b)
{
	if(NULL == ll->head) {
		ll->tail = ll->head = b;
		ll->head->next = ll->head->prev = NULL;
	}
	else {
		bam_record_t *rec=NULL;

		ll->tail->next = b;
		ll->tail->next->prev = ll->tail;
		ll->tail->next->next = NULL;
		ll->tail = ll->tail->next;

		rec = ll->tail;
		while(NULL != rec->prev) {

			if(rec->prev->b->core.tid <= rec->b->core.tid ||
					(rec->prev->b->core.tid == rec->b->core.tid &&
					 rec->prev->b->core.pos <= rec->b->core.pos)) {
				break;
			}
			// swap
			bam_record_t *prev = rec->prev;
			rec->prev = prev->prev;
			prev->next = rec->next;
			rec->next = prev;
			prev->prev = rec;
		}
	}
	ll->size++;
}

void bam_record_ll_add1(bam_record_ll_t *ll, bam1_t *b, int32_t fp_i, node_t *node)
{
	bam_record_t *bam_record = NULL;
	bam_record = srma_malloc(sizeof(bam_record_t), __func__, "bam_record");
	bam_record->b = b;
	bam_record->node = node;
	bam_record->fp_i = fp_i;
	bam_record->prev = bam_record->next = NULL;
	bam_record_ll_add(ll, bam_record);
}

bam_record_t *bam_record_ll_remove(bam_record_ll_t *ll)
{
	bam_record_t *ret = ll->head;
	if(NULL == ll || NULL == ll->head) {
		return NULL;
	}
	ll->head = ll->head->next;
	ll->size--;
	if(NULL == ll->head) {
		assert(0 == ll->size); // DEBUG
		ll->tail = NULL;
	}
	else {
		ll->head->prev = NULL;
	}
	ret->next = ret->prev = NULL;
	return ret;
}

void bam_record_ll_clear(bam_record_ll_t *ll)
{
	while(ll->head != NULL) {
		bam_record_t *rec = ll->head;
		ll->head = ll->head->next;
		bam_record_free(rec);
	}
	ll->head = ll->tail = NULL;
	ll->size = 0;
}

void bam_record_ll_destroy(bam_record_ll_t *ll)
{
	bam_record_ll_clear(ll);
	free(ll);
}

// read in bam structure
static inline bam1_t *srma_sam_io_read(srma_sam_io_t *s, int32_t fp_i)
{
	bam1_t *b = NULL;

	b = srma_calloc(1, sizeof(bam1_t), __func__, "b");

	switch(s->fps_in_type[fp_i]) {
		case SRMA_SAM_IO_TYPE_SAM:
		case SRMA_SAM_IO_TYPE_BAM:
			if(0 < samread(s->fps_in[fp_i], b)) return b;
			break;
		case SRMA_SAM_IO_TYPE_BAM_ITER:
			// DEBUG
			assert(s->fps_in[fp_i]->type & 0x01); // should be a BAM
			if(NULL != s->bam_iters[fp_i] && 0 <= bam_iter_read(s->fps_in[fp_i]->x.bam, s->bam_iters[fp_i], b)) return b;
			break;
		default:
			srma_error(__func__, "SAM/BAM type", Exit, OutOfRange);
	}
	bam_destroy1(b);
	return NULL;
}

// zero-based co-ordinates
srma_sam_io_t *srma_sam_io_init(char **fn_inputs, int32_t fn_inputs_num, char **fn_outputs, int32_t fn_outputs_num, int32_t use_ranges)
{
	int32_t i;
	srma_sam_io_t *s= NULL;

	assert(fn_inputs_num == fn_outputs_num || fn_outputs_num == 1);

	s = srma_malloc(sizeof(srma_sam_io_t), __func__, "s");

	s->fps_in_num = fn_inputs_num;
	s->fps_in = srma_malloc(sizeof(samfile_t*)*s->fps_in_num, __func__, "s->fps_in");
	s->fps_in_type = srma_malloc(sizeof(int32_t)*s->fps_in_num, __func__, "s->fps_in_type");
	s->fps_out_num = fn_outputs_num;
	s->fps_out = srma_malloc(sizeof(samfile_t*)*s->fps_out_num, __func__, "s->fps_out");

	for(i=0;i<s->fps_in_num;i++) {
		s->fps_in_type[i] = (1 == __is_sam(fn_inputs[i])) ? SRMA_SAM_IO_TYPE_SAM : SRMA_SAM_IO_TYPE_BAM;
		s->fps_in[i] = samopen(fn_inputs[i], __is_sam(fn_inputs[i]) ? "r" : "rb", 0);
		if(NULL == s->fps_in[i]) {
			srma_error(__func__, fn_inputs[i], Exit, OpenFileError);
		}
		if(1 == s->fps_in_num || 1 < s->fps_out_num) { // from a single file or to multiple files
			s->fps_out[i] = samopen(fn_outputs[i], __is_sam(fn_outputs[i]) ? "wh" : "wb", s->fps_in[i]->header);
			if(NULL == s->fps_out[i]) {
				srma_error(__func__, fn_outputs[i], Exit, OpenFileError);
			}
		}

		// TODO: merged header ?
		// TODO: program records ? 
		// TODO: set coordinate sorted ?
	}

	if(1 < s->fps_in_num && 1 == s->fps_out_num) { // to a single file and from multiple files
		// TODO: use a merged header instead of just the first one 
		srma_error(__func__, "NOT IMPLEMENTED", Exit, OutOfRange);
		s->fps_out[0] = samopen(fn_outputs[0], __is_sam(fn_outputs[0]) ? "wh" : "wb", s->fps_in[0]->header);
	}

	if(0 == use_ranges) {
		s->bam_iters = NULL;
		s->bam_indexes = NULL;
	}
	else {
		s->bam_indexes = srma_malloc(sizeof(bam_index_t*)*s->fps_in_num, __func__, "s->bam_indexes");
		s->bam_iters = srma_malloc(sizeof(bam_iter_t)*s->fps_in_num, __func__, "s->bam_indexes");
		for(i=0;i<s->fps_in_num;i++) {
			s->bam_iters[i] = NULL;
			switch(s->fps_in_type[i]) {
				case SRMA_SAM_IO_TYPE_BAM:
					// load the index
					s->fps_in_type[i] = SRMA_SAM_IO_TYPE_BAM_ITER; // fall through
				case SRMA_SAM_IO_TYPE_BAM_ITER:
					s->bam_indexes[i] = bam_index_load(fn_inputs[i]);
					if(NULL == s->bam_indexes[i]) {
						srma_error(__func__, fn_inputs[i], Exit, OpenFileError);
					}
					break;
				case SRMA_SAM_IO_TYPE_SAM:
					srma_error(__func__, "Range querying is not possible with SAM files", Exit, OutOfRange);
					break;
				default:
					srma_error(__func__, "SAM/BAM type", Exit, OutOfRange);
			}
		}
	}

	// empty buffer
	s->buffer = bam_record_ll_init();

	return s;
}

inline int32_t srma_sam_io_has_next(srma_sam_io_t *s) 
{
	return (0 ==  s->buffer->size) ? 0 : 1;
}

bam_record_t *srma_sam_io_get_next(srma_sam_io_t *s)
{
	if(1 == srma_sam_io_has_next(s)) {
		bam_record_t *ret = bam_record_ll_remove(s->buffer);
		// add one to buffer from the same file
		if(NULL != s->fps_in[ret->fp_i]) { // not closed
			bam1_t *b = srma_sam_io_read(s, ret->fp_i);
			if(NULL != b) { // has more
				bam_record_ll_add1(s->buffer, b, ret->fp_i, NULL);
			}
		}
		return ret;
	}
	else {
		return NULL;
	}
}

void srma_sam_io_init_buffer(srma_sam_io_t *s)
{
	int32_t i;
	// initialize the buffer
	for(i=0;i<s->fps_in_num;i++) {
		bam1_t *b = NULL;
		b = srma_sam_io_read(s, i);
		if(NULL != b) {
			bam_record_ll_add1(s->buffer, b, i, NULL);
		}
	}
}

// zero-based
void srma_sam_io_change_range(srma_sam_io_t *s, int32_t sequence_index, int32_t position_start, int32_t position_end)
{
	int32_t i;

	// query
	for(i=0;i<s->fps_in_num;i++) {
		switch(s->fps_in_type[i]) {
			case SRMA_SAM_IO_TYPE_BAM:
				s->fps_in_type[i] = SRMA_SAM_IO_TYPE_BAM_ITER; // fall through
			case SRMA_SAM_IO_TYPE_BAM_ITER:
				assert(NULL != s->bam_indexes[i]); // DEBUG
				s->bam_iters[i] = bam_iter_query(s->bam_indexes[i], sequence_index, position_start, position_end); 
				break;
			case SRMA_SAM_IO_TYPE_SAM:
				// ignore
				break;
			default:
				srma_error(__func__, "SAM/BAM type", Exit, OutOfRange);
		}
	}

	// initialize the buffer
	srma_sam_io_init_buffer(s);
}


// clear buffer and destroy the iterators
// TODO: only need to if haven't closed before
void srma_sam_io_close_inputs(srma_sam_io_t *s)
{
	int32_t i;

	// clear buffer
	bam_record_ll_clear(s->buffer);

	// query
	for(i=0;i<s->fps_in_num;i++) {
		if(SRMA_SAM_IO_TYPE_BAM_ITER == s->fps_in_type[i] && NULL != s->bam_iters[i]) { // BAM
			bam_iter_destroy(s->bam_iters[i]);
			s->bam_iters[i] = NULL;
		}
	}
}

void srma_sam_io_write(srma_sam_io_t *s, bam_record_t *b)
{
	if(1 == s->fps_out_num) {
		if(samwrite(s->fps_out[0], b->b) < 0) {
			srma_error(__func__, "samwrite", Exit, WriteFileError);
		}
	}
	else {
		if(samwrite(s->fps_out[b->fp_i], b->b) < 0) {
			srma_error(__func__, "samwrite", Exit, WriteFileError);
		}
	}
}

void srma_sam_io_free(srma_sam_io_t *s)
{
	int32_t i;
	// destroy buffer and iterators
	srma_sam_io_close_inputs(s);
	bam_record_ll_destroy(s->buffer);
	// close the input/output files as well as the indexes
	for(i=0;i<s->fps_in_num;i++) {
		samclose(s->fps_in[i]);
		switch(s->fps_in_type[i]) {
			case SRMA_SAM_IO_TYPE_SAM:
			case SRMA_SAM_IO_TYPE_BAM: 
				break;
			case SRMA_SAM_IO_TYPE_BAM_ITER:
				bam_index_destroy(s->bam_indexes[i]);
				break;
			default:
				srma_error(__func__, "SAM/BAM type", Exit, OutOfRange);
		}
	}
	for(i=0;i<s->fps_out_num;i++) {
		samclose(s->fps_out[i]);
	}
	free(s->fps_in);
	free(s->fps_out);
	free(s->fps_in_type);
	free(s->bam_indexes);
	free(s->bam_iters);
	free(s);
}
