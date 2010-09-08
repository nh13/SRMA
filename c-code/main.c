#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>
#include <config.h>

#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_util.h"
#include "srma_sam_io.h"
#include "srma_faidx.h"
#include "sw_align.h"
#include "graph.h"
#include "main.h"

#define THREAD_GRAPH_BLOCK_SIZE 256

static pthread_mutex_t thread_graph_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t thread_align_mutex = PTHREAD_MUTEX_INITIALIZER;

static void srma_core_output_progress(bam_record_t *bam_record, faidx_t *fai, int64_t records_processed, int32_t *max_msg_length)
{
	static char str[1024] = "\0";
	int32_t l;
	if(sprintf(str, "\rRecords processed: %lld (last %s:%u-%u)", (long long int)records_processed, srma_fai_name(fai, bam_record->b->core.tid), bam_record->b->core.pos+1, bam_calend(&bam_record->b->core, bam1_cigar(bam_record->b))+1) < 0) {
		srma_error(__func__, "Could not output progress", Exit, OutOfRange);
	}
	l = strlen(str);
	if((*max_msg_length) < l) {
		(*max_msg_length) = l;
	}
	else {
		while(l < (*max_msg_length)) {
			str[l]=' ';
			l++;
		}
		str[l]='\0';
	}
	if(fputs(str, stderr) == EOF) {
		srma_error(__func__, "Could not output progress", Exit, WriteFileError);
	}
}

static bam_record_t *srma_core_get_next_record(srma_sam_io_t *srma_sam_io, faidx_t *fai, srma_range_t *range_in, ref_t *ref)
{
	bam_record_t *bam_record = NULL;
	int32_t end;

	if(1 != srma_sam_io_has_next(srma_sam_io)) {
		return NULL;
	}

	bam_record = srma_sam_io_get_next(srma_sam_io);
	// check if we need more reference sequence
	end = bam_calend(&bam_record->b->core, bam1_cigar(bam_record->b));
	if(bam_record->b->core.pos < ref->beg && ref->end < end) {
		srma_fai_fetch(fai, ref, range_in->tid, bam_record->b->core.pos, end);
	}
	else if(bam_record->b->core.pos < ref->beg) {
		srma_fai_fetch(fai, ref, range_in->tid, bam_record->b->core.pos, ref->end);
	}
	else if(ref->end < end) {
		srma_fai_fetch(fai, ref, range_in->tid, ref->beg, end);
	}
	return bam_record;
}

typedef struct {
	bam_record_ll_t *to_graph_list;
	bam_record_ll_t *to_align_list;
	graph_t *graph;
	ref_t *ref;
	srma_ranges_t *ranges_out;
} thread_graph_data_t;

// Note: this must be used on the original bam record
static int srma_core_start_contained(bam_record_t *bam_record, srma_ranges_t *ranges_out)
{
	srma_range_t *range_out = NULL;
	range_out = srma_ranges_peek(ranges_out);
	if(NULL == range_out) { // no more ranges
		return 0;
	}

	while(range_out->tid < bam_record->b->core.tid
			|| (range_out->tid == bam_record->b->core.tid 
				&& range_out->end < bam_record->b->core.pos)) { // find a new range
		range_out = srma_ranges_poll(ranges_out);
		if(NULL == range_out) {
			return 0;
		}
	}	
	while(bam_record->b->core.tid < range_out->tid
			|| (range_out->tid == bam_record->b->core.tid 
				&& bam_record->b->core.pos < range_out->beg)) { // before the next range
		return 0;
	}
	else {
		return 1;
	}
}

static void srma_core_add_to_graph_worker(bam_record_ll_t *to_graph_list, bam_record_ll_t *to_align_list, graph_t *graph, ref_t *ref, srma_ranges_t *ranges_out, int32_t use_threads)
{
	int32_t i;

	do {
		bam_record_t *start=NULL, *cur= NULL;

		// --- SYNC ON --- 
		if(1 == use_threads) pthread_mutex_lock(&thread_graph_mutex);
		start = to_graph_list->iter;
		if(start == NULL) {
			// No more to process
			if(1 == use_threads) pthread_mutex_unlock(&thread_graph_mutex);
			break;
			// --- SYNC OFF --- 
		}
		else {
			// shift iterator
			for(i=0;i<THREAD_GRAPH_BLOCK_SIZE && NULL != to_graph_list->iter;i++) {
				to_graph_list->iter = to_graph_list->iter->next;
			}
			if(1 == use_threads) pthread_mutex_unlock(&thread_graph_mutex);
			// --- SYNC OFF --- 

			// process
			cur = start;
			for(i=0;i<THREAD_GRAPH_BLOCK_SIZE && NULL != cur;i++) {
				// remember to save the start node
				cur->node = graph_add_sam(graph, cur->b, ref, use_threads);
				assert(NULL != cur->node); // DEBUG
				cur = cur->next;
			}

			// add to align_list
			// --- SYNC ON --- 
			if(1 == use_threads) pthread_mutex_lock(&thread_align_mutex);
			cur = start;
			for(i=0;i<THREAD_GRAPH_BLOCK_SIZE && NULL != cur;i++) {
				if(graph->contig == cur->b->core.tid + 1 && 1 == srma_core_start_contained(cur, ranges_out)) {
					bam_record_ll_add1(to_align_list, cur->b, cur->fp_i, cur->node); // IMPORTANT: must remove from to_graph_list later 
				}
				else {
					// destroy bam 
					bam_destroy1(cur->b);
				}
				// remove potentially sensitive data
				cur->b = NULL;
				cur->fp_i = -1;
				cur->node = NULL;
				cur = cur->next;
			}
			if(1 == use_threads) pthread_mutex_unlock(&thread_align_mutex);
			// --- SYNC OFF --- 
		}
	} while(1);
}

void *srma_core_add_to_graph_thread(void *arg)
{
	thread_graph_data_t *data = (thread_graph_data_t*)arg;
	bam_record_ll_t *to_graph_list = data->to_graph_list;
	bam_record_ll_t *to_align_list = data->to_align_list;
	graph_t *graph = data->graph;
	ref_t *ref = data->ref;
	srma_ranges_t *ranges_out = data->ranges_out;

	srma_core_add_to_graph_worker(to_graph_list, to_align_list, graph, ref, ranges_out, 1);

	return arg;
}

static void srma_core_add_to_graph(srma_opt_t *opt, graph_t *graph, bam_record_ll_t *to_graph_list, bam_record_ll_t *to_align_list, ref_t *ref, srma_ranges_t *ranges_out)
{
	pthread_t *threads = NULL;
	thread_graph_data_t *data = NULL;
	int32_t i;
	void *status=NULL;

	if(0 < to_graph_list->size) {
		if(0 == to_align_list->size &&
				graph->contig != to_graph_list->head->b->core.tid) { // move to a new contig
			// prune and move to a new contig
			graph_prune(graph, to_graph_list->head->b->core.tid, to_graph_list->head->b->core.pos+1, opt->offset);
		}

		to_graph_list->iter = to_graph_list->head; // IMPORTANT: must be set

		if(1 < opt->num_threads) {
			// prepare thread data
			data = srma_malloc(sizeof(thread_graph_data_t)*opt->num_threads, __func__, "data");
			for(i=0;i<opt->num_threads;i++) {
				data[i].to_graph_list = to_graph_list;
				data[i].to_align_list = to_align_list;
				data[i].graph = graph;
				data[i].ref = ref;
				data[i].ranges_out = ranges_out;
			}

			threads = srma_malloc(sizeof(pthread_t)*opt->num_threads, __func__, "threads");

			// start threads
			for(i=0;i<opt->num_threads;i++) {
				if(0 != pthread_create(&threads[i], NULL, srma_core_add_to_graph_thread, &data[i])) {
					srma_error(__func__, "pthread_create", Exit, ThreadError);
				}
			}

			// wait for threads to finish
			for(i=0;i<opt->num_threads;i++) {
				if(0 != pthread_join(threads[i], &status)) {
					srma_error(__func__, "pthread_join", Exit, ThreadError);
				}
			}

			// free thread data
			free(threads);
			free(data);
		}
		else {
			srma_core_add_to_graph_worker(to_graph_list, to_align_list, graph, ref, ranges_out, 0);
		}

		// clear graph list
		bam_record_ll_clear(to_graph_list);
	}
}	

typedef struct {
	srma_opt_t *opt;
	bam_record_ll_t *to_align_list;
	graph_t *graph;
	cov_cutoffs_t *cov_cutoffs;
} thread_align_data_t;

static void srma_core_align_worker(srma_opt_t *opt, bam_record_ll_t *to_align_list, graph_t *graph, cov_cutoffs_t *cov_cutoffs, int32_t use_threads)
{
	sw_heap_t *heap=NULL;
	int32_t i;

	heap = sw_heap_init(SRMA_SW_HEAP_MIN, opt->max_heap_size);

	do {
		bam_record_t *start=NULL, *cur= NULL;

		// --- SYNC ON --- 
		if(1 == use_threads) pthread_mutex_lock(&thread_align_mutex);
		start = to_align_list->iter;
		if(start == NULL) {
			// No more to process
			if(1 == use_threads) pthread_mutex_unlock(&thread_align_mutex);
			break;
			// --- SYNC OFF --- 
		}
		else {
			// shift iterator
			for(i=0;i<THREAD_GRAPH_BLOCK_SIZE && NULL != to_align_list->iter;i++) {
				to_align_list->iter = to_align_list->iter->next;
			}
			if(1 == use_threads) pthread_mutex_unlock(&thread_align_mutex);
			// --- SYNC OFF --- 

			// process
			cur = start;
			for(i=0;i<THREAD_GRAPH_BLOCK_SIZE && NULL != cur;i++) {
				// remember to save the start node
				cur->b = sw_align(graph, cur->b, cur->node, heap, NULL, opt->offset, cov_cutoffs, opt->correct_bases, opt->use_qualities, opt->max_total_coverage, opt->max_heap_size);
				cur = cur->next;
			}
		}
	} while(1);

	sw_heap_free(heap);
}

void *srma_core_align_thread(void *arg)
{
	thread_align_data_t *data = (thread_align_data_t*)arg;
	srma_opt_t *opt = data->opt;
	bam_record_ll_t *to_align_list = data->to_align_list;
	graph_t *graph = data->graph;
	cov_cutoffs_t *cov_cutoffs = data->cov_cutoffs;

	srma_core_align_worker(opt, to_align_list, graph, cov_cutoffs, 1);

	return arg;
}

static void srma_core_align(srma_opt_t *opt, graph_t *graph, faidx_t *fai, bam_record_ll_t *to_align_list, bam_record_ll_t *to_output_list, cov_cutoffs_t *cov_cutoffs, srma_sam_io_t *srma_sam_io, int32_t flush, int64_t *records_processed, int32_t *max_msg_length)
{
	pthread_t *threads = NULL;
	thread_align_data_t *data = NULL;
	int32_t i;
	void *status = NULL;
	bam_record_t *bam_record = NULL, *bam_record_last = NULL;;

	if(0 < to_align_list->size) {
		//bam_calend(&b->core, bam1_cigar(b))
		if(0 != flush || bam_calend(&to_align_list->head->b->core, bam1_cigar(to_align_list->head->b)) + opt->offset < to_align_list->tail->b->core.pos) {

			to_align_list->iter = to_align_list->head; // IMPORTANT: must be set

			if(1 < opt->num_threads) {
				// prepare thread data
				data = srma_malloc(sizeof(thread_align_data_t)*opt->num_threads, __func__, "data");
				for(i=0;i<opt->num_threads;i++) {
					data[i].opt = opt;
					data[i].to_align_list = to_align_list;
					data[i].graph = graph;
					data[i].cov_cutoffs = cov_cutoffs;
				}

				threads = srma_malloc(sizeof(pthread_t)*opt->num_threads, __func__, "threads");

				// start threads
				for(i=0;i<opt->num_threads;i++) {
					if(0 != pthread_create(&threads[i], NULL, srma_core_align_thread, &data[i])) {
						srma_error(__func__, "pthread_create", Exit, ThreadError);
					}
				}

				// wait for threads to finish
				for(i=0;i<opt->num_threads;i++) {
					if(0 != pthread_join(threads[i], &status)) {
						srma_error(__func__, "pthread_join", Exit, ThreadError);
					}
				}

				free(threads);
				free(data);
			}
			else {
				srma_core_align_worker(opt, to_align_list, graph, cov_cutoffs, 0);
			}

			// add aligned bams to the output queue
			bam_record = to_align_list->head;
			while(bam_record != to_align_list->iter) {
				bam_record_last = bam_record;
				// add to the output queue
				bam_record_ll_add(to_output_list, bam_record_ll_remove(to_align_list));
				// move to next
				bam_record = to_align_list->head;
				// output the # that have been processed
				(*records_processed)++;
			}

			// prune the graph
			if(NULL != bam_record_last) {
				srma_core_output_progress(bam_record_last, fai, (*records_processed), max_msg_length);
				if(0 < to_align_list->size) {
					graph_prune(graph, to_align_list->head->b->core.tid, to_align_list->head->b->core.pos+1, opt->offset);
				}
				else {
					graph_prune(graph, bam_record_last->b->core.tid, bam_record_last->b->core.pos+1, opt->offset); 
				}
			}
		}
	}
	if(0 != flush && NULL != bam_record_last) {
		srma_core_output_progress(bam_record_last, fai, (*records_processed), max_msg_length);
	}
	bam_record_last = NULL;

	// output alignments
	while(0 < to_output_list->size) {
		bam_record = to_output_list->head; // peek
		// alignment could have moved (+OFFSET), with another moving (-OFFSET) 
		if(0 != flush || 
				bam_record->b->core.tid < graph->contig || // other alignments will be less than
				(bam_record->b->core.pos + 1 + 2*opt->offset) < graph->position_start) { // other alignments will not be less than
			bam_record = bam_record_ll_remove(to_output_list);
			// write to file
			srma_sam_io_write(srma_sam_io, bam_record);
			// free bam_record
			bam_record_free(bam_record);
		}
		else { // other alignments could be less than
			break;
		}
	}
}

static void srma_core(srma_opt_t *opt)
{
	clock_t start_time, end_time;
	cov_cutoffs_t *cov_cutoffs = NULL;
	faidx_t *fai = NULL;
	bam_record_ll_t *to_graph_list = NULL, *to_align_list = NULL, *to_output_list = NULL;;
	srma_ranges_t *ranges_in=NULL, *ranges_out=NULL;
	srma_range_t *range_in=NULL;
	srma_sam_io_t *srma_sam_io = NULL;
	ref_t *ref = NULL;
	graph_t *graph = NULL;
	bam_record_t *bam_record = NULL;
	int32_t prev_tid = -1, prev_pos = -1; // TODO: reset when we move to a new range
	int64_t records_processed=0;
	int32_t max_msg_length=0;

	start_time = clock();

	// reference sequence
	fai = fai_load(opt->fn_ref);
	if(NULL == fai) {
		srma_error(__func__, "Problem opening FASTA reference and FASTA index", Exit, OpenFileError);
	}

	// range/ranges
	if(NULL != opt->range && NULL != opt->fn_ranges) {
		srma_error(__func__, "-R and -Z were both specified", Exit, CommandLineArgument);
	}
	else if(NULL == opt->range && NULL == opt->fn_ranges) {
		// add ranges from the fasta index
		ranges_in = srma_ranges_init_from_fai(fai);
		ranges_out = srma_ranges_init_from_fai(fai);
		// init inputs/ouputs/indexes
		srma_sam_io = srma_sam_io_init(opt->fn_inputs, opt->fn_inputs_num, opt->fn_outputs, opt->fn_outputs_num, opt->fn_output_header, 0);
	}
	else {
		// Use ranges
		ranges_in = srma_ranges_init(fai, opt->fn_ranges, opt->range, opt->offset);
		ranges_out = srma_ranges_init(fai, opt->fn_ranges, opt->range, 0);
		// init inputs/ouputs/indexes
		srma_sam_io = srma_sam_io_init(opt->fn_inputs, opt->fn_inputs_num, opt->fn_outputs, opt->fn_outputs_num, opt->fn_output_header, 1);
	}

	// init
	cov_cutoffs = cov_cutoffs_init(opt->min_allele_coverage, opt->min_allele_prob);
	to_output_list = bam_record_ll_init();

	// Get first range
	range_in = srma_ranges_peek(ranges_in);
	while(NULL != range_in) {

		// init structures
		prev_tid = prev_pos=-1;
		to_graph_list = bam_record_ll_init();
		to_align_list = bam_record_ll_init();
		graph = graph_init();
		ref = srma_malloc(sizeof(ref_t), __func__, "ref");
		ref->ref = NULL;

		// move to first range
		srma_sam_io_change_range(srma_sam_io, range_in->tid, range_in->beg, range_in->end);
		// Get first reference sequence
		srma_fai_fetch(fai, ref, range_in->tid, range_in->beg, range_in->end); 

		// get first record
		bam_record = srma_core_get_next_record(srma_sam_io, fai, range_in, ref);

		// initial prune
		if(NULL != bam_record) { 
			graph_prune(graph, bam_record->b->core.tid, bam_record->b->core.pos+1, opt->offset);
		}

		while(NULL != bam_record) { // while records to process in this interval
			if(bam_record->b->core.flag & BAM_FUNMAP) {
				// removed from output
				bam_record_free(bam_record);
				bam_record = srma_core_get_next_record(srma_sam_io, fai, range_in, ref);
				continue;
			}
			else if(bam_record->b->core.qual < opt->min_mapq) {
				// removed from output
				bam_record_free(bam_record);
				bam_record = srma_core_get_next_record(srma_sam_io, fai, range_in, ref);
				continue;
			}

			// check sorted
			if(bam_record->b->core.tid < prev_tid ||
					(bam_record->b->core.tid == prev_tid && bam_record->b->core.pos < prev_pos)) {
				fprintf(stderr, "QNAME=%s\n", bam1_qname(bam_record->b));
				srma_error(__func__, "SAM/BAM file is not co-ordinate sorted", Exit, OutOfRange);
			}
			prev_tid = bam_record->b->core.tid;
			prev_pos = bam_record->b->core.pos;

			// process graph
			if(opt->max_queue_size <= to_graph_list->size) {
				// add all records left to the graph
				srma_core_add_to_graph(opt, graph, to_graph_list, to_align_list, ref, ranges_out);
			}

			// add to the graph list
			bam_record_ll_add(to_graph_list, bam_record);

			// align
			if(opt->max_queue_size <= to_align_list->size) {
				srma_core_align(opt, graph, fai, to_align_list, to_output_list, cov_cutoffs, srma_sam_io, 0, &records_processed, &max_msg_length);
			}

			bam_record = srma_core_get_next_record(srma_sam_io, fai, range_in, ref);
		}
		// process graph
		if(0 < to_graph_list->size) {
			// add all records left to the graph
			srma_core_add_to_graph(opt, graph, to_graph_list, to_align_list, ref, ranges_out);
		}
		// process alignments
		if(0 < to_align_list->size) {
			// align all records left
			srma_core_align(opt, graph, fai, to_align_list, to_output_list, cov_cutoffs, srma_sam_io, 1, &records_processed, &max_msg_length);
		}

		// free
		graph_free(graph);
		graph=NULL;
		bam_record_ll_destroy(to_graph_list);
		to_graph_list=NULL;
		bam_record_ll_destroy(to_align_list);
		to_align_list=NULL;
		if(NULL != ref->ref) {
			free(ref->ref);
		}
		free(ref);
		ref=NULL;
		
		// destroy input buffer etc.
		srma_sam_io_close_inputs(srma_sam_io);

		// next input range
		range_in = srma_ranges_poll(ranges_in);
	}
	// flush the output queue
	while(0 < to_output_list->size) {
		bam_record = bam_record_ll_remove(to_output_list);
		// write to file
		srma_sam_io_write(srma_sam_io, bam_record);
		// free bam_record
		bam_record_free(bam_record);
	}
	fprintf(stderr, "\n"); // check return value ? 

	// free
	srma_ranges_free(ranges_in);
	srma_ranges_free(ranges_out);
	bam_record_ll_destroy(to_output_list);
	fai_destroy(fai);
	srma_sam_io_free(srma_sam_io);
	cov_cutoffs_free(cov_cutoffs);

	// timing 
	end_time = clock();
	fprintf(stderr, "Elapsed time: %.2f sec\n", (float)(end_time - start_time) / CLOCKS_PER_SEC); 
}

static inline char **add_file(char **fns, int32_t fns_num, char *fn)
{
	fns = srma_realloc(fns, sizeof(char*)*(1+fns_num), __func__, "fns");
	fns[fns_num] = srma_strdup(fn, __func__);
	return fns;
}

static int print_usage(srma_opt_t *opt)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: srma (short-read micro re-aligner)\n");
	fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
	fprintf(stderr, "Usage:   srma [options]\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options (required):\n");
	fprintf(stderr, "         -i FILE     input sam/bam\n");
	fprintf(stderr, "         -o FILE     output sam/bam\n");
	fprintf(stderr, "         -r FILE     reference FASTA file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Options (optional):\n");
	fprintf(stderr, "         -O INT      alignment offset [%d]\n", opt->offset);
	fprintf(stderr, "         -m INT      minimum mapping quality [%d]\n", opt->min_mapq);
	fprintf(stderr, "         -p FLOAT    minimum allele probability conditioned on coverage [%lf]\n", opt->min_allele_prob);
	fprintf(stderr, "         -c INT      minimum hapoloid coverage for the consensus [%d]\n", opt->min_allele_coverage);
	fprintf(stderr, "         -t INT      maximum total coverage over a reference base [%d]\n", opt->max_total_coverage);
	fprintf(stderr, "         -R STRING   an genomic range to consider [%s]\n", opt->range);
	fprintf(stderr, "         -Z FILE     input genomic ranges to consider [%s]\n", opt->fn_ranges);
	fprintf(stderr, "         -C INT      correct aligned bases [%d]\n", opt->correct_bases);
	fprintf(stderr, "         -q INT      use sequence qualities to weight alignments [%d]\n", opt->use_qualities);
	fprintf(stderr, "         -H INT      maximum heap size [%d]\n", opt->max_heap_size);
	fprintf(stderr, "         -Q INT      maximum queue size [%d]\n", opt->max_queue_size);
	fprintf(stderr, "         -n INT      number of threads [%d]\n", opt->num_threads);
	fprintf(stderr, "         -b FILE     copy the header in FILE to the output SAM/BAM (multiple input, single output only)\n");
	fprintf(stderr, "         -h          print this message\n");
	fprintf(stderr, "\n");

	return 1;
}

int main(int argc, char *argv[])
{

	srma_opt_t opt;
	int c, i;

	// required
	opt.fn_inputs = NULL;
	opt.fn_inputs_num = 0;
	opt.fn_outputs = NULL;
	opt.fn_outputs_num = 0;
	opt.fn_output_header = NULL;
	opt.fn_ref = NULL;
	// optional
	opt.offset = 20;
	opt.min_mapq = 0;
	opt.min_allele_prob = 0.1;
	opt.min_allele_coverage = 3;
	opt.max_total_coverage = 100;
	opt.fn_ranges = NULL;
	opt.range = NULL;
	opt.correct_bases = 0;
	opt.use_qualities = 1;
	opt.max_heap_size = 8192;
	opt.max_queue_size = 65536;
	opt.num_threads = 1;

	while(0 <= (c = getopt(argc, argv, "i:o:r:O:m:p:c:t:R:Z:C:qHQ:n:bh"))) {
		switch(c) {
			case 'i':
				opt.fn_inputs = add_file(opt.fn_inputs, opt.fn_inputs_num, optarg); 
				opt.fn_inputs_num++;
				break;
			case 'o':
				opt.fn_outputs = add_file(opt.fn_outputs, opt.fn_outputs_num, optarg);
				opt.fn_outputs_num++;
				break;
			case 'r':
				opt.fn_ref = srma_strdup(optarg, __func__); break;
			case 'O':
				opt.offset = atoi(optarg); break;
			case 'm':
				opt.min_mapq = atoi(optarg); break;
			case 'p':
				opt.min_allele_prob = atof(optarg); break;
			case 'c':
				opt.min_allele_coverage = atoi(optarg); break;
			case 't':
				opt.max_total_coverage = atoi(optarg); break;
			case 'R':
				opt.range = srma_strdup(optarg, __func__); break;
			case 'Z':
				opt.fn_ranges = srma_strdup(optarg, __func__); break;
			case 'C':
				opt.correct_bases = atoi(optarg); break;
			case 'q':
				opt.use_qualities = atoi(optarg); break;
			case 'H':
				opt.max_heap_size = atoi(optarg); break;
			case 'Q':
				opt.max_queue_size = atoi(optarg); break;
			case 'n':
				opt.num_threads = atoi(optarg); break;
			case 'b':
				opt.fn_output_header = srma_strdup(optarg, __func__); break;
			case 'h':
			default:
				return print_usage(&opt);
		}
	}
	if(argc != optind || 1 == argc) {
		return print_usage(&opt);
	}

	for(i=0;i<argc;i++) {
		if(0 < i) fputc(' ', stderr);
		fprintf(stderr, "%s", argv[i]);
	}
	fprintf(stderr, "\n");

	// Check args
	// TODO thorough bounds checking
	if(0 == opt.fn_inputs_num) {
		srma_error(__func__, "No inputs (-i)", Exit, CommandLineArgument);
	}
	if(0 == opt.fn_outputs_num) { // TODO could have stdout
		srma_error(__func__, "No outputs (-o)", Exit, CommandLineArgument);
	}
	if(NULL == opt.fn_ref) {
		srma_error(__func__, "No reference FASTA (-r)", Exit, CommandLineArgument);
	}
	if(NULL != opt.range && NULL != opt.fn_ranges) {
		srma_error(__func__, "-R and -Z were both specified", Exit, CommandLineArgument);
	}
	if(opt.fn_inputs_num != opt.fn_outputs_num && 1 != opt.fn_outputs_num) {
		srma_error(__func__, "The same number of inputs and outputs must be specified, or only one output", Exit, CommandLineArgument);
	}
	if(opt.fn_inputs_num > 1 && 1 == opt.fn_outputs_num && NULL == opt.fn_output_header) {
		srma_error(__func__, "Multiple inputs and a single output requires the -b option", Exit, CommandLineArgument);
	}
	if(1 < opt.num_threads) {
		fprintf(stderr, "***** Warning: multiple threads may not increase performance,    ******\n"); 
		fprintf(stderr, "***** and may actually derease peformance.  Try running multiple ******\n");
		fprintf(stderr, "***** processes using the -R option.                             ******\n");
	}

	srma_core(&opt);

	for(i=0;i<opt.fn_inputs_num;i++) {
		free(opt.fn_inputs[i]);
	}
	free(opt.fn_inputs);
	for(i=0;i<opt.fn_outputs_num;i++) {
		free(opt.fn_outputs[i]);
	}
	free(opt.fn_outputs);
	free(opt.fn_output_header);
	free(opt.fn_ref);
	free(opt.fn_ranges);
	free(opt.range);

	return 0;
}
