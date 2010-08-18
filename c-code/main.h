#ifndef MAIN_H_
#define MAIN_H_

typedef struct {
	char **fn_inputs;
	int32_t fn_inputs_num;
	char **fn_outputs;
	int32_t fn_outputs_num;
	char *fn_ref;

	int32_t offset;
	int32_t min_mapq;
	double min_allele_prob;
	int32_t min_allele_coverage;
	int32_t max_total_coverage;
	char *fn_ranges;
	char *range;
	int32_t correct_bases;
	int32_t use_qualities;
	int32_t max_heap_size;
	int32_t max_queue_size;
	int32_t num_threads;
} srma_opt_t;

#endif
