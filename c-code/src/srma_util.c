#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_util.h"

static void bam_aln_left_justify(bam_aln_t *aln)
{
	int32_t i, prev_del, prev_ins, start_del, end_del, start_ins, end_ins;

	i = prev_del = prev_ins = 0;
	start_del = end_del = start_ins = end_ins = -1;

	while(i<aln->length) {
		assert(0 == prev_ins || 0 == prev_del);

		if(ALN_GAP == aln->read[i]) {
			if(0 == prev_del) {
				start_del = i;
			}
			prev_del = 1;
			end_del = i;
			prev_ins = 0;
			start_ins = -1;
			end_ins = -1;
			i++;
		}
		else if(ALN_GAP == aln->ref[i]) {
			if(0 == prev_ins) {
				start_ins = i;
			}
			prev_ins = 1;
			end_ins = i;
			prev_del = 0;
			start_del = -1;
			end_del = -1;
			i++;
		}
		else {
			if(1 == prev_del) {
				assert(0 < start_del);
				assert(start_del <= end_del);
				start_del--;
				while(0 <= start_del && // Bases remaining to examine 
						aln->read[start_del] != ALN_GAP && // Hit another deletion 
						aln->ref[start_del] != ALN_GAP && // Hit an insertion 
						aln->ref[start_del] == aln->ref[end_del]) { // src ref base matches dest ref base 
					assert(ALN_GAP != aln->ref[start_del]);
					assert(ALN_GAP != aln->ref[end_del]);
					assert(ALN_GAP != aln->read[start_del]);
					assert(ALN_GAP == aln->read[end_del]);
					aln->read[end_del] = aln->read[start_del];
					aln->read[start_del] = ALN_GAP;
					start_del--;
					end_del--;
				}
				end_del++; // We decremented when we exited the loop 
				i = end_del;
				assert(ALN_GAP != aln->read[i]);
				assert(ALN_GAP != aln->ref[i]);
			}
			else if(1 == prev_ins) {
				assert(start_ins <= end_ins);
				start_ins--;
				while(0 <= start_ins && // Bases remaining to examine 
						aln->read[start_ins] != ALN_GAP && // Hit another deletion 
						aln->ref[start_ins] != ALN_GAP && // Hit an insertion 
						aln->read[start_ins] == aln->read[end_ins]) { // src aln->read base matches dest aln->read base 
					assert(ALN_GAP != aln->read[start_ins]);
					assert(ALN_GAP != aln->read[end_ins]);
					assert(ALN_GAP != aln->ref[start_ins]);
					assert(ALN_GAP == aln->ref[end_ins]);
					aln->ref[end_ins] = aln->ref[start_ins];
					aln->ref[start_ins] = ALN_GAP;
					start_ins--;
					end_ins--;
				}
				end_ins++; // We decremented when we exited the loop 
				i = end_ins;
				assert(ALN_GAP != aln->read[i]);
				assert(ALN_GAP != aln->ref[i]);
			}
			else {
				i++;
			}
			prev_del = 0;
			prev_ins = 0;
			start_del = -1;
			end_del = -1;
			start_ins = -1;
			end_ins = -1;
		}
	}
}


bam_aln_t *bam_aln_init(bam1_t *bam, ref_t *ref)
{

	bam_aln_t *aln = NULL;
	int32_t i, j, aln_index, read_index, ref_index;
	int32_t op, l;

	aln = srma_malloc(sizeof(bam_aln_t), __func__, "aln");

	// Get alignment length
	aln->length = aln->positions_length = 0;
	for(i=0;i<bam->core.n_cigar;i++) {
		op = bam1_cigar(bam)[i] & BAM_CIGAR_MASK; // operation
		l = bam1_cigar(bam)[i] >> BAM_CIGAR_SHIFT; // length
		switch(op) {
			case BAM_CMATCH:
			case BAM_CINS:
				aln->length += l;
				aln->positions_length += l;
				break;
			case BAM_CDEL:
				aln->length += l;
				break;
			default:
				break; // ignore
		}
	}
	aln->read = srma_calloc((aln->length+1), sizeof(uint8_t), __func__, "aln->read");
	aln->ref = srma_calloc((aln->length+1), sizeof(uint8_t), __func__, "aln->ref");

	aln_index = read_index = 0;
	ref_index = bam->core.pos - ref->beg; // zero-based 
	assert(ref->beg <= bam->core.pos); // DEBUG
	for(i=0;i<bam->core.n_cigar;i++) {
		op = bam1_cigar(bam)[i] & BAM_CIGAR_MASK; // operation
		l = bam1_cigar(bam)[i] >> BAM_CIGAR_SHIFT; // length
		switch(op) {
			case BAM_CMATCH:
				for(j=0;j<l;j++) {
					/*
					if(ref_index < 0 || strlen(ref->ref) <= ref_index) {
						fprintf(stderr, "0 <= %d <= %d\n", ref_index, (int)strlen(ref->ref));
					}
					assert(0 <= ref_index); // DEBUG
					assert(ref_index < strlen(ref->ref)); // DEBUG
					*/
					aln->ref[aln_index] = ref->ref[ref_index];
					aln->read[aln_index] = bam_nt16_rev_table[bam1_seqi(bam1_seq(bam), read_index)];
					ref_index++;
					read_index++;
					aln_index++;
				}
				break;
			case BAM_CINS:
				for(j=0;j<l;j++) {
					aln->ref[aln_index] = ALN_GAP;
					aln->read[aln_index] = bam_nt16_rev_table[bam1_seqi(bam1_seq(bam), read_index)];
					read_index++;
					aln_index++;
				}
				break;
			case BAM_CDEL:
				for(j=0;j<l;j++) {
					//assert(0 <= ref_index); // DEBUG
					//assert(ref_index < strlen(ref->ref)); // DEBUG
					aln->ref[aln_index] = ref->ref[ref_index];
					aln->read[aln_index] = ALN_GAP;
					ref_index++;
					aln_index++;
				}
				break;
			case BAM_CSOFT_CLIP:
				read_index++;
				break;
			case BAM_CHARD_CLIP:
				break;
			default:
				srma_error(__func__, "op", Exit, OutOfRange);
		}
	}
	aln->read[aln->length] = aln->ref[aln->length] = '\0';

	// left justify indels
	if(read_index < aln_index || ref_index < aln_index) {
		bam_aln_left_justify(aln);
	}

	// set positions
	aln->positions = srma_calloc(aln->positions_length, sizeof(uint32_t), __func__, "aln->positions");
	aln->positions_index = srma_calloc(aln->positions_length, sizeof(uint32_t), __func__, "aln->positions_index");
	aln_index = read_index = 0;
	ref_index = -1;
	while(aln_index < aln->length) {
		// Previous not an insertion
		if(0 == aln_index || ALN_GAP != aln->ref[aln_index-1]) {
			ref_index++;
		}
		if(ALN_GAP != aln->read[aln_index]) { // Not a deletion
			aln->positions[read_index] = ref_index;
			aln->positions_index[read_index] = aln_index;
			read_index++;
		}
		aln_index++;
	}

	return aln;
}

void bam_aln_free(bam_aln_t *aln)
{
	free(aln->read);
	free(aln->ref);
	free(aln->positions);
	free(aln->positions_index);
	free(aln);
}

inline int32_t srma_char2qual(char qual)
{
	int32_t q = (((int32_t)qual) - 33);
	if(q < 0) {
		return 1;
	}
	else if(255 < q) {
		return 255;
	}
	else {
		return q;
	}
}

inline char srma_qual2char(int32_t qual)
{
	return (char)(qual + 33);
}

cov_cutoffs_t *cov_cutoffs_init(int32_t min_allele_coverage, double min_allele_prob)
{
	cov_cutoffs_t *c = NULL;
	int32_t **bin_coeff = NULL;
	int32_t bin_coeff_rows = 0;
	int32_t i, cur_coverage, last_coverage;

	c = srma_malloc(sizeof(cov_cutoffs_t), __func__, "c");

	c->min_allele_coverage = srma_malloc(sizeof(int32_t), __func__, "c->min_allele_coverage");
	c->min_allele_coverage[0] = 0; // dummy for 0
	c->length = 1;
	c->max_coverage = 0;

	bin_coeff = srma_malloc(sizeof(int32_t*), __func__, "bin_coeff");
	bin_coeff[0] = NULL;
	bin_coeff_rows = 1;

	cur_coverage = 1;
	last_coverage = 0;

	fprintf(stderr, "Allele coverage cutoffs:\n");

	while(last_coverage < min_allele_coverage) {
		int32_t bin_coeff_cols = 0;

		bin_coeff_rows++;
		bin_coeff = srma_realloc(bin_coeff, sizeof(int32_t*)*bin_coeff_rows, __func__, "bin_coeff");
		bin_coeff[bin_coeff_rows-1] = NULL;

		// dynamic programming for binomial co-efficients
		int32_t i;
		bin_coeff_cols = cur_coverage+1;
		bin_coeff[bin_coeff_rows-1] = srma_malloc(sizeof(int32_t)*bin_coeff_cols, __func__, "bin_coeff[bin_coeff_rows-1]");
		for(i=0;i<=cur_coverage;i++) {
			if(0 == i || cur_coverage == i) {
				bin_coeff[bin_coeff_rows-1][i] = 1;
			}
			else {
				bin_coeff[bin_coeff_rows-1][i] = bin_coeff[cur_coverage-1][i-1] + bin_coeff[cur_coverage-1][i];
			}
		}

		double p = 0.0;
		double p2 = pow(0.5, cur_coverage);
		int32_t ctr = -1; // will always be incremented
		do {
			ctr++;
			p += p2 * bin_coeff[cur_coverage][ctr];
		} while(p < min_allele_prob);

		last_coverage = ctr;
		c->length++;
		c->min_allele_coverage = srma_realloc(c->min_allele_coverage, c->length*sizeof(int32_t), __func__, "c->min_allele_coverage");
		if(min_allele_coverage < ctr) {


			c->min_allele_coverage[c->length-1] = min_allele_coverage; 
		}
		else {
			c->min_allele_coverage[c->length-1] = ctr; 
		}
		fprintf(stderr, "coverage: %d\tminimum allele coverage: %d\n", cur_coverage, c->min_allele_coverage[c->length-1]);
		c->max_coverage = cur_coverage;
		cur_coverage++;
	}
	fprintf(stderr, "coverage: >%d\tminimum allele coverage: %d\n", c->max_coverage, min_allele_coverage);

	for(i=0;i<bin_coeff_rows;i++) {
		free(bin_coeff[i]);
	}
	free(bin_coeff);

	return c;
}

int32_t cov_cutoffs_get(cov_cutoffs_t *c, int32_t coverage)
{
	if(coverage < 0) { 
		return 0;
	}
	if(c->max_coverage < coverage) {
		return c->min_allele_coverage[c->max_coverage];
	}
	else {
		return c->min_allele_coverage[coverage];
	}
}

void cov_cutoffs_free(cov_cutoffs_t *c)
{
	free(c->min_allele_coverage);
	free(c);
}
