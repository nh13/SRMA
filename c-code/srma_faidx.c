#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <ctype.h>

#include "samtools/faidx.c" // since we want access to variables within faidx.c
#include "samtools/khash.h"
#include "srma_error.h"
#include "srma_alloc.h"
#include "srma_faidx.h"


// modified from samtools/faidx.c
static inline int srma_ranges_parse_range(srma_ranges_t *r, const faidx_t *fai, char *str, srma_range_t *range, int32_t offset)
{
	int32_t i, j, k, l;
	char *s = NULL, *p = NULL;
	faidx1_t val;
	khiter_t iter;
	khash_t(s) *h = NULL;
	int32_t tid, beg, end;

	beg = end = -1;
	h = fai->hash;
	l = strlen(str);
	p = s = srma_malloc(sizeof(char)*(l+1), __func__, s);
	// remove commas (',')
	for (i = k = 0; i != l; ++i) 
		if (str[i] != ',' && !isspace(str[i])) s[k++] = str[i];
	s[k] = '\0';
	// skip to colon
	for (i = 0; i != k; ++i) if (s[i] == ':') break;
	s[i] = '\0';
	// get the reference id
	iter = kh_get(s, h, s); // get the reference id
	if(iter == kh_end(h)) { // not found
		fprintf(stderr, "Reference sequence [%s] not found in the FASTA index", s);
		free(s);
		return 0;
	}
	tid = -1;
	for(j=0;j<fai->n;j++) {
		if(0 == strcmp(s, fai->name[j])) {
			tid = j;
			break;
		}
	}
	if(tid < 0) { // we should fin it
		srma_error(__func__, "tid < 0", Exit, OutOfRange);
	}
	val = kh_value(h, iter);
	if (i == k) { /* dump the whole sequence */
		beg = 0; end = val.len-1;
	} else {
		for (p = s + i + 1; i != k; ++i) if (s[i] == '-') break;
		// TODO: make sure atoi works here
		beg = atoi(p) -1;
		if (i < k) {
			p = s + i + 1;
			// TODO: make sure atoi works here
			end = atoi(p) - 1;
		} else end = val.len-1;
		// expand range based on offset
		beg -= offset;
		end += offset;
	}
	// free 
	free(s);
	// bounds check
	if(beg < 0) beg = 0;
	else if(beg >= val.len) beg = val.len-1;
	if(end < 0) end = 0;
	else if(end >= val.len) end = val.len-1;
	if (beg > end) return 0; // parse error

	srma_ranges_add(r, tid, beg, end);

	return 1;
}

static void srma_ranges_init_range(srma_ranges_t *r, const faidx_t *fai, char *range, int32_t offset)
{
	if(0 == srma_ranges_parse_range(r, fai, range, &r->ranges[0], offset)) {
		srma_error(__func__, range, Exit, OutOfRange);
	}
}

static void srma_ranges_init_ranges(srma_ranges_t *r, const faidx_t *fai, char *fn_ranges, int32_t offset)
{
	FILE *fp = NULL;
	char line[1024]="\0";

	if(0 == (fp = fopen(fn_ranges, "r"))) {
		srma_error(__func__, fn_ranges, Exit, OpenFileError);
	}

	while(NULL != fgets(line, 1024, fp)) {
		if(0 < strlen(line)) {
			line[strlen(line)-1]='\0';
		}
		if(0 == srma_ranges_parse_range(r, fai, line, &r->ranges[r->length-1], offset)) {
			srma_error(__func__, line, Exit, OutOfRange);
		}
	}

	fclose(fp);
}

srma_ranges_t *srma_ranges_init(const faidx_t *fai, char *fn_ranges, char *range, int32_t offset)
{
	srma_ranges_t *r = NULL;
	
	r = srma_malloc(sizeof(srma_ranges_t), __func__, "r");
	r->i=0;
	r->length=0;
	r->ranges = NULL;
	if(NULL != fn_ranges) {
		srma_ranges_init_ranges(r, fai, fn_ranges, offset);
	}
	else if(NULL != range) {
		srma_ranges_init_range(r, fai, range, offset);
	}
	else {
		return r;
	}
	if(0 == r->length) {
		srma_error(__func__, "Did not add any ranges", Exit, OutOfRange);
	}
	// TODO: check increasing order, and overlaps
	return r;
}

void srma_ranges_add(srma_ranges_t *r, int32_t tid, int32_t beg, int32_t end)
{
	int32_t i, j;
	// find insertion point
	// could binary search, but not so important to be efficient here
	for(i=0;i<r->length;i++) {
		// is next greater than
		if(tid < r->ranges[i].tid || (tid == r->ranges[i].tid && beg < r->ranges[i].beg)) {
			break;
		}
	}
	// insert
	r->length++;
	r->ranges = srma_realloc(r->ranges, sizeof(srma_range_t)*r->length, __func__, "r->ranges");
	for(j=r->length-1;i<j;j--) { // shift down
		r->ranges[j].tid = r->ranges[j-1].tid;
		r->ranges[j].beg = r->ranges[j-1].beg;
		r->ranges[j].end = r->ranges[j-1].end;
	}
	r->ranges[i].tid = tid;
	r->ranges[i].beg = beg;
	r->ranges[i].end = end;
	// check for overlapping ranges and merge
	for(i=r->length-1,j=0;0<i;i--) {
		if(r->ranges[i].tid == r->ranges[i-1].tid &&
				r->ranges[i].beg <= r->ranges[i-1].end) {  // overlapping
			// shift down
			if(r->ranges[i-1].end < r->ranges[i].end) { // not wholly contained
				r->ranges[i-1].end = r->ranges[i].end;
			}
			// count how many overlap
			j++;
		}
	}
	// reallocate memory, removing the end j items
	if(0 < j) { 
		r->length-=j;
		assert(0 < r->length);
		r->ranges = srma_realloc(r->ranges, sizeof(srma_range_t)*r->length, __func__, "r->ranges");
	}
}

srma_ranges_t *srma_ranges_init_from_fai(faidx_t *fai)
{
	int32_t i;
	khiter_t iter;
	srma_ranges_t *r;
	faidx1_t val;
	
	r = srma_malloc(sizeof(srma_ranges_t), __func__, "r");
	r->i=0;
	r->length=0;
	r->ranges = NULL;

	for(i=0;i<fai->n;i++) {
		iter = kh_get(s, fai->hash, fai->name[i]);
		if(iter == kh_end(fai->hash)) srma_error(__func__, "iter == kh_end(fai->hash)", Exit, OutOfRange);
		val = kh_value(fai->hash, iter);
		srma_ranges_add(r, i, 0, val.len-1); // zero-based
	}

	return r;
}

inline srma_range_t *srma_ranges_peek(srma_ranges_t *r)
{
	if(NULL != r && r->i < r->length) {
		return &r->ranges[r->i];
	}
	else {
		return NULL;
	}
}

inline srma_range_t *srma_ranges_poll(srma_ranges_t *r)
{
	if(NULL != r && r->i+1 < r->length) { // more left ?
		r->i++;
		return &r->ranges[r->i];
	}
	else {
		return NULL;
	}
}

void srma_ranges_free(srma_ranges_t *r)
{
	if(NULL != r) {
		free(r->ranges);
	}
	free(r);
}

// want control of malloc and bound checking
void srma_fai_fetch(const faidx_t *fai, ref_t *ref, int32_t tid, int32_t beg, int32_t end)
{
	int l;
	char c;
	khiter_t iter;
	faidx1_t val;
	char *seq=NULL;

	// Adjust position
	iter = kh_get(s, fai->hash, fai->name[tid]);
	if(iter == kh_end(fai->hash)) srma_error(__func__, "iter == kh_end(fai->hash)", Exit, OutOfRange);
	val = kh_value(fai->hash, iter);
	if(beg == -1) beg = 0;
	if(end == -1) end = val.len-1;

	// Validate inputs
	if(end < beg) srma_error(__func__, "end < beg", Exit, OutOfRange);
	if(beg < 0) srma_error(__func__, "beg < 0", Exit, OutOfRange);
	else if(val.len <= beg) srma_error(__func__, "val.len <= beg", Exit, OutOfRange);
	if(end < 0) srma_error(__func__, "end < 0", Exit, OutOfRange);
	else if(val.len <= end) srma_error(__func__, "val.len <= end", Exit, OutOfRange);
	
	free(ref->ref);
	ref->ref=NULL;
	ref->tid = tid;
	ref->beg = beg;
	ref->end = end;

	// Now retrieve the sequence 
	l = 0;
	seq = srma_malloc(sizeof(char)*(end - beg + 2), __func__, "seq");
	razf_seek(fai->rz, val.offset + beg / val.line_blen * val.line_len + beg % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < end - beg + 1)
		if (isgraph(c)) seq[l++] = c;
	seq[l] = '\0';
	if(l != end - beg + 1) srma_error(__func__, "len != end - beg + 1", Exit, OutOfRange);

	ref->ref = seq;
}

inline char *srma_fai_name(faidx_t *fai, int tid)
{
	return fai->name[tid];
}
