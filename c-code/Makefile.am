bin_PROGRAMS = srma

srma_SOURCES = \
			   srma_error.c srma_error.h \
			   srma_alloc.c srma_alloc.h \
			   srma_util.c srma_util.h \
			   srma_sam_io.c srma_sam_io.h \
			   srma_faidx.c srma_faidx.h  \
			   sw_node.c sw_node.h \
			   sw_heap.c sw_heap.h \
			   sw_align.c sw_align.h \
			   node.c node.h \
			   graph.c graph.h \
			   main.c main.h

nodist_srma_SOURCES = \
			   samtools/bam_import.c samtools/bam_import.h \
			   samtools/bam_index.c samtools/bam_index.h \
			   samtools/bam_aux.c samtools/bam_aux.h \
			   samtools/bam_pileup.c samtools/bam_pileup.h \
			   samtools/bam.c samtools/bam.h \
			   samtools/sam.c samtools/sam.h \
			   samtools/bgzf.c    samtools/bgzf.h \
			   samtools/razf.c samtools/razf.h \
			   samtools/kstring.c samtools/kstring.h \
			   samtools/sam_header.c   samtools/sam_header.h 

srma_LDADD = -lz

EXTRA_DIST = LICENSE \
			 autogen.sh 

docdir = ${datadir}/doc/${PACKAGE}
dist_doc_DATA = LICENSE 
