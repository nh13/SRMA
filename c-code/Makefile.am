bin_PROGRAMS = srma

srma_SOURCES = \
			   src/srma_error.c src/srma_error.h \
			   src/srma_alloc.c src/srma_alloc.h \
			   src/srma_util.c src/srma_util.h \
			   src/srma_sam_io.c src/srma_sam_io.h \
			   src/srma_faidx.c src/srma_faidx.h  \
			   src/sw_node.c src/sw_node.h \
			   src/sw_heap.c src/sw_heap.h \
			   src/sw_align.c src/sw_align.h \
			   src/node.c src/node.h \
			   src/graph.c src/graph.h \
			   src/main.c src/main.h

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
