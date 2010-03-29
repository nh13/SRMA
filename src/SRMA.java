/*
 * LICENSE to be determined
 */
package srma;

import srma.Align;

import net.sf.samtools.*;
import net.sf.samtools.util.*;
import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.*;

import java.io.*;
import java.util.*;

public class SRMA extends CommandLineProgram { 

    @Usage public final String USAGE = getStandardUsagePreamble() + "Short read micro assembler.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file.")
        public File INPUT=null;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output SAM or BAM file.", optional=true)
        public File OUTPUT=null;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference FASTA file.")
        public File REFERENCE=null;
    @Option(doc="The alignment offset.", optional=true)
        public int OFFSET=20;
    @Option(doc="The minimum haploid coverage for the consensus.", optional=true)
        public int COVERAGE=0;
    @Option(doc="The file containing ranges to examine.", optional=true)
        public File RANGES=null;
    @Option(doc="A range to examine.", optional=true)
        public String RANGE=null;

    private List<ReferenceSequence> referenceSequences = null;
    private LinkedList<SAMRecord> toProcessSAMRecordList = null;
    private LinkedList<Node> toProcesSAMRecordNodeList = null;
    private PriorityQueue<SAMRecord> toOutputSAMRecordPriorityQueue = null;
    private SAMFileReader in = null;
    private SAMFileHeader header = null;
    private SAMFileWriter out = null;
    private Graph graph = null;
    private CloseableIterator<SAMRecord> recordIter = null;

    // for RANGES
    private boolean useRanges = false;
    private Ranges ranges = null;
    private Iterator<Range> rangeIterator = null;
    // for outputting within RANGES
    private int outputRangesIndex = -1;
    private Range outputRange = null;

    public static void main(final String[] args) {
        new SRMA().instanceMain(args);
    }
    /*
     * Current assumptions:
     * - can fit entire partial order graph in memory
     * */
    protected int doWork() 
    {
        int ctr=0;
        int prevReferenceIndex=-1, prevAlignmentStart=-1;

        try { 

            referenceSequences = new ArrayList();

            // Check input files
            IoUtil.assertFileIsReadable(INPUT);
            IoUtil.assertFileIsReadable(REFERENCE);

            // Initialize basic input/output files
            this.toProcessSAMRecordList = new LinkedList<SAMRecord>();
            this.toProcesSAMRecordNodeList = new LinkedList<Node>();
            this.toOutputSAMRecordPriorityQueue = new PriorityQueue(40, new SAMRecordCoordinateComparator()); 
            this.in = new SAMFileReader(INPUT, true);
            this.header = this.in.getFileHeader();
            if(null == OUTPUT) { // to STDOUT as a SAM
                this.out = new SAMFileWriterFactory().makeSAMWriter(this.header, true, System.out);
            }
            else { // to BAM file
                this.out = new SAMFileWriterFactory().makeSAMOrBAMWriter(this.header, true, OUTPUT);
            }

            // Get references
            this.getReferences(REFERENCE);

            // Get ranges
            if(null == RANGES && null == RANGE) {
                this.useRanges = false;
                // initialize SAM iter
                this.recordIter = this.in.iterator();
            }
            else if(null != RANGES && null != RANGE) {
                throw new Exception("RANGES and RANGE were both specified.\n");
            }
            else {
                this.useRanges = true;
                if(null != RANGES) {
                    IoUtil.assertFileIsReadable(RANGES);
                    this.ranges = new Ranges(RANGES, this.referenceSequences);
                }
                else {
                    this.ranges = new Ranges(RANGE, this.referenceSequences);
                }

                // initialize this.recordIter
                this.rangeIterator = ranges.iterator();
                if(!this.rangeIterator.hasNext()) {
                    return 0;
                }

                // initialize SAM iter
                Range r = this.rangeIterator.next();
                this.recordIter = this.in.query(this.referenceSequences.get(r.referenceIndex).getName(),
                        r.startPosition,
                        r.endPosition,
                        false);

                // initialize output ranges
                this.outputRangesIndex = 0;
                this.outputRange = this.ranges.get(this.outputRangesIndex);
            }

            // Initialize graph
            this.graph = new Graph(this.header, this.referenceSequences);

            SAMRecord rec = this.getNextSAMRecord();
            while(null != rec) {
                if(!rec.getReadUnmappedFlag()) { // only mapped reads
                    Node recNode = null;

                    // Make sure that it is sorted
                    if(rec.getReferenceIndex() < prevReferenceIndex || (rec.getReferenceIndex() == prevReferenceIndex && rec.getAlignmentStart() < prevAlignmentStart)) {
                        throw new Exception("SAM/BAM file is not co-ordinate sorted.");
                    }
                    prevReferenceIndex = rec.getReferenceIndex();
                    prevAlignmentStart = rec.getAlignmentStart();

                    // Add only if it is from the same contig
                    if(this.graph.contig != rec.getReferenceIndex()+1) {
                        // Process the rest of the reads
                        ctr = this.processList(ctr, false, false);
                    }

                    // Add to the graph 
                    try {
                        recNode = this.graph.addSAMRecord(rec);
                    } catch (Graph.GraphException e) {
                        if(Graph.GraphException.NOT_IMPLEMENTED != e.type) {
                            throw e;
                        }
                    }
                    if(this.useRanges) {
                        // Partition by the alignment start
                        if(this.recordAlignmentStartContained(rec)) {
                            this.toProcessSAMRecordList.add(rec);
                            this.toProcesSAMRecordNodeList.add(recNode);
                        }
                    }
                    else {
                        this.toProcessSAMRecordList.add(rec);
                        this.toProcesSAMRecordNodeList.add(recNode);
                    }

                    // Process the available reads
                    ctr = this.processList(ctr, true, false);
                }

                // get new record
                rec = this.getNextSAMRecord();
            }
            // Process the rest of the reads
            ctr = this.processList(ctr, true, true);


            // Close input/output files
            this.in.close();
            this.out.close();

            // Newline to end it all
            System.err.println("");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        return 0;
    }

    private void getReferences(File REFERENCE)
    {
        ReferenceSequenceFileFactory referenceSequenceFileFactory = new ReferenceSequenceFileFactory();
        ReferenceSequenceFile referenceSequenceFile = referenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
        ReferenceSequence referenceSequence = null;
        do {
            referenceSequence = referenceSequenceFile.nextSequence();
            if(null != referenceSequence) {
                this.referenceSequences.add(referenceSequence);
            }
        } while(null != referenceSequence);
    }

    private SAMRecord getNextSAMRecord()
    {
        if(this.recordIter.hasNext()) {
            return this.recordIter.next();
        }
        else if(this.useRanges) {
            do {
                if(this.rangeIterator.hasNext()) {
                    // close previous iterator
                    this.recordIter.close();
                    // get new range
                    Range r = this.rangeIterator.next();
                    // seek in the SAM file
                    this.recordIter = this.in.query(this.referenceSequences.get(r.referenceIndex).getName(),
                            r.startPosition,
                            r.endPosition,
                            false);
                }
                else {
                    this.recordIter.close();
                    return null;
                }
            } while(false == this.recordIter.hasNext());

            return this.recordIter.next();
        }
        else {
            this.recordIter.close();
            return null;
        }
    }

    private int processList(int ctr, boolean prune, boolean finish)
        throws Exception
    {
        SAMRecord curSAMRecord = null;

        // Process alignments
        while(0 < this.toProcessSAMRecordList.size() && ((!finish && this.toProcessSAMRecordList.getFirst().getAlignmentEnd() + this.OFFSET < this.toProcessSAMRecordList.getLast().getAlignmentStart()) || finish)) {
            curSAMRecord = this.toProcessSAMRecordList.removeFirst();
            Node curSAMRecordNode = this.toProcesSAMRecordNodeList.removeFirst();
            ctr++;
            System.err.print("\rctr:" + ctr + " AL:" + curSAMRecord.getAlignmentStart() + ":" + curSAMRecord.getAlignmentEnd() + ":" + curSAMRecord.toString());
            if(prune) {
                this.graph.prune(curSAMRecord.getReferenceIndex(), curSAMRecord.getAlignmentStart(), this.OFFSET); 
            }

            // Align - this will overwrite/change the alignment
            curSAMRecord = Align.align(this.graph, curSAMRecord, curSAMRecordNode, this.referenceSequences, OFFSET, COVERAGE);
            // Add to a heap/priority-queue to assure output is sorted
            this.toOutputSAMRecordPriorityQueue.add(curSAMRecord);
        }

        // Output alignments
        while(0 < this.toOutputSAMRecordPriorityQueue.size()) {
            curSAMRecord = this.toOutputSAMRecordPriorityQueue.peek();
            if(finish || curSAMRecord.getAlignmentStart() < graph.position_start) { // other alignments will not be less than
                this.out.addAlignment(this.toOutputSAMRecordPriorityQueue.poll());
            }
            else { // other alignments could be less than
                break;
            }
        }

        return ctr;
    }

    private boolean recordAlignmentStartContained(SAMRecord rec) 
    {
        int recAlignmentStart = -1;

        if(this.ranges.size() <= this.outputRangesIndex) { // no more ranges
            return false;
        }

        recAlignmentStart = rec.getAlignmentStart();
        while(this.outputRange.endPosition < recAlignmentStart) { // find a new range
            this.outputRangesIndex++;
            if(this.ranges.size() <= this.outputRangesIndex) {
                return false;
            }
            this.outputRange = this.ranges.get(this.outputRangesIndex);
        }
        if(recAlignmentStart < this.outputRange.startPosition) {
            // not within range
            return false;
        }
        else {
            // must be within range
            return true;
        }
    }
}
