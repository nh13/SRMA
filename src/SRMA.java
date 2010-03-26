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

    public static final String OFFSET_SHORT_NAME = "O"; 
    public static final String COVERAGE_SHORT_NAME = "C";
    public static final String RANGES_SHORT_NAME = "X";

    @Usage public final String USAGE = getStandardUsagePreamble() + "Prints a SAM or BAM file to the screen.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file to view.")
        public File INPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference FASTA file.")
        public File REFERENCE;
    @Option(shortName=OFFSET_SHORT_NAME, doc="The alignment offset.", optional=true)
        public int OFFSET=0;
    @Option(shortName=COVERAGE_SHORT_NAME, doc="The minimum haploid coverage for the consensus.", optional=true)
        public int COVERAGE=1;
    @Option(shortName=RANGES_SHORT_NAME, doc="The file containing ranges to examine.", optional=true)
        public File RANGES=null;

    private List<ReferenceSequence> referenceSequences = null;
    private LinkedList<SAMRecord> samRecordList = null;
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

        referenceSequences = new ArrayList();

        // Check input files
        IoUtil.assertFileIsReadable(INPUT);
        IoUtil.assertFileIsReadable(REFERENCE);

        // Initialize basic input/output files
        this.samRecordList = new LinkedList<SAMRecord>();
        this.in = new SAMFileReader(INPUT, true);
        this.header = this.in.getFileHeader();
        this.out = new SAMFileWriterFactory().makeSAMWriter(this.header, true, System.out);

        // Get references
        this.getReferences(REFERENCE);

        // Get ranges
        if(null == RANGES) {
            this.useRanges = false;
            // initialize SAM iter
            this.recordIter = this.in.iterator();
        }
        else {
            this.useRanges = true;
            // check file
            IoUtil.assertFileIsReadable(RANGES);

            // initialize this.recordIter
            this.ranges = new Ranges(RANGES, this.referenceSequences);
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

        System.err.println("");
        try { // Go through each SAM record
            SAMRecord rec = this.getNextSAMRecord();
            while(null != rec) {
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
                    this.graph.addSAMRecord(rec);
                } catch (Graph.GraphException e) {
                    if(Graph.GraphException.NOT_IMPLEMENTED != e.type) {
                        throw e;
                    }
                }
                if(this.useRanges) {
                    // Partition by the alignment start
                    if(this.recordAlignmentStartContained(rec)) {
                        System.err.println("ADDING " + rec.getAlignmentStart());
                        this.samRecordList.add(rec);
                    }
                }
                else {
                    this.samRecordList.add(rec);
                }

                // Process the available reads
                ctr = this.processList(ctr, true, false);

                // get new record
                rec = this.getNextSAMRecord();
            }
            // Process the rest of the reads
            ctr = this.processList(ctr, true, true);

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        // Close input/output files
        this.in.close();
        this.out.close();

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
        while(0 < this.samRecordList.size() && 
                ((!finish && this.samRecordList.getFirst().getAlignmentEnd() + this.OFFSET < this.samRecordList.getLast().getAlignmentStart()) || finish)) {
            SAMRecord curSAMRecord = this.samRecordList.removeFirst();
            ctr++;
            System.err.print("\rctr:" + ctr + " AL:" + curSAMRecord.getAlignmentStart() + ":" + curSAMRecord.getAlignmentEnd() + ":" + curSAMRecord.toString());
            if(prune) {
                this.graph.prune(curSAMRecord.getReferenceIndex(), curSAMRecord.getAlignmentStart(), this.OFFSET); 
            }
            this.out.addAlignment(Align.Align(this.graph, curSAMRecord, OFFSET, COVERAGE));
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
