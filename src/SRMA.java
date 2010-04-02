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

/* Documentation:
 * - allow if MINIMUM_ALLELE_FREQUENCY or MINIMUM_ALLELE_COVERAGE
 * - requires a ".dict" file for the FASTA
 * */

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
    @Option(doc="The minimum allele frequency for the conensus.", optional=true)
        public double MINIMUM_ALLELE_FREQUENCY=0.0;
    @Option(doc="The minimum haploid coverage for the consensus.", optional=true)
        public int MINIMUM_ALLELE_COVERAGE=0;
    @Option(doc="The file containing ranges to examine.", optional=true)
        public File RANGES=null;
    @Option(doc="A range to examine.", optional=true)
        public String RANGE=null;

    private final static int SRMA_OUTPUT_CTR = 100;
    private int maxOutputStringLength = 0;

    ReferenceSequenceFile referenceSequenceFile = null; 
    private ReferenceSequence referenceSequence = null;
    private SAMSequenceDictionary referenceDictionary = null;

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
    // for inputting within RANGES
    private Ranges inputRanges = null;
    private Iterator<Range> inputRangesIterator = null;
    private Range inputRange = null;
    // for outputting within RANGES
    private Ranges outputRanges = null;
    private Iterator<Range> outputRangesIterator = null;
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
            ReferenceSequenceFileFactory referenceSequenceFileFactory = new ReferenceSequenceFileFactory();
            this.referenceSequenceFile = referenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
            this.referenceDictionary = this.referenceSequenceFile.getSequenceDictionary();
            if(null == this.referenceDictionary) {
                // Try manually
                String dictionaryName = new String(REFERENCE.getAbsolutePath());
                dictionaryName += ".dict";
                System.err.println("DICT: " + dictionaryName);
                final File dictionary = new File(dictionaryName);
                if (dictionary.exists()) {
                    IoUtil.assertFileIsReadable(dictionary);
                    final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
                    final SAMFileHeader header = codec.decode(new AsciiLineReader(new FileInputStream(dictionary)),
                            dictionary.toString());
                    if (header.getSequenceDictionary() != null && header.getSequenceDictionary().size() > 0) {
                        this.referenceDictionary = header.getSequenceDictionary();
                    }
                }
                else {
                    throw new Exception("Could not open sequence dictionary file: " + dictionaryName);
                }
            }

            // Get ranges
            if(null == RANGES && null == RANGE) {
                this.useRanges = false;
                // initialize SAM iter
                this.recordIter = this.in.iterator();
                this.referenceSequence = this.referenceSequenceFile.nextSequence();
            }
            else if(null != RANGES && null != RANGE) {
                throw new Exception("RANGES and RANGE were both specified.\n");
            }
            else {
                this.useRanges = true;
                if(null != RANGES) {
                    IoUtil.assertFileIsReadable(RANGES);
                    this.inputRanges = new Ranges(RANGES, this.referenceDictionary, OFFSET);
                    this.outputRanges = new Ranges(RANGES, this.referenceDictionary);
                }
                else {
                    this.inputRanges = new Ranges(RANGE, this.referenceDictionary, OFFSET);
                    this.outputRanges = new Ranges(RANGE, this.referenceDictionary);
                }

                this.inputRangesIterator = this.inputRanges.iterator();
                this.outputRangesIterator = this.outputRanges.iterator();
                if(!this.inputRangesIterator.hasNext()) {
                    return 0;
                }

                this.inputRange = this.inputRangesIterator.next();

                do {
                    this.referenceSequence = this.referenceSequenceFile.nextSequence();
                } while(null != this.referenceSequence && this.referenceSequence.getContigIndex() < this.inputRange.referenceIndex);
                if(null == this.referenceSequence) {
                    throw new Exception("Premature EOF in the reference sequence");
                }
                else if(this.referenceSequence.getContigIndex() != this.inputRange.referenceIndex) {
                    throw new Exception("Could not find the reference sequence");
                }

                this.recordIter = this.in.query(this.referenceDictionary.getSequence(this.inputRange.referenceIndex).getSequenceName(),
                        this.inputRange.startPosition,
                        this.inputRange.endPosition,
                        false);
                this.outputRange = this.outputRangesIterator.next();
            }

            // Initialize graph
            this.graph = new Graph(this.header);

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
                        recNode = this.graph.addSAMRecord(rec, this.referenceSequence);
                    } catch (Graph.GraphException e) {
                        if(Graph.GraphException.NOT_IMPLEMENTED != e.type) {
                            throw e;
                        }
                    }
                    if(this.useRanges) {
                        // Partition by the alignment start
                        if(this.recordAlignmentStartContained(rec)) { // only add if it will be outputted
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
                else {
                    // TODO
                    // Print this out somehow in some order somewhere
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

    private SAMRecord getNextSAMRecord()
        throws Exception
    {
        if(this.recordIter.hasNext()) {
            return this.recordIter.next();
        }
        else if(this.useRanges) {
            do {
                if(this.inputRangesIterator.hasNext()) {
                    // close previous iterator
                    this.recordIter.close();
                    // get new range
                    this.inputRange = this.inputRangesIterator.next();
                    // seek in the SAM file
                    this.recordIter = this.in.query(this.referenceDictionary.getSequence(this.inputRange.referenceIndex).getSequenceName(),
                            this.inputRange.startPosition,
                            this.inputRange.endPosition,
                            false);
                }
                else {
                    this.recordIter.close();
                    return null;
                }
            } while(false == this.recordIter.hasNext());
            while(null != this.referenceSequence && this.referenceSequence.getContigIndex() < this.inputRange.referenceIndex) {
                this.referenceSequence = this.referenceSequenceFile.nextSequence();
            }
            if(null == this.referenceSequence) {
                throw new Exception("Premature EOF in the reference sequence");
            }
            else if(this.referenceSequence.getContigIndex() != this.inputRange.referenceIndex) {
                throw new Exception("Could not find the reference sequence");
            }

            return this.recordIter.next();
        }
        else {
            this.recordIter.close();
            return null;
        }
    }

    private void outputProgress(SAMRecord rec, int ctr)
    {
        if(0 == (ctr % SRMA_OUTPUT_CTR)) {
            // TODO: enforce column width ?
            String outputString = new String("ctr:" + ctr + " AL:" + rec.getAlignmentStart() + ":" + rec.getAlignmentEnd() + ":" + rec.toString());
            int outputStringLength = outputString.length();
            if(this.maxOutputStringLength < outputStringLength) {
                this.maxOutputStringLength = outputStringLength;
            }
            System.err.print("\r" + outputString);
            int i;
            for(i=outputStringLength;i < this.maxOutputStringLength;i++) { // pad with blanks
                System.err.print(" ");
            }
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
            if(prune) {
                this.graph.prune(curSAMRecord.getReferenceIndex(), curSAMRecord.getAlignmentStart(), this.OFFSET); 
            }
            this.outputProgress(curSAMRecord, ctr);
            ctr++;

            // Align - this will overwrite/change the alignment
            curSAMRecord = Align.align(this.graph, 
                    curSAMRecord, 
                    curSAMRecordNode, 
                    this.referenceSequence,
                    OFFSET, 
                    MINIMUM_ALLELE_FREQUENCY,
                    MINIMUM_ALLELE_COVERAGE); 
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

        if(null == this.outputRange) { // no more ranges
            return false;
        }

        recAlignmentStart = rec.getAlignmentStart();
        while(this.outputRange.endPosition < recAlignmentStart) { // find a new range
            if(!this.outputRangesIterator.hasNext()) { // no more ranges
                this.outputRange = null;
                return false;
            }
            this.outputRange = this.outputRangesIterator.next();
        }
        if(recAlignmentStart < this.outputRange.startPosition) { // before range
            // not within range
            return false;
        }
        else {
            // must be within range
            return true;
        }
    }
}
