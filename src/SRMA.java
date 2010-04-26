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

//import java.lang.Runtime;
import java.io.*;
import java.util.*;
import java.lang.Math;

/* Documentation:
 * */

public class SRMA extends CommandLineProgram { 

    public final String PROGRAM_VERSION="0.1.5";
    @Usage (programVersion=PROGRAM_VERSION)
        public final String USAGE = getStandardUsagePreamble() + "Short read micro re-aligner.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file.")
        public File INPUT=null;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output SAM or BAM file.", optional=true)
        public File OUTPUT=null;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference FASTA file.")
        public File REFERENCE=null;
    @Option(doc="The alignment offset.", optional=true)
        public int OFFSET=20;
    @Option(doc="The minimum mapping quality.", optional=true)
        public int MIN_MAPQ=0;
    @Option(doc="The minimum allele probability conditioned on coverage (for the binomial quantile).", optional=true)
        public double MINIMUM_ALLELE_PROBABILITY=0.1;
    @Option(doc="The minimum haploid coverage for the consensus.", optional=true)
        public int MINIMUM_ALLELE_COVERAGE=3;
    @Option(doc="The file containing ranges to examine.", optional=true)
        public File RANGES=null;
    @Option(doc="A range to examine.", optional=true)
        public String RANGE=null;
    @Option(doc="Correct bases.", optional=true)
        public boolean CORRECT_BASES=false;
    @Option(doc="Use sequence qualities", optional=true)
        public boolean USE_SEQUENCE_QUALITIES=true;
    @Option(doc="Whether to suppress job-progress info on System.err", optional=true)
        public boolean QUIET_STDERR=false;
    @Option(doc="The maximum number of nodes on the heap before re-alignment is ignored", optional=true)
        public int MAX_HEAP_SIZE = 8192;

    private long startTime;
    private long endTime;


    private final static int SRMA_OUTPUT_CTR = 100;
    private int maxOutputStringLength = 0;
    private String maxOutputString = null;

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
    private AlleleCoverageCutoffs alleleCoverageCutoffs = null;

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

    protected String[] customCommandLineValidation()
    {
        return super.customCommandLineValidation();
    }


    /*
     * Current assumptions:
     * - can fit entire partial order graph in memory
     * */
    protected int doWork() 
    {
        int ctr=0;
        int prevReferenceIndex=-1, prevAlignmentStart=-1;

        // initialize
        this.alleleCoverageCutoffs = new AlleleCoverageCutoffs(MINIMUM_ALLELE_COVERAGE, MINIMUM_ALLELE_PROBABILITY, QUIET_STDERR);

        try { 
            this.startTime = System.nanoTime();

            // Check input files
            IoUtil.assertFileIsReadable(INPUT);
            IoUtil.assertFileIsReadable(REFERENCE);

            // Initialize basic input/output files
            this.toProcessSAMRecordList = new LinkedList<SAMRecord>();
            this.toProcesSAMRecordNodeList = new LinkedList<Node>();
            this.toOutputSAMRecordPriorityQueue = new PriorityQueue<SAMRecord>(40, new SAMRecordCoordinateComparator()); 
            this.in = new SAMFileReader(INPUT, true);
            this.header = this.in.getFileHeader();

            // Add SRMA to the header
            SAMProgramRecord programRecord = this.header.getProgramRecord("srma");
            String programVersion = new String(PROGRAM_VERSION);
            if(null == programRecord) { // create a new one
                programRecord = new SAMProgramRecord("srma");
                programRecord.setProgramVersion(programVersion);
                this.header.addProgramRecord(programRecord);
            }
            else if(0 != programVersion.compareTo(programRecord.getProgramVersion())) { // new version, but srma exists
                programVersion = new String("srma-" + PROGRAM_VERSION); // append "srma-" so we know it was srma
                programRecord = this.header.createProgramRecord();
                programRecord.setProgramVersion(programVersion);
            }

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
                if(rec.getReadUnmappedFlag()) { 
                    // TODO
                    // Print this out somehow in some order somewhere
                }
                else if(rec.getMappingQuality() < MIN_MAPQ) {
                    // TODO
                    // Print this out somehow in some order somewhere
                }
                else {
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
                        ctr = this.processList(programRecord, ctr, false, true);

                        // Get new reference sequence
                        while(null != this.referenceSequence && this.referenceSequence.getContigIndex() < rec.getReferenceIndex()) {
                            this.referenceSequence = this.referenceSequenceFile.nextSequence();
                        }
                        if(null == this.referenceSequence) {
                            throw new Exception("Premature EOF in the reference sequence");
                        }
                        else if(this.referenceSequence.getContigIndex() != rec.getReferenceIndex()) {
                            throw new Exception("Could not find the reference sequence");
                        }
                    }

                    // Add to the graph 
                    recNode = this.graph.addSAMRecord(rec, this.referenceSequence);

                    if(null != recNode) { // successfully added
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
                    }

                    // Process the available reads
                    ctr = this.processList(programRecord, ctr, true, false);
                }

                // get new record
                rec = this.getNextSAMRecord();
            }
            // Process the rest of the reads
            ctr = this.processList(programRecord, ctr, true, true);


            // Close input/output files
            this.in.close();
            this.out.close();

            this.endTime = System.nanoTime();

            // to end it all
            if(!QUIET_STDERR) {
                System.err.println("");
                System.err.println("SRMA complete");
                // Memory
                double totalMemory = (double)Runtime.getRuntime().totalMemory();
                double totalMemoryLog2 = Math.log(totalMemory) / Math.log(2.0);
                if(totalMemoryLog2 < 10) {
                    System.err.println("Total memory usage: " + (int)totalMemory + "B");
                } 
                else if(totalMemoryLog2 < 20) {
                    System.err.println("Total memory usage: " + (Math.round(100 * totalMemory / Math.pow(2, 10)) / 100) + "KB");
                }
                else if(totalMemoryLog2 < 30) {
                    System.err.println("Total memory usage: " + (Math.round(100 * totalMemory / Math.pow(2, 20)) / 100) + "MB");
                }
                else {
                    System.err.println("Total memory usage: " + (Math.round(100 * totalMemory / Math.pow(2, 30)) / 100) + "GB");
                }
                // Run time
                long seconds = (this.endTime - this.startTime) / 1000000000;
                long hours = seconds / 3600; seconds -= hours * 3600; 
                long minutes = seconds / 60; seconds -= minutes* 60; 
                System.err.println("Total execution time: " + hours + "h : " + minutes + "m : " + seconds + "s");
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        // this is annoying
        QUIET = true;

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
            
        if(QUIET_STDERR) {
            return;
        }
        else if(0 == (ctr % SRMA_OUTPUT_CTR) || ctr < 0) {
            // TODO: enforce column width ?
            String outputString = new String("Records processsed: " + ctr + " (" + rec.getReferenceName() + ":" + rec.getAlignmentStart() + "-" + rec.getAlignmentEnd() + ")");
            int outputStringLength = outputString.length();
            if(this.maxOutputStringLength < outputStringLength) {
                this.maxOutputStringLength = outputStringLength;
                int i;
                maxOutputString = new String("");
                for(i=outputStringLength;i < this.maxOutputStringLength;i++) { // pad with blanks
                    this.maxOutputString += " ";
                }
            }
            System.err.print(this.maxOutputString + "\r" + outputString);
        }
    }

    private int processList(SAMProgramRecord programRecord, int ctr, boolean prune, boolean finish)
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
                    programRecord,
                    OFFSET, 
                    this.alleleCoverageCutoffs,
                    CORRECT_BASES,
                    USE_SEQUENCE_QUALITIES,
                    MAX_HEAP_SIZE);
            // Add to a heap/priority-queue to assure output is sorted
            this.toOutputSAMRecordPriorityQueue.add(curSAMRecord);
        }
        if(finish && null != curSAMRecord) {
            this.outputProgress(curSAMRecord, ctr);
        }
        curSAMRecord = null;


        // Output alignments
        while(0 < this.toOutputSAMRecordPriorityQueue.size()) {
            curSAMRecord = this.toOutputSAMRecordPriorityQueue.peek();
            // alignment could have moved (+OFFSET), with another moving (-OFFSET) 
            if(finish || curSAMRecord.getAlignmentStart() + 2*OFFSET < graph.position_start) { // other alignments will not be less than
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
