/*
 * LICENSE to be determined
 */
package srma;

import srma.Align;
import srma.ThreadPoolLinkedList;
import srma.AlignRecordComparator;
import srma.SAMRecordIO;

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

    public final String PROGRAM_VERSION="0.1.14";
    @Usage (programVersion=PROGRAM_VERSION)
        public final String USAGE = getStandardUsagePreamble() + "Short read micro re-aligner.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file.")
        public List<File> INPUT = new ArrayList<File>();
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output SAM or BAM file.", optional=true)
        public List<File> OUTPUT = new ArrayList<File>();
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
    @Option(doc="The maximum total coverage over a position to consider, otherwise ignore the re-alignment.", optional=true)
        public int MAXIMUM_TOTAL_COVERAGE=100;
    @Option(doc="The file containing ranges to examine.", optional=true)
        public File RANGES=null;
    @Option(doc="A range to examine.", optional=true)
        public String RANGE=null;
    @Option(doc="Correct bases.", optional=true)
        public boolean CORRECT_BASES=false;
    @Option(doc="Use sequence qualities", optional=true)
        public boolean USE_SEQUENCE_QUALITIES=true;
    @Option(doc="Only align read pairs within this range (ex. -1000:1000); read #1 is assumed to be 5' of read #2", optional=true)
        public String INSERT_SIZE_RANGE=null;
    @Option(doc="Whether to suppress job-progress info on System.err", optional=true)
        public boolean QUIET_STDERR=false;
    @Option(doc="The maximum number of nodes on the heap before re-alignment is ignored", optional=true)
        public int MAX_HEAP_SIZE = 8192;
    @Option(doc="The maximum number of SAM records in the queue before re-alignment", optional=true)
        public int MAX_QUEUE_SIZE = 65536;
    @Option(doc="The number of threads for parallel processing", optional=true)
        public int NUM_THREADS = 1;

    private long startTime;
    private long endTime;

    private final static int SRMA_OUTPUT_CTR = 100;
    private int maxOutputStringLength = 0;

    ReferenceSequenceFile referenceSequenceFile = null; 
    private ReferenceSequence referenceSequence = null;
    private SAMSequenceDictionary referenceDictionary = null;

    private ThreadPoolLinkedList toAddToGraphList = null;
    private ThreadPoolLinkedList toAlignList = null;
    private PriorityQueue<AlignRecord> toOutputQueue = null;

    private SAMRecordIO io = null;

    private Graph graph = null;
    private AlleleCoverageCutoffs alleleCoverageCutoffs = null;
    private int InsertSizeRangeLow, InsertSizeRangeHigh;

    // for RANGES
    // for inputting within RANGES
    private Ranges inputRanges = null;
    private Iterator<Range> inputRangesIterator = null;
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
        AlignRecord rec = null;

        try { 
            this.startTime = System.nanoTime();

            if(1 < this.NUM_THREADS) {
                System.err.println("** Warning: option NUM_THREADS currently may not increase performance significantly. **");
                System.err.println("**  Try running multiple processes with RANGE if the speed does not increase.      **");
            }

            IoUtil.assertFileIsReadable(REFERENCE);

            if(0 == this.OUTPUT.size() && !QUIET) {
                throw new Exception("Please use option 'QUIET' when outputting to stdout.");
            }

            // Get references
            this.referenceSequenceFile = new IndexedFastaSequenceFile(REFERENCE);
            if(!this.referenceSequenceFile.isIndexed()) {
                throw new Exception("Reference sequence file was not indexed.");
            }
            this.referenceDictionary = this.referenceSequenceFile.getSequenceDictionary();
            if(null == this.referenceDictionary) {
                throw new Exception("Could not find FASTA dictionary file.");
            }

            // Get ranges
            if(null != RANGES && null != RANGE) {
                throw new Exception("RANGES and RANGE were both specified.\n");
            }
            else if(null == RANGES && null == RANGE) {
                this.inputRanges = new Ranges(this.referenceDictionary);
                this.outputRanges = new Ranges(this.referenceDictionary);
                this.io = new SAMRecordIO(INPUT, OUTPUT, PROGRAM_VERSION, false, this.referenceDictionary);
            }
            else {
                if(null != RANGES) {
                    IoUtil.assertFileIsReadable(RANGES);
                    this.inputRanges = new Ranges(RANGES, this.referenceDictionary, OFFSET);
                    this.outputRanges = new Ranges(RANGES, this.referenceDictionary, 0);
                }
                else {
                    this.inputRanges = new Ranges(RANGE, this.referenceDictionary, OFFSET);
                    this.outputRanges = new Ranges(RANGE, this.referenceDictionary, 0);
                }
                this.io = new SAMRecordIO(INPUT, OUTPUT, PROGRAM_VERSION, true, this.referenceDictionary);
            }

            this.inputRangesIterator = this.inputRanges.iterator();
            this.outputRangesIterator = this.outputRanges.iterator();
            if(!this.inputRangesIterator.hasNext()) {
                return 0;
            }
            this.outputRange = this.outputRangesIterator.next();

            // init
            this.alleleCoverageCutoffs = new AlleleCoverageCutoffs(MINIMUM_ALLELE_COVERAGE, MINIMUM_ALLELE_PROBABILITY, QUIET_STDERR);
            this.toOutputQueue = new PriorityQueue<AlignRecord>(40, new AlignRecordComparator()); 
            this.InsertSizeRangeSet();

            while(this.inputRangesIterator.hasNext()) {
                Range inputRange = this.inputRangesIterator.next();

                int prevReferenceIndex=-1;
                int prevAlignmentStart=-1;

                this.toAddToGraphList = new ThreadPoolLinkedList();
                this.toAlignList = new ThreadPoolLinkedList();
                this.graph = new Graph();
                this.io.query(this.referenceDictionary.getSequence(inputRange.referenceIndex).getSequenceName(),
                        inputRange.startPosition, 
                        inputRange.endPosition);

                // Get first record
                if(null == rec) {
                    rec = this.getNextAlignRecord();
                }

                // skip this range
                if(null != rec &&
                        rec.record.getReferenceIndex() != inputRange.referenceIndex) {
                    continue;
                }

                // TODO: could retrieve just the specified bases
                if(null == this.referenceSequence || inputRange.referenceIndex != this.referenceSequence.getContigIndex()) {
                    this.referenceSequence = this.referenceSequenceFile.getSequence(this.referenceDictionary.getSequence(inputRange.referenceIndex).getSequenceName());
                    if(null == this.referenceSequence) {
                        throw new Exception("Premature EOF in the reference sequence");
                    }
                    else if(this.referenceSequence.getContigIndex() != inputRange.referenceIndex) {
                        throw new Exception("Could not find the reference sequence");
                    }
                }
               
                if(null != rec) {
                    // do an initial prune
                    this.graph.prune(rec.record.getReferenceIndex(), rec.record.getAlignmentStart(), 0);
                }

                while(null != rec) {
                    if(rec.record.getReadUnmappedFlag()) { 
                        // TODO
                        // Print this out somehow in some order somewhere
                        rec = this.getNextAlignRecord();
                        continue;
                    }
                    else if(rec.record.getMappingQuality() < MIN_MAPQ) {
                        // TODO
                        // Print this out somehow in some order somewhere
                        rec = this.getNextAlignRecord();
                        continue;
                    }
                    else {
                        // Make sure that it is sorted
                        if(rec.record.getReferenceIndex() < prevReferenceIndex || (rec.record.getReferenceIndex() == prevReferenceIndex && rec.record.getAlignmentStart() < prevAlignmentStart)) {
                            throw new Exception("SAM/BAM file is not co-ordinate sorted.");
                        }
                        prevReferenceIndex = rec.record.getReferenceIndex();
                        prevAlignmentStart = rec.record.getAlignmentStart();
                    }

                    // process graph 
                    if(this.MAX_QUEUE_SIZE <= this.toAddToGraphList.size()) {
                        this.processToAddToGraphList();
                    }

                    // Add the current record to the graph list
                    this.toAddToGraphList.add(rec);

                    // align
                    if(this.MAX_QUEUE_SIZE <= this.toAlignList.size()) { 
                        ctr = this.processToAlignList(ctr, false);
                    }

                    // get new record
                    rec = this.getNextAlignRecord();
                }
                // process graph 
                this.processToAddToGraphList();
                // align
                ctr = this.processToAlignList(ctr, true);
                // nullify
                this.toAddToGraphList = null;
                this.toAlignList = null;
                this.graph = null;
            }

            // Output any alignments left in the queue 
            while(0 < this.toOutputQueue.size()) {
                rec = this.toOutputQueue.poll();
                this.io.output(rec);
            }
            // Close input/output files
            this.io.closeAll();

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
                else {
                    System.err.println("Total memory usage: " + (Math.round(100 * totalMemory / Math.pow(2, 20)) / 100) + "MB");
                }
                // Run time
                long seconds = (this.endTime - this.startTime) / 1000000000;
                long hours = seconds / 3600; seconds -= hours * 3600; 
                long minutes = seconds / 60; seconds -= minutes* 60; 
                System.err.println("Total execution time: " + hours + "h : " + minutes + "m : " + seconds + "s");
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("Please report bugs to srma-help@lists.sourceforge.net");
            System.exit(1);
        }

        // this is annoying
        QUIET = true;

        return 0;
    }

    private AlignRecord getNextAlignRecord()
        throws Exception
    {
        if(this.io.hasNextAlignRecord()) {
            // HERE: if we have only a subsequence of a contig, we need to check bounds
            return this.io.getNextAlignRecord();
        }
        else {
            return null;
        }
    }

    private void outputProgress(SAMRecord rec, int ctr)
    {
        if(QUIET_STDERR) {
            return;
        }
        else {
            // TODO: enforce column width ?
            int i;
            String outputString = new String("Records processsed: " + ctr + " (last " + rec.getReferenceName() + ":" + rec.getAlignmentStart() + "-" + rec.getAlignmentEnd() + ")");
            int outputStringLength = outputString.length();
            if(this.maxOutputStringLength < outputStringLength) {
                this.maxOutputStringLength = outputStringLength;
            }
            System.err.print("\r" + outputString);
            for(i=outputStringLength;i<this.maxOutputStringLength;i++) {
                System.err.print(" ");
            }
        }
    }

    private void processToAddToGraphList()
        throws Exception
    {
        // Process alignments
        if(0 < this.toAddToGraphList.size()) { 

            int i, size;
            LinkedList<Thread> threads = null;
            LinkedList<LinkedList<AlignRecord>> toAddToGraphThreadLists = null;
            LinkedList<LinkedList<AlignRecord>> toAlignThreadLists = null;


            if(0 == this.toAlignList.size() 
                    && this.graph.contig != this.toAddToGraphList.getFirst().record.getReferenceIndex()+1) 
            {
                // Move to a new contig
                this.graph.prune(this.toAddToGraphList.getFirst().record.getReferenceIndex(),
                        this.toAddToGraphList.getFirst().record.getAlignmentStart(),
                        0);
            }

            // Get the records for the threads 
            toAddToGraphThreadLists = toAddToGraphList.getThreadLists(this.NUM_THREADS, this.graph.contig);

            // Create threads
            threads = new LinkedList<Thread>();
            toAlignThreadLists = new LinkedList<LinkedList<AlignRecord>>();
            for(i=0;i<this.NUM_THREADS;i++) {
                toAlignThreadLists.add(new LinkedList<AlignRecord>());
                threads.add(new GraphThread(i, 
                            toAddToGraphThreadLists.get(i).listIterator(),
                            toAlignThreadLists.get(i)));
            }

            // Start
            for(i=0;i<this.NUM_THREADS;i++) {
                threads.get(i).start();
            }

            // Join
            for(i=0;i<this.NUM_THREADS;i++) {
                threads.get(i).join();
            }

            // Copy records to be re-aligned
            List<ListIterator<AlignRecord>> iters = new LinkedList<ListIterator<AlignRecord>>();
            for(i=size=0;i<this.NUM_THREADS;i++) {
                size += toAlignThreadLists.get(i).size();
                iters.add(toAlignThreadLists.get(i).listIterator());
            }

            for(i=0;0<size;i++) {
                if(this.NUM_THREADS <= i) {
                    i=0;
                }
                if(iters.get(i).hasNext()) {
                    AlignRecord rec = iters.get(i).next();

                    if(null != rec.node 
                            && this.graph.contig == rec.record.getReferenceIndex()+1
                            && this.recordAlignmentStartContained(rec.record))
                    {
                        if(this.withinInsertSizeRange(rec.record)) {
                            // align then output
                            this.toAlignList.add(rec);
                        }
                        else {
                            // just otuput
                            this.toOutputQueue.add(rec);
                        }
                    }
                    size--;
                }
            }
        }
    }

    private int processToAlignList(int ctr, boolean flush)
        throws Exception
    {
        SAMRecord lastSAMRecord = null;

        // Process available alignments
        if(0 < this.toAlignList.size()) { 

            // Check if we are in bounds 
            if(flush || 
                    this.toAlignList.getFirst().record.getAlignmentEnd() + OFFSET < this.toAlignList.getLast().record.getAlignmentStart()) 
            { 
                int i, size;
                LinkedList<Thread> threads = null;
                LinkedList<LinkedList<AlignRecord>> toAlignThreadLists = null;

                // HERE
                /*
                   for(i=this.graph.position_start;i<=this.graph.position_end;i++) {
                   System.err.println(i + "\t"
                   + this.graph.getCoverage(i));
                   }
                   */

                // Get thread data
                if(flush) {
                    toAlignThreadLists = this.toAlignList.getThreadLists(this.NUM_THREADS, this.graph.contig);
                }
                else {
                    toAlignThreadLists = this.toAlignList.getAlignRecordThreadLists(this.NUM_THREADS, 
                            this.graph.contig,
                            this.toAlignList.getLast().record.getAlignmentStart() - this.OFFSET);
                }

                // Create threads
                threads = new LinkedList<Thread>();
                for(i=0;i<this.NUM_THREADS;i++) {
                    threads.add(new AlignThread(i, 
                                this.io.programRecord, 
                                toAlignThreadLists.get(i).listIterator()));
                }

                // Start
                for(i=0;i<this.NUM_THREADS;i++) {
                    threads.get(i).start();
                }

                // Join
                for(i=0;i<this.NUM_THREADS;i++) {
                    threads.get(i).join();
                }

                // Output the alignments
                List<ListIterator<AlignRecord>> iters = new LinkedList<ListIterator<AlignRecord>>();
                for(i=size=0;i<this.NUM_THREADS;i++) {
                    size += toAlignThreadLists.get(i).size();
                    iters.add(toAlignThreadLists.get(i).listIterator());
                }
                for(i=0;0<size;i++) {
                    if(this.NUM_THREADS <= i) {
                        i=0;
                    }
                    if(iters.get(i).hasNext()) {
                        AlignRecord rec = iters.get(i).next();
                        lastSAMRecord = rec.record;
                        this.toOutputQueue.add(rec);
                        ctr++;

                        size--;
                    }
                }

                // Prune the graph
                if(null != lastSAMRecord) {
                    this.outputProgress(lastSAMRecord, ctr);
                    if(0 < toAlignList.size()) {
                        this.graph.prune(toAlignList.getFirst().record.getReferenceIndex(), toAlignList.getFirst().record.getAlignmentStart(), this.OFFSET);
                    }
                    else {
                        this.graph.prune(lastSAMRecord.getReferenceIndex(), lastSAMRecord.getAlignmentStart(), this.OFFSET);
                    }
                }
            }
        }
        if(flush && null != lastSAMRecord) {
            this.outputProgress(lastSAMRecord, ctr);
        }
        lastSAMRecord = null;

        // Output alignments
        while(0 < this.toOutputQueue.size()) {
            AlignRecord rec = this.toOutputQueue.peek();
            // alignment could have moved (+OFFSET), with another moving (-OFFSET) 
            if(rec.record.getReferenceIndex() + 1 != graph.contig || // different contig
                    rec.record.getAlignmentStart() + 2*OFFSET < graph.position_start) { // other alignments will not be less than
                this.io.output(this.toOutputQueue.poll());
                    }
            else { // other alignments could be less than
                break;
            }
        }

        return ctr;
    }

    private boolean recordAlignmentStartContained(SAMRecord rec) 
    {
        int recReferenceIndex = -1; 
        int recAlignmentStart = -1;

        if(null == this.outputRange) { // no more ranges
            return false;
        }

        recReferenceIndex = rec.getReferenceIndex();
        recAlignmentStart = rec.getAlignmentStart();
        while(this.outputRange.referenceIndex < recReferenceIndex 
                || (this.outputRange.referenceIndex == recReferenceIndex 
                    && this.outputRange.endPosition < recAlignmentStart)) 
        { // find a new range
            if(!this.outputRangesIterator.hasNext()) { // no more ranges
                this.outputRange = null;
                return false;
            }
            this.outputRange = this.outputRangesIterator.next();
        }
        if(recReferenceIndex < this.outputRange.referenceIndex 
                || (recReferenceIndex == this.outputRange.referenceIndex
                    && recAlignmentStart < this.outputRange.startPosition)) { // before range
            // not within range
            return false;
                    }
        else {
            // must be within range
            return true;
        }
    }

    private void InsertSizeRangeSet() 
        throws Exception
    {
        String low, high;
        int colon;

        if(null == INSERT_SIZE_RANGE) {
            InsertSizeRangeLow = Integer.MIN_VALUE;
            InsertSizeRangeHigh = Integer.MAX_VALUE;
            return;
        }

        colon = INSERT_SIZE_RANGE.indexOf(":");
        if(colon < 0 || INSERT_SIZE_RANGE.length() <= colon) {
            throw new Exception("INSERT_SIZE_RANGE was improperly specified");
        }

        InsertSizeRangeLow = Integer.parseInt(INSERT_SIZE_RANGE.substring(0, colon));
        InsertSizeRangeHigh = Integer.parseInt(INSERT_SIZE_RANGE.substring(colon+1));

        if(InsertSizeRangeHigh < InsertSizeRangeLow) {
            throw new Exception("INSERT_SIZE_RANGE was improperly specified");
        }
    }

    private boolean withinInsertSizeRange(SAMRecord rec)
    {
        int isize;
        if(null == this.INSERT_SIZE_RANGE) {
            return true;
        }
        if(!rec.getReadPairedFlag()
                || rec.getMateUnmappedFlag()) {
            // not paired
            // or only one end maps
            return true;
                }
        // assumes: paired and both mapped
        isize = rec.getInferredInsertSize();
        if(rec.getFirstOfPairFlag()) { // first record, 5'
            if(InsertSizeRangeLow <= isize && isize <= InsertSizeRangeHigh) {
                return true;
            }
            else {
                return false;
            }
        } 
        else { // second record, 3'
            if(InsertSizeRangeLow <= -isize && -isize <= InsertSizeRangeHigh) {
                return true;
            }
            else {
                return false;
            }
        }
    }

    private class AlignThread extends Thread {

        private int threadID;
        private SAMProgramRecord programRecord;
        private ListIterator<AlignRecord> iter;

        public AlignThread(int threadID,
                SAMProgramRecord programRecord,
                ListIterator<AlignRecord> iter)
        {
            this.threadID = threadID;
            this.programRecord = programRecord;
            this.iter = iter;
        }

        public void run() 
        {
            try {

                // Align each record
                while(iter.hasNext()) {
                    AlignRecord curAlignRecord = iter.next();
                    SAMRecord curSAMRecord = curAlignRecord.record;
                    Node curNode = curAlignRecord.node;

                    // Align - this will overwrite/change the alignment
                    Align.align(graph,
                            curSAMRecord,
                            curNode,
                            referenceSequence,
                            this.programRecord,
                            OFFSET,
                            alleleCoverageCutoffs,
                            CORRECT_BASES,
                            USE_SEQUENCE_QUALITIES,
                            MAXIMUM_TOTAL_COVERAGE,
                            MAX_HEAP_SIZE);
                }

            } catch (Exception e) {
                e.printStackTrace();
                System.err.println("Please report bugs to srma-help@lists.sourceforge.net");
                System.exit(1);
            }
        }
    }

    private class GraphThread extends Thread {
        private int threadID;
        private ListIterator<AlignRecord> iterAlignRecords;
        List<AlignRecord> toAlignThreadList;

        public GraphThread(int threadID,
                ListIterator<AlignRecord> iterAlignRecords,
                List<AlignRecord> toAlignThreadList)
        {
            this.threadID = threadID;
            this.iterAlignRecords = iterAlignRecords;
            this.toAlignThreadList = toAlignThreadList;
        }

        public void run()
        {
            while(this.iterAlignRecords.hasNext()) {
                // Get record
                AlignRecord rec = this.iterAlignRecords.next();

                Node recNode = null;

                if(graph.contig != rec.record.getReferenceIndex()+1) {
                    break;
                }

                // Add to the graph 
                try {
                    recNode = graph.addSAMRecord(rec.record, referenceSequence);
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Please report bugs to srma-help@lists.sourceforge.net");
                    System.exit(1);
                }

                // Keep track of start node
                rec.setNode(recNode);
                toAlignThreadList.add(rec);
            }
        }
    }
}
