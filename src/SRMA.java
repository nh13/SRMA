/*
 * LICENSE to be determined
 */
package srma;

import srma.Align;
import srma.ThreadPoolLinkedList;

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

    public final String PROGRAM_VERSION="0.1.6";
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
    private String maxOutputString = null;

    ReferenceSequenceFile referenceSequenceFile = null; 
    private ReferenceSequence referenceSequence = null;
    private SAMSequenceDictionary referenceDictionary = null;

    private ThreadPoolLinkedList<SAMRecord> toAddToGraphSAMRecordList = null;
    private ThreadPoolLinkedList<SAMRecord> toProcessSAMRecordList = null;
    private ThreadPoolLinkedList<Node> toProcessSAMRecordNodeList = null;
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
        this.maxOutputString = new String("");
        this.alleleCoverageCutoffs = new AlleleCoverageCutoffs(MINIMUM_ALLELE_COVERAGE, MINIMUM_ALLELE_PROBABILITY, QUIET_STDERR);

        try { 
            this.startTime = System.nanoTime();

            if(1 < this.NUM_THREADS) {
                System.err.println("** Warning: option NUM_THREADS currently does not increase performance significantly **");
            }

            // Check input files
            IoUtil.assertFileIsReadable(INPUT);
            IoUtil.assertFileIsReadable(REFERENCE);

            // Initialize basic input/output files
            this.toAddToGraphSAMRecordList = new ThreadPoolLinkedList<SAMRecord>(NUM_THREADS);
            this.toProcessSAMRecordList = new ThreadPoolLinkedList<SAMRecord>(NUM_THREADS);
            this.toProcessSAMRecordNodeList = new ThreadPoolLinkedList<Node>(NUM_THREADS);
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
            this.referenceSequenceFile = new IndexedFastaSequenceFile(REFERENCE);
            if(!this.referenceSequenceFile.isIndexed()) {
                throw new Exception("Reference sequence file was not indexed.");
            }
            this.referenceDictionary = this.referenceSequenceFile.getSequenceDictionary();

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

                this.referenceSequence = this.referenceSequenceFile.getSequence(this.referenceDictionary.getSequence(this.inputRange.referenceIndex).getSequenceName());
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
            // Continue while either
            // - there is input
            // - there are records to add to the graph
            // - there are records to re-align
            while(null != rec ||
                    0 < this.toAddToGraphSAMRecordList.size() ||
                    0 < this.toProcessSAMRecordList.size()) 
            {
                /*
                   System.err.println("WHILE LOOP (" + 
                   (null != rec) +
                   ") (" +
                   (0 < this.toAddToGraphSAMRecordList.size()) +
                   ") (" +
                   (0 < this.toProcessSAMRecordList.size()) +
                   ")");
                   */
                if(null != rec) {
                    if(rec.getReadUnmappedFlag()) { 
                        // TODO
                        // Print this out somehow in some order somewhere
                        rec = this.getNextSAMRecord();
                        continue;
                    }
                    else if(rec.getMappingQuality() < MIN_MAPQ) {
                        // TODO
                        // Print this out somehow in some order somewhere
                        rec = this.getNextSAMRecord();
                        continue;
                    }
                    else {
                        // Make sure that it is sorted
                        if(rec.getReferenceIndex() < prevReferenceIndex || (rec.getReferenceIndex() == prevReferenceIndex && rec.getAlignmentStart() < prevAlignmentStart)) {
                            throw new Exception("SAM/BAM file is not co-ordinate sorted.");
                        }
                        prevReferenceIndex = rec.getReferenceIndex();
                        prevAlignmentStart = rec.getAlignmentStart();
                    }
                }

                // Process the graph if either:
                // - no more input
                // - we are moving to a new contig
                // - we have reached the queue size
                // RANGES? TODO
                if(null == rec ||
                        this.graph.contig != rec.getReferenceIndex()+1 
                        || this.MAX_QUEUE_SIZE <= this.toAddToGraphSAMRecordList.size())

                {
                    // add all to the graph, up to the next contig
                    // Threaded add to the graph
                    this.processGraph();
                }

                if(null != rec) {
                    // If we are moving to a new contig, force processing
                    if(this.graph.contig != rec.getReferenceIndex()+1) {
                        // Process the rest of the reads
                        ctr = this.processList(programRecord, ctr, true, this.NUM_THREADS);

                        // Get new reference sequence
                        this.referenceSequence = this.referenceSequenceFile.getSequence(this.referenceDictionary.getSequence(rec.getReferenceIndex()).getSequenceName());

                        if(null == this.referenceSequence) {
                            throw new Exception("Premature EOF in the reference sequence");
                        }
                        else if(this.referenceSequence.getContigIndex() != rec.getReferenceIndex()) {
                            throw new Exception("Could not find the reference sequence");
                        }
                    }

                    // Add the current record to the graph addition list
                    this.toAddToGraphSAMRecordList.add(rec);
                }

                // Process the available reads
                if(null == rec &&
                        0 == this.toAddToGraphSAMRecordList.size()) 
                {
                    // flush
                    ctr = this.processList(programRecord, ctr, true, this.NUM_THREADS);
                }
                else {
                    // there may be more to come
                    ctr = this.processList(programRecord, ctr, false, this.NUM_THREADS);
                }

                // get new record
                if(null != rec) {
                    rec = this.getNextSAMRecord();
                }
            }

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


            this.referenceSequence = this.referenceSequenceFile.getSequence(this.referenceDictionary.getSequence(this.inputRange.referenceIndex).getSequenceName());
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
        else {
            // TODO: enforce column width ?
            String outputString = new String("Records processsed: " + ctr + " (last " + rec.getReferenceName() + ":" + rec.getAlignmentStart() + "-" + rec.getAlignmentEnd() + ")");
            int outputStringLength = outputString.length();
            if(this.maxOutputStringLength < outputStringLength) {
                int i;
                for(i=this.maxOutputStringLength;i< outputStringLength;i++) { // pad with blanks
                    this.maxOutputString += " ";
                }
                this.maxOutputStringLength = outputStringLength;
            }
            System.err.print("\r" + outputString + this.maxOutputString);
        }
    }

    private void processGraph()
        throws Exception
    {
        // Process alignments
        if(0 < this.toAddToGraphSAMRecordList.size()) { 

            int i;
            LinkedList<Thread> threads = null;
            LinkedList<LinkedList<Node>> toProcessSAMRecordNodeListThread = null;

            if(0 == this.toProcessSAMRecordList.size() 
                    && this.graph.contig != this.toAddToGraphSAMRecordList.getFirst().getReferenceIndex()+1) {
                // Move to a new contig
                this.graph.prune(this.toAddToGraphSAMRecordList.getFirst().getReferenceIndex(),
                        this.toAddToGraphSAMRecordList.getFirst().getAlignmentStart(),
                        0);
                    }

            // Create threads
            threads = new LinkedList<Thread>();
            toProcessSAMRecordNodeListThread = new LinkedList<LinkedList<Node>>();
            for(i=0;i<this.NUM_THREADS;i++) {
                toProcessSAMRecordNodeListThread.add(new LinkedList<Node>());
                threads.add(new GraphThread(i, 
                            this.graph,
                            this.referenceSequence,
                            this.toAddToGraphSAMRecordList.listIterator(i),
                            toProcessSAMRecordNodeListThread.get(i)));
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
            int first = toAddToGraphSAMRecordList.getFirstIndex();
            while(0 < toAddToGraphSAMRecordList.size()) {
                SAMRecord rec = toAddToGraphSAMRecordList.getFirst();
                if(this.graph.contig != rec.getReferenceIndex()+1) {
                    break;
                }
                rec = toAddToGraphSAMRecordList.removeFirst();

                Node recNode = toProcessSAMRecordNodeListThread.get(first).removeFirst();
                if(null != recNode &&
                        (!this.useRanges || this.recordAlignmentStartContained(rec))) 
                {
                    this.toProcessSAMRecordList.add(rec);
                    this.toProcessSAMRecordNodeList.add(recNode);
                }
                first++;
                if(this.NUM_THREADS <= first) {
                    first = 0;
                }
            }

            // Rebalance
            this.toAddToGraphSAMRecordList.rebalance();
        }
    }

    private int processList(SAMProgramRecord programRecord, int ctr, boolean flush, int NUM_THREADS)
        throws Exception
    {
        SAMRecord lastSAMRecord = null;

        // Check if we should process
        if(0 == this.toProcessSAMRecordList.size() ||
                (!flush && this.toProcessSAMRecordList.size() < this.MAX_QUEUE_SIZE)) 
        {
            return ctr;
        }

        // Process alignments
        if(0 < this.toProcessSAMRecordList.size()) { 


            // Check if we are in bounds 
            if(flush || 
                    this.toProcessSAMRecordList.getFirst().getAlignmentEnd() + OFFSET < this.toProcessSAMRecordList.getLast().getAlignmentStart()) 
            { 
                int i;
                LinkedList<Thread> threads = null;
                List<Integer> numProcessedSAMRecordList = null;
                int numProcessed = 0; 

                // Create threads
                threads = new LinkedList<Thread>();
                numProcessedSAMRecordList = new LinkedList<Integer>();

                for(i=0;i<NUM_THREADS;i++) {
                    numProcessedSAMRecordList.add(new Integer(0));
                    threads.add(new AlignThread(i, 
                                programRecord, 
                                flush,
                                this.toProcessSAMRecordList.getLast().getAlignmentStart(),
                                this.toProcessSAMRecordList.listIterator(i),
                                this.toProcessSAMRecordNodeList.listIterator(i),
                                numProcessedSAMRecordList));
                }

                // Start
                for(i=0;i<NUM_THREADS;i++) {
                    threads.get(i).start();
                }

                // Join
                for(i=numProcessed=0;i<NUM_THREADS;i++) {
                    threads.get(i).join();
                }

                // Get the number processed
                for(i=numProcessed=0;i<NUM_THREADS;i++) {
                    numProcessed += numProcessedSAMRecordList.get(i); 
                }


                // Output the alignments
                for(i=0;i<numProcessed;i++) {
                    lastSAMRecord = this.toProcessSAMRecordList.removeFirst();
                    this.toOutputSAMRecordPriorityQueue.add(lastSAMRecord);
                    this.toProcessSAMRecordNodeList.removeFirst();
                    ctr++;
                }

                // Prune the graph
                if(null != lastSAMRecord) {
                    this.outputProgress(lastSAMRecord, ctr);
                    if(0 < toProcessSAMRecordList.size()) {
                        this.graph.prune(toProcessSAMRecordList.getFirst().getReferenceIndex(), toProcessSAMRecordList.getFirst().getAlignmentStart(), this.OFFSET);
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
        while(0 < this.toOutputSAMRecordPriorityQueue.size()) {
            SAMRecord curSAMRecord = this.toOutputSAMRecordPriorityQueue.peek();
            // alignment could have moved (+OFFSET), with another moving (-OFFSET) 
            if(flush || curSAMRecord.getAlignmentStart() + 2*OFFSET < graph.position_start) { // other alignments will not be less than
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

    private class AlignThread extends Thread {

        private int threadID;
        private SAMProgramRecord programRecord;
        private boolean flush;
        private int lastAlignmentStart;
        private ListIterator<SAMRecord> iterSAMRecords;
        private ListIterator<Node> iterNodes;
        private List<Integer> numProcessedSAMRecordList;

        public AlignThread(int threadID,
                SAMProgramRecord programRecord,
                boolean flush,
                int lastAlignmentStart,
                ListIterator<SAMRecord> iterSAMRecords,
                ListIterator<Node> iterNodes,
                List<Integer> numProcessedSAMRecordList)
        {
            this.threadID = threadID;
            this.programRecord = programRecord;
            this.flush = flush;
            this.lastAlignmentStart = lastAlignmentStart;
            this.iterSAMRecords = iterSAMRecords;
            this.iterNodes = iterNodes;
            this.numProcessedSAMRecordList = numProcessedSAMRecordList;
        }

        public void run() 
        {
            try {
                int prevReferenceIndex=-1;
                int prevAlignmentStart=-1;
                int numProcessed = 0;
                // Do stuff
                while(iterSAMRecords.hasNext()) {
                    SAMRecord curSAMRecord = iterSAMRecords.next();

                    if(!flush && this.lastAlignmentStart <= curSAMRecord.getAlignmentEnd() + OFFSET) { 
                        break;
                    }

                    if(curSAMRecord.getReferenceIndex() < prevReferenceIndex || (curSAMRecord.getReferenceIndex() == prevReferenceIndex && curSAMRecord.getAlignmentStart() < prevAlignmentStart)) {

                        System.err.println(curSAMRecord.getReferenceIndex() + " < " + prevReferenceIndex + " " + curSAMRecord.getAlignmentStart() + " < " + prevAlignmentStart);
                        throw new Exception("SRMA bug: thread SAMRecords are not co-ordinate sorted.");
                    }

                    prevReferenceIndex = curSAMRecord.getReferenceIndex();
                    prevAlignmentStart = curSAMRecord.getAlignmentStart();

                    // Align - this will overwrite/change the alignment
                    Align.align(graph,
                            curSAMRecord,
                            iterNodes.next(),
                            referenceSequence,
                            this.programRecord,
                            OFFSET,
                            alleleCoverageCutoffs,
                            CORRECT_BASES,
                            USE_SEQUENCE_QUALITIES,
                            MAXIMUM_TOTAL_COVERAGE,
                            MAX_HEAP_SIZE);
                    numProcessed++;
                }

                // Store the number of records numProcessed by this thread
                synchronized(numProcessedSAMRecordList) {
                    numProcessedSAMRecordList.set(threadID, numProcessed);
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
        private Graph graph;
        private ReferenceSequence referenceSequence;
        private ListIterator<SAMRecord> iterSAMRecords;
        List<Node> toProcessSAMRecordNodeListThread;

        public GraphThread(int threadID,
                Graph graph,
                ReferenceSequence referenceSequence,
                ListIterator<SAMRecord> iterSAMRecords,
                List<Node> toProcessSAMRecordNodeListThread)
        {
            this.threadID = threadID;
            this.graph = graph;
            this.referenceSequence = referenceSequence;
            this.iterSAMRecords = iterSAMRecords;
            this.toProcessSAMRecordNodeListThread = toProcessSAMRecordNodeListThread;
        }

        public void run()
        {
            while(this.iterSAMRecords.hasNext()) {
                // Get record
                SAMRecord rec = this.iterSAMRecords.next();

                if(this.graph.contig != rec.getReferenceIndex()+1) {
                    break;
                }

                // Add to the graph 
                Node recNode = null;
                try {
                    recNode = this.graph.addSAMRecord(rec, this.referenceSequence);
                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Please report bugs to srma-help@lists.sourceforge.net");
                    System.exit(1);
                }

                // Keep track of start node
                toProcessSAMRecordNodeListThread.add(recNode);
            }
        }
    }
}
