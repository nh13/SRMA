/*
 * LICENSE to be determined
 */
package srma;

import srma.Align;

import net.sf.samtools.*;
import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.*;

import java.io.*;
import java.util.*;

public class SRMA extends CommandLineProgram { 

    public static final String OFFSET_SHORT_NAME = "O"; 
    public static final String COVERAGE_SHORT_NAME = "C";

    @Usage public final String USAGE = getStandardUsagePreamble() + "Prints a SAM or BAM file to the screen.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file to view.")
        public File INPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference FASTA file.")
        public File REFERENCE;
    @Option(shortName=OFFSET_SHORT_NAME, doc="The alignment offset.")
        public int OFFSET=0;
    @Option(shortName=COVERAGE_SHORT_NAME, doc="The minimum haploid coverage for the consensus")
        public int COVERAGE=1;

    public static void main(final String[] args) {
        new SRMA().instanceMain(args);
    }

    /*
     * Current assumptions:
     * - can fit entire partial order graph in memory
     * */
    protected int doWork() 
    {
        List<ReferenceSequence> referenceSequences = new ArrayList();

        IoUtil.assertFileIsReadable(INPUT);

        // Get references
        IoUtil.assertFileIsReadable(REFERENCE);
        ReferenceSequenceFileFactory referenceSequenceFileFactory = new ReferenceSequenceFileFactory();
        ReferenceSequenceFile referenceSequenceFile = referenceSequenceFileFactory.getReferenceSequenceFile(REFERENCE);
        ReferenceSequence referenceSequence = null;
        do {
            referenceSequence = referenceSequenceFile.nextSequence();
            if(null != referenceSequence) {
                referenceSequences.add(referenceSequence);
            }
        } while(null != referenceSequence);

        // Initialize graph, input/output files
        LinkedList<SAMRecord> list = new LinkedList<SAMRecord>();
        SAMFileReader in = new SAMFileReader(INPUT, true);
        final SAMFileHeader header = in.getFileHeader();
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, System.out);
        Graph graph = new Graph(header, referenceSequences);

        // Go through each SAM record
        System.err.println("");
        try {
            //PrintStream graphOut = new PrintStream("graph.txt");
            for (final SAMRecord rec : in) {
                // TODO: Make sure that it is sorted
                // Add only if it is from the same contig
                if(graph.contig != rec.getReferenceIndex()+1) {
                    // Process the rest of the reads
                    while(0 < list.size()) {
                        SAMRecord curSAMRecord = list.removeFirst();
                        System.err.print("\rAL:" + curSAMRecord.getAlignmentStart() + ":" + curSAMRecord.getAlignmentEnd() + ":" + curSAMRecord.toString());
                        out.addAlignment(Align.Align(graph, curSAMRecord, OFFSET, COVERAGE));
                    }
                }

                // Add to the graph 
                try {
                    graph.addSAMRecord(rec);
                } catch (Graph.GraphException e) {
                    if(Graph.GraphException.NOT_IMPLEMENTED != e.type) {
                        throw e;
                    }
                }
                list.add(rec);

                //System.err.println("Printing GRAPH:");
                //graph.print(System.err);
                //graphOut.println("GRAPH"); // HERE
                //graph.print(graphOut); // HERE

                // HERE
                System.err.print("\rIN:" + rec.getAlignmentStart() + ":" + rec.getAlignmentEnd() + ":" + rec.toString());

                // TODO: check if we should process ... 
                while(0 < list.size() && list.getFirst().getAlignmentEnd() + OFFSET < list.getLast().getAlignmentStart()) {
                    SAMRecord curSAMRecord = list.removeFirst();
                    System.err.print("\rAL:" + curSAMRecord.getAlignmentStart() + ":" + curSAMRecord.getAlignmentEnd() + ":" + curSAMRecord.toString());
                    graph.prune(curSAMRecord.getReferenceIndex(), curSAMRecord.getAlignmentStart(), OFFSET); 
                    out.addAlignment(Align.Align(graph, curSAMRecord, OFFSET, COVERAGE));
                }
            }
            // Process the rest of the reads
            while(0 < list.size()) {
                SAMRecord curSAMRecord = list.removeFirst();
                System.err.print("\rAL:" + curSAMRecord.getAlignmentStart() + ":" + curSAMRecord.getAlignmentEnd() + ":" + curSAMRecord.toString());
                //graphOut.println("HERE2\t" + list.getFirst().getAlignmentEnd() + ":" + list.getLast().getAlignmentStart());
                out.addAlignment(Align.Align(graph, curSAMRecord, OFFSET, COVERAGE));
            }
            //graphOut.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        // Close files
        in.close();
        out.close();

        // HERE

        return 0;
    }
}
