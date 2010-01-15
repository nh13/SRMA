/*
 * LICENSE to be determined
 */
package srma;

import srma.Align;

import net.sf.samtools.*;
import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.*;

import java.io.File;
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
     * - single contig
     * - can fit entire partial order graph in memory
     * - SAM entries are on the forward strand
     * */
    protected int doWork() {
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
        SAMFileReader in = new SAMFileReader(INPUT);
        final SAMFileHeader header = in.getFileHeader();
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, System.out);
        Graph graph = new Graph(header, referenceSequences);

        // Go through each SAM record
        int ctr = 0;
        try {
            for (final SAMRecord rec : in) {
                // TODO: Make sure that it is sorted
                // Add only if it is from the same contig
                if(graph.contig != rec.getReferenceIndex()) {
                    // Process the rest of the reads
                    while(0 < list.size()) {
                        out.addAlignment(Align.Align(graph, list.removeFirst(), OFFSET, COVERAGE));
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
                
                // HERE
                System.err.println("Printing GRAPH:");
                graph.print(System.err);

                // TODO: check if we should process ... 
                while(0 < list.size() && graph.position_start + OFFSET <= list.getFirst().getAlignmentStart()) {
                    ctr++;
                    System.err.println("ctr="+ctr);
                    SAMRecord curSAMRecord = list.removeFirst();
                    graph.prune(curSAMRecord.getAlignmentStart() - OFFSET);
                    out.addAlignment(Align.Align(graph, curSAMRecord, OFFSET, COVERAGE));
                }
            }
            // Process the rest of the reads
            System.err.println("Finishing!");
            while(0 < list.size()) {
                ctr++;
                System.err.println("ctr="+ctr);
                out.addAlignment(Align.Align(graph, list.removeFirst(), OFFSET, COVERAGE));
            }
        } catch (Exception e) {
            System.err.println(e.toString());
            e.printStackTrace();
            System.exit(1);
        }
        // Close files
        in.close();
        out.close();

        return 0;
    }
}
