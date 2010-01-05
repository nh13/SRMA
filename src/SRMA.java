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
        SAMFileReader in = new SAMFileReader(INPUT);
        final SAMFileHeader header = in.getFileHeader();

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

        Graph graph = new Graph(header, referenceSequences);

        for (final SAMRecord rec : in) {
            // TODO: Make sure that it is sorted
            try {
                // Add only if it is from the same contig
                if(graph.contig == rec.getReferenceIndex()) {
                    graph.addSAMRecord(rec);
                }

                // Process
                // TODO:

                // Add...
                if(graph.contig != rec.getReferenceIndex()) {
                    throw new Exception("Error.  Multiple contigs not supported.");
                    //graph.addSAMRecord(rec);
                }
            } catch (Exception e) {
                System.err.println(e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }
        in.close();

        // DEBUGGING
        //graph.print();

        /* Align sam records */
        in = new SAMFileReader(INPUT);
        final SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, System.out);
        for (final SAMRecord rec : in) {
            // TODO: Make sure that it is sorted
            try {
                // Align the data
                //out.addAlignment(rec); // HERE
                Align.Align(graph, rec, OFFSET, COVERAGE);
                out.addAlignment(rec);
            } catch (Exception e) {
                System.err.println(e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }
        in.close();
        out.close();
        
        /*
           final SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, System.out);
           for (final SAMRecord rec : in) {
           if (System.out.checkError()) {
           return 0;
           }

           out.addAlignment(rec);
           }
           out.close();
           */

        return 0;

    }
}
