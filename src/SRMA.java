/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import net.sf.picard.cmdline.*;
import net.sf.picard.io.IoUtil;
import net.sf.picard.reference.*;

import java.io.File;
import java.util.*;

public class SRMA extends CommandLineProgram { 

    @Usage public final String USAGE = getStandardUsagePreamble() + "Prints a SAM or BAM file to the screen.";
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The SAM or BAM file to view.")
        public File INPUT;
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference FASTA file.")
        public File REFERENCE;

    public static void main(final String[] args) {
        new SRMA().instanceMain(args);
    }

    protected int doWork() {
        List<ReferenceSequence> referenceSequences = new ArrayList();
                
        IoUtil.assertFileIsReadable(INPUT);
        final SAMFileReader in = new SAMFileReader(INPUT);
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
        LinkedList<SAMRecord> recordQueue = new LinkedList<SAMRecord>();

        for (final SAMRecord rec : in) {
            // TODO: Make sure that it is sorted
            try {
                // Add it to the queue
                recordQueue.add(rec);

                // Add only if it is from the same contig
                if(graph.contig == rec.getReferenceIndex()) {
                    graph.addSAMRecord(rec);
                }

                // Process
                // TODO:

                // Add...
                if(graph.contig != rec.getReferenceIndex()) {
                    graph.addSAMRecord(rec);
                }
            } catch (Exception e) {
                System.err.println(e.toString());
                e.printStackTrace();
                System.exit(1);
            }
        }

        // DEBUGGING
        graph.print();

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
