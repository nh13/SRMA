/*
 * LICENSE to be determined
 */
package srma;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;

import java.util.*;
//import srma.Node;

public class Graph {
    int contig;
    int position_start;
    int position_end;
    List<Node> referenceNodes;
    SAMFileHeader header;
    List<ReferenceSequence> sequences;

    // init
    // destroy
    // add bam
    // insert match
    // fix reverse insertion
    // bases length

    public Graph(SAMFileHeader header, List<ReferenceSequence> sequences)
    {
        this.header = header;
        this.sequences = sequences;
        this.contig = 0;
        this.position_start = -1;
        this.position_end = -1;
        referenceNodes = new ArrayList<Node>(); 
    }

    public void update(SAMRecord record) throws Exception
    {
        Alignment alignment;

        // Get the alignment
        alignment = this.getAlignment(record);

        //System.out.println("refr:" + new String(alignment.reference));
        //System.out.println("read:" + new String(alignment.read));
        
        // TODO
    }

    private void insertMatchNode(char base, int contig, int position, Node node)
    {
        // TODO
    }

    private void fixReverseInsertion(Node node)
    {
        // TODO
    }


    private Alignment getAlignment(SAMRecord record) throws Exception
    {
        Cigar cigar;
        List<CigarElement> cigarElements;
        Iterator<CigarElement> iter;
        CigarElement cigarElement;
        byte read[];
        byte reference[];
        int cigarElementLength;
        CigarOperator cigarElementOperator;
        int readIndex, referenceIndex;
        int i, index;
        byte referenceBases[];
        byte readBases[];
        int alignmentStart, alignmentLength;
        
        alignmentStart = record.getAlignmentStart();
        cigar = record.getCigar();
        readBases = record.getReadBases();

        // Get alignment length
        alignmentLength = 0;
        cigarElements = cigar.getCigarElements();
        iter = cigarElements.iterator();
        while(iter.hasNext()) {
            alignmentLength += iter.next().getLength();
        }

        read = new byte[alignmentLength];
        reference = new byte[alignmentLength];

        // Get reference bases
        referenceIndex = record.getReferenceIndex();
        if(referenceIndex < 0) {
            throw new Exception("Reference index out of range: " + referenceIndex);
        }
        referenceBases = sequences.get(referenceIndex).getBases();

        // Copy over alignment
        iter = cigarElements.iterator();
        index = readIndex = referenceIndex = 0;
        while(iter.hasNext()) {
            cigarElement = iter.next();
            cigarElementLength = cigarElement.getLength();
            cigarElementOperator = cigarElement.getOperator();
            for(i=0;i<cigarElementLength;i++) {
                switch(cigarElementOperator) {
                    case M:
                        reference[index] = referenceBases[alignmentStart - 1 + referenceIndex];
                        read[index] = readBases[readIndex]; 
                        referenceIndex++;
                        readIndex++;
                        break;
                    case D:
                        reference[index] = referenceBases[alignmentStart - 1 + referenceIndex];
                        read[index] = Alignment.GAP;
                        referenceIndex++;
                        break;
                    case I:
                        reference[index] = Alignment.GAP; 
                        read[index] = readBases[readIndex]; 
                        readIndex++;
                        break;
                    default:
                        throw new Exception("Illegal Cigar Operator: " + cigarElementOperator);
                }
                index++;
            }
        }

        return new Alignment(read, reference);
    }

    // Inner class
    public class Alignment {
        public static final int GAP = '-';
        byte read[];
        byte reference[];

        public Alignment(byte read[], byte reference[]) 
        {
            this.read = read;
            this.reference = reference;
        }
    }
}
