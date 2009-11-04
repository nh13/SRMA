/*
 * LICENSE to be determined
 */
package srma;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;

import java.util.*;
import srma.Node;

public class Graph {
    int contig;
    int position_start;
    int position_end;
    List<Node> referenceNodes;
    SAMFileHeader header;
    List<ReferenceSequence> sequences;

    public Graph(SAMFileHeader header, List<ReferenceSequence> sequences)
    {
        this.header = header;
        this.sequences = sequences;
        this.contig = 0;
        this.position_start = -1;
        this.position_end = -1;
        referenceNodes = new ArrayList<Node>(); 
    }

    public void addSAMRecord(SAMRecord record) throws Exception
    {
        Alignment alignment;

        // Get the alignment
        alignment = new Alignment(record);

        /*
        System.out.println("refr:" + new String(alignment.reference));
        System.out.println("read:" + new String(alignment.read));
        */

        int i, j, ref_i, seq_i;
        Node prev=null, cur=null;

        for(i=ref_i=seq_i=0;i<alignment.length;i++) { // go through the alignment
            // TODO: lower/upper case?
            if(alignment.read[i] == alignment.reference[i]) { // match
                cur = insertMatch(alignment.read[i],
                        record.getReferenceIndex(),
                        record.getAlignmentStart() + ref_i,
                        prev);
                ref_i++; seq_i++;
            }
            else if(alignment.reference[i] == Alignment.GAP) { // insertion
                assert Alignment.GAP != alignment.read[i];
                for(cur=null,j=0;null != prev && j < prev.next.size();j++) {
                    Node tmpNode = prev.next.get(j);
                    if(Node.NodeType.MATCH != tmpNode.type) {
                        cur = tmpNode;
                    }
                }
                if(null == cur) { // No such insertion
                    cur = new Node(alignment.read[i],
                            Node.NodeType.INSERTION,
                            prev.contig,
                            prev.position,
                            prev);
                }
                seq_i++;
            }
            else if(alignment.read[i] == Alignment.GAP) { // deletion
                assert Alignment.GAP != alignment.read[i];
                ref_i++;
                cur = prev;
                // Skip over deletion.  Note: impossible to end with a deletion.
                while(Alignment.GAP == alignment.read[i]) {
                    assert i < alignment.length;
                    assert Alignment.GAP != alignment.reference[i];
                    i++;
                    ref_i++;
                }
                i--;
            }
            else { // mismatch
                cur = new Node(alignment.read[i],
                            Node.NodeType.MISMATCH,
                            record.getReferenceIndex(),
                            record.getAlignmentStart() + ref_i,
                            prev);
            }
            assert null != cur;
            if(null != prev 
                    && prev != cur) {
                if(!prev.next.contains(cur)) {
                    prev.next.add(cur);
                }
                if(!cur.prev.contains(cur)) {
                    cur.prev.add(prev);
                }
                fixReverseInsertion(cur);
            }
            prev = cur;
        }
    }

    private Node insertMatch(byte base, int contig, int position, Node node)
    {
        // TODO
        return null;
    }

    private void fixReverseInsertion(Node node)
    {
        // TODO
    }

    // Inner class
    public class Alignment {
        public static final int GAP = '-';
        int length;
        byte read[];
        byte reference[];

        public Alignment(SAMRecord record) throws Exception
        {
            Cigar cigar;
            List<CigarElement> cigarElements;
            Iterator<CigarElement> iter;
            CigarElement cigarElement;
            int cigarElementLength;
            CigarOperator cigarElementOperator;
            int readIndex, referenceIndex;
            int i, index;
            byte referenceBases[];
            byte readBases[];
            int alignmentStart;

            alignmentStart = record.getAlignmentStart();
            cigar = record.getCigar();
            readBases = record.getReadBases();

            // Get alignment length
            this.length = 0;
            cigarElements = cigar.getCigarElements();
            iter = cigarElements.iterator();
            while(iter.hasNext()) {
                this.length += iter.next().getLength();
            }

            this.read = new byte[this.length];
            this.reference = new byte[this.length];

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
                            this.reference[index] = referenceBases[alignmentStart - 1 + referenceIndex];
                            this.read[index] = readBases[readIndex]; 
                            referenceIndex++;
                            readIndex++;
                            break;
                        case D:
                            this.reference[index] = referenceBases[alignmentStart - 1 + referenceIndex];
                            this.read[index] = Alignment.GAP;
                            referenceIndex++;
                            break;
                        case I:
                            this.reference[index] = Alignment.GAP; 
                            this.read[index] = readBases[readIndex]; 
                            readIndex++;
                            break;
                        default:
                            throw new Exception("Illegal Cigar Operator: " + cigarElementOperator);
                    }
                    index++;
                }
            }
        }
    }
}
