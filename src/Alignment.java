/*
 * LICENSE to be determined
 */
package srma;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;

import java.util.*;

public class Alignment {
    public static final int GAP = '-';
    int length;
    byte read[];
    byte reference[];

    public Alignment(SAMRecord record, List<ReferenceSequence> sequences) throws Exception
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
                // TODO: make sure these are upper case
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
