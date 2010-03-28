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
    int positions[]; // one for each read base
    int positionsIndex[]; // one for each read base, index into read/reference

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
        int positionsLength;

        alignmentStart = record.getAlignmentStart();
        cigar = record.getCigar();
        readBases = record.getReadBases();

        // Get alignment length
        this.length = positionsLength = 0;
        cigarElements = cigar.getCigarElements();
        iter = cigarElements.iterator();
        while(iter.hasNext()) {
            cigarElement = iter.next();
            this.length += cigarElement.getLength();
            switch(cigarElement.getOperator()) {
                case M:
                case I:
                    positionsLength += cigarElement.getLength();
                    break;
                default:
                    break;
            }
        }

        this.read = new byte[this.length];
        this.reference = new byte[this.length];
        this.positions = new int[positionsLength];
        this.positionsIndex = new int[positionsLength];

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
                        this.positions[readIndex] = alignmentStart + referenceIndex;
                        this.positionsIndex[readIndex] = index;
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
                        this.positions[readIndex] = alignmentStart + referenceIndex - 1;
                        this.positionsIndex[readIndex] = index;
                        readIndex++;
                        break;
                    default:
                        throw new Exception("Illegal Cigar Operator: " + cigarElementOperator);
                }
                index++;
            }
        }

        // left justify
        this.leftJustify();

    }

    private void leftJustify()
    {
        // Left-justify alignment
        int i;
        int prevDel, prevIns, startDel, endDel, startIns, endIns;

        i = prevDel = prevIns = 0;
        startDel = endDel = startIns = endIns = -1;

        while(i<this.length) {
            assert (0 == prevIns || 0 == prevDel);

            if(Alignment.GAP == this.read[i]) {
                if(0 == prevDel) {
                    startDel = i;
                }
                prevDel = 1;
                endDel = i;
                prevIns = 0;
                startIns = -1;
                endIns = -1;
                i++;
            }
            else if(Alignment.GAP == this.reference[i]) {
                if(0 == prevIns) {
                    startIns = i;
                }
                prevIns = 1;
                endIns = i;
                prevDel = 0;
                startDel = -1;
                endDel = -1;
                i++;
            }
            else {
                if(1 == prevDel) {
                    assert (0 < startDel);
                    assert (startDel <= endDel);
                    startDel--;
                    while(0 <= startDel && // Bases remaining to examine 
                            this.read[startDel] != Alignment.GAP && // Hit another deletion 
                            this.reference[startDel] != Alignment.GAP && // Hit an insertion 
                            this.reference[startDel] == this.reference[endDel]) { // src ref base matches dest ref base 
                        assert (Alignment.GAP != this.reference[startDel]);
                        assert (Alignment.GAP != this.reference[endDel]);
                        assert (Alignment.GAP != this.read[startDel]);
                        assert (Alignment.GAP == this.read[endDel]);
                        this.read[endDel] = this.read[startDel];
                        this.read[startDel] = Alignment.GAP;
                        startDel--;
                        endDel--;
                            }
                    endDel++; // We decremented when we exited the loop 
                    i = endDel;
                    assert (Alignment.GAP != this.read[i]);
                    assert (Alignment.GAP != this.reference[i]);
                }
                else if(1 == prevIns) {
                    assert (startIns <= endIns);
                    startIns--;
                    while(0 <= startIns && // Bases remaining to examine 
                            this.read[startIns] != Alignment.GAP && // Hit another deletion 
                            this.reference[startIns] != Alignment.GAP && // Hit an insertion 
                            this.read[startIns] == this.read[endIns]) { // src this.read base matches dest this.read base 
                        assert (Alignment.GAP != this.read[startIns]);
                        assert (Alignment.GAP != this.read[endIns]);
                        assert (Alignment.GAP != this.reference[startIns]);
                        assert (Alignment.GAP == this.reference[endIns]);
                        this.reference[endIns] = this.reference[startIns];
                        this.reference[startIns] = Alignment.GAP;
                        startIns--;
                        endIns--;
                            }
                    endIns++; // We decremented when we exited the loop 
                    i = endIns;
                    assert (Alignment.GAP != this.read[i]);
                    assert (Alignment.GAP != this.reference[i]);
                }
                else {
                    i++;
                }
                prevDel = 0;
                prevIns = 0;
                startDel = -1;
                endDel = -1;
                startIns = -1;
                endIns = -1;
            }
        }
    }
}
