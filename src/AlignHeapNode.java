/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class AlignHeapNode {
    // These should be in a util class
    public enum Space { NTSPACE, COLORSPACE }
    public static final byte COLORSPACE_ADAPTOR = 'A';
    public static final byte DNA[] = {'A', 'C', 'G', 'T', 'N'};

    Node node;
    int readOffset; // # of bases from the beginning of the read
    int score; // alignment score
    int startPosition; // zero based
    Space space;

    public AlignHeapNode(AlignHeapNode prev,
            Node curNode,
            byte base,
            byte qual,
            Space space) throws Exception 
    {
        this.node = curNode;
        this.space = space;

        if(null == prev) { // first base
            if(Space.COLORSPACE == space) {
                base = colorSpaceNextBase(COLORSPACE_ADAPTOR, base);
            }
            this.readOffset = 0;
            this.score = 0;
            this.startPosition = curNode.position;
        }
        else {
            assert(prev.space == space);
            if(Space.COLORSPACE == space) {
                base = colorSpaceNextBase(prev.node.base, base);
            }
            this.readOffset = prev.readOffset + 1;
            this.score = prev.score;
            this.startPosition = prev.startPosition;
        }

        this.score += (base == curNode.base) ? 0 : -1*CHAR2QUAL(qual); 
    }

    // FUNCTIONS BELOW SHOULD GO IN A UTIL CLASS

    public int CHAR2QUAL(byte qual)
    {
        return (((int)qual) - 33);
    }

    public byte colorSpaceNextBase(byte base, byte color) throws Exception
    {
        int start=0, by=0, result=0;

        switch(base) {
            case 'A':
                start=0; by=1; break;
            case 'C':
                start=1; by=-1; break;
            case 'G':
                start=2; by=1; break;
            case 'T':
                start=3; by=-1; break;
            case 'N':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        switch(color) {
            case '0':
                result = start; break;
            case '1':
                result = start + by; break;
            case '2':
                result = start + 2*by; break;
            case '3':
                result = start + 3*by; break;
            case '4':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        if(result < 0) {
            return DNA[4 - ( (-1*result) % 4)];
        }
        else {
            return DNA[result % 4];
        }
    }
}
