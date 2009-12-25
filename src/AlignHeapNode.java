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
    public static final byte COLORS[] = {'0', '1', '2', '3', '4'};

    AlignHeapNode prev; // previous
    Node node;
    int readOffset; // # of bases from the beginning of the read
    int score; // alignment score
    int startPosition; // one-based
    int alignmentLength; // alignment length
    Space space;

    /* 
     * Creates a heap node
     * @param prev The previous node, null otherwise.
     * @param curNode The current node in the graph.
     * @param base The base (or color) in the read.
     * @param qual The base (or color) quality in the read.
     * @param space The space of the alignment.
     * */
    public AlignHeapNode(AlignHeapNode prev,
            Node curNode,
            byte base,
            byte qual,
            Space space) throws Exception 
    {
        if(null == curNode) {
            throw new Exception("Error.  curNode was null");
        }

        this.node = curNode;
        this.space = space;

        if(null == prev) { // first base
            if(Space.COLORSPACE == space) {
                base = colorSpaceNextBase(COLORSPACE_ADAPTOR, base);
            }
            this.readOffset = 0;
            this.score = 0;
            this.startPosition = curNode.position;
            this.alignmentLength = 1;
            this.prev = null;
        }
        else {
            assert(prev.space == space);
            if(Space.COLORSPACE == space) {
                base = colorSpaceNextBase(prev.node.base, base);
            }
            this.readOffset = prev.readOffset + 1;
            this.score = prev.score;
            this.startPosition = prev.startPosition;
            this.alignmentLength = prev.alignmentLength + 1;
            this.prev = prev;
        }

        this.score += (base == curNode.base) ? 0 : -1*CHAR2QUAL(qual); 
    }

    // FUNCTIONS BELOW SHOULD GO IN A UTIL CLASS

    public static int CHAR2QUAL(byte qual)
    {
        return (((int)qual) - 33);
    }

    public static byte colorSpaceNextBase(byte base, byte color) throws Exception
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
    
    public static byte colorSpaceEncode(byte b1, byte b2) throws Exception
    {
        int start=0, by=0, result=0;

        switch(b1) {
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

        switch(b2) {
            case 'A':
                result = start; break;
            case 'C':
                result = start + by; break;
            case 'G':
                result = start + 2*by; break;
            case 'T':
                result = start + 3*by; break;
            case 'N':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        if(result < 0) {
            return COLORS[4 - ( (-1*result) % 4)];
        }
        else {
            return COLORS[result % 4];
        }
    }

    public static void normalizeColorSpaceRead(byte read[])
        throws Exception
    {
        if(read.length < 2) {
            throw new Exception("Read was too short");
        }
        else if(read[0] != COLORSPACE_ADAPTOR) {
            byte base = colorSpaceNextBase(read[0], read[1]);
            read[0] = COLORSPACE_ADAPTOR;
            read[1] = colorSpaceEncode(read[0], base);
        }
        // Remove adaptor
        int i;
        byte tmp[] = new byte[read.length-1];
        for(i=0;i<read.length-1;i++) {
            tmp[i] = read[i+1];
        }
        read = tmp;
    }
}
