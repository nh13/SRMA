/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class AlignHeapNode {
    AlignHeapNode prev; // previous
    Node node;
    int readOffset; // # of bases from the beginning of the read
    int score; // alignment score
    int startPosition; // one-based
    int alignmentLength; // alignment length
    SRMAUtil.Space space;

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
            SRMAUtil.Space space) throws Exception 
    {
        if(null == curNode) {
            throw new Exception("Error.  curNode was null");
        }

        this.node = curNode;
        this.space = space;

        if(null == prev) { // first base
            if(SRMAUtil.Space.COLORSPACE == space) {
                base = SRMAUtil.colorSpaceNextBase(SRMAUtil.COLORSPACE_ADAPTOR, base);
            }
            this.readOffset = 0;
            this.score = 0;
            this.startPosition = curNode.position;
            this.alignmentLength = 1;
            this.prev = null;
        }
        else {
            assert(prev.space == space);
            if(SRMAUtil.Space.COLORSPACE == space) {
                base = SRMAUtil.colorSpaceNextBase(prev.node.base, base);
            }
            if(prev.node.position + 1 == this.node.position || this.node.type == Node.INSERTION) {
                // match/mismatch or insertion
                this.readOffset = prev.readOffset + 1;
            }
            else {
                // Deletion
                this.readOffset = prev.readOffset;
            }
            this.score = prev.score;
            this.startPosition = prev.startPosition;
            this.alignmentLength = prev.alignmentLength + 1;
            this.prev = prev;
        }
        this.score += (base == curNode.base) ? 0 : -1*SRMAUtil.CHAR2QUAL(qual); 
    }
}
