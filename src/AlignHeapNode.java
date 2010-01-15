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
            char base,
            char qual,
            SRMAUtil.Space space) throws Exception 
    {
        if(null == curNode) {
            throw new Exception("Error.  curNode was null");
        }

        this.node = curNode;
        this.space = space;

        if(null == prev) { // first base
            if(SRMAUtil.Space.COLORSPACE == space) {
                base = SRMAUtil.colorSpaceNextBase(SRMAUtil.COLORSPACE_ADAPTOR, (char)base);
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
                base = SRMAUtil.colorSpaceNextBase((char)prev.node.base, (char)base);
            }
            this.readOffset = prev.readOffset + 1;
            this.score = prev.score;
            this.startPosition = prev.startPosition;
            this.alignmentLength = this.node.position - prev.node.position + this.node.offset - prev.node.offset;
            this.prev = prev;
        }
        // HERE
        // System.err.println("base="+(char)base+" curNode.base="+(char)curNode.base+" score="+((base == curNode.base) ? 0 : -1*SRMAUtil.CHAR2QUAL(qual)));
        this.score += (base == curNode.base) ? 0 : -1*SRMAUtil.CHAR2QUAL(qual); 
    }
}
