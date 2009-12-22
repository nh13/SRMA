/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import srma.*;

public class AlignHeapNodeComparator implements Comparator {
    private AlignHeap.HeapType type;

    public AlignHeapNodeComparator(AlignHeap.HeapType type) {
        this.type = type;
    }

    public int compare(Object o1, Object o2) 
    {
        AlignHeapNode a = (AlignHeapNode)o1;
        AlignHeapNode b = (AlignHeapNode)o2;

        // compare: contig, position, read offset, type, and base
        if(null == a ||
                a.node.contig < b.node.contig ||
                (a.node.contig == b.node.contig &&
                 a.node.position < b.node.position) ||
                (a.node.contig == b.node.contig &&
                 a.node.position == b.node.position &&
                 a.readOffset < b.readOffset) ||
                (a.node.contig == b.node.contig &&
                 a.node.position == b.node.position &&
                 a.readOffset < b.readOffset &&
                 a.node.type < b.node.type) ||
                (a.node.contig == b.node.contig &&
                 a.node.position == b.node.position &&
                 a.readOffset < b.readOffset &&
                 a.node.type < b.node.type &&
                 a.node.base < b.node.base)) {
            return (AlignHeap.HeapType.MINHEAP == type) ? -1 : 1;
                 }
        else if(null != b &&
                a.node.contig == b.node.contig &&
                a.node.position == b.node.position &&
                a.readOffset == b.readOffset &&
                a.node.type == b.node.type &&
                a.node.base == b.node.base) {
            return 0;
                }
        else {
            return (AlignHeap.HeapType.MINHEAP == type) ? 1 : -1;
        }
    }
}
