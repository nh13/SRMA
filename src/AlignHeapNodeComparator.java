/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import srma.*;

public class AlignHeapNodeComparator implements Comparator<AlignHeapNode> {
    private AlignHeap.HeapType type;

    /*
     * MAX heap sorts by maximum genomic coordinate (use with '-' strand)
     * MIN heap sorts by minimum genomic coordinate (use with '+' strand)
     * */
    public AlignHeapNodeComparator(AlignHeap.HeapType type) {
        this.type = type;
    }

    public int compare(AlignHeapNode a, AlignHeapNode b)
    {
        // sort by:
        // - MIN/MAX genomic coordinate
        // - MIN read offset
        // - min node type
        // - min base
        // - min score

        // contig
        if(a.node.contig < b.node.contig) {
            return (AlignHeap.HeapType.MINHEAP == this.type) ? -1 : 1;
        }
        else if(a.node.contig > b.node.contig) {
            return (AlignHeap.HeapType.MINHEAP == this.type) ? 1 : -1;
        }
        // position
        if(a.node.position < b.node.position) {
            return (AlignHeap.HeapType.MINHEAP == this.type) ? -1 : 1;
        }
        else if(a.node.position > b.node.position) {
            return (AlignHeap.HeapType.MINHEAP == this.type) ? 1 : -1;
        }
        // readOffset
        if(a.readOffset < b.readOffset) {
            return -1;
        }
        else if(a.readOffset > b.readOffset) {
            return 1;
        }
        // type
        if(a.node.type < b.node.type) {
            return -1;
        }
        else if(a.node.type > b.node.type) {
            return 1;
        }
        // base
        if(a.node.base < b.node.base) {
            return -1;
        }
        else if(a.node.base > b.node.base) {
            return 1;
        }
        // score
        if(a.score < b.score) {
            return -1;
        }
        else if(a.score > b.score) {
            return 1;
        }
        // same
        return 0;
    }
}
