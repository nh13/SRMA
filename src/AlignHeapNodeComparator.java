/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import srma.*;

public class AlignHeapNodeComparator implements Comparator {
    private AlignHeap.HeapType type;


    /*
     * MAX heap sorts by maximum genomic coordinate (use with '-' strand)
     * MIN heap sorts by minimum genomic coordinate (use with '+' strand)
     * */
    public AlignHeapNodeComparator(AlignHeap.HeapType type) {
        this.type = type;
    }

    public int compare(Object o1, Object o2) 
    {
        AlignHeapNode a = (AlignHeapNode)o1;
        AlignHeapNode b = (AlignHeapNode)o2;

        int coordinate;

        // sort by:
        // - MIN/MAX genomic coordinate
        // - MIN reaad offset
        // - min node type
        // - min base
        // - min score

        if(null == a ||
                a.node.contig < b.node.contig ||
                (a.node.contig == b.node.contig &&
                 a.node.position < b.node.position)) {
            coordinate = -1;
        }
        else if(null != b &&
                a.node.contig == b.node.contig &&
                a.node.position == b.node.position) {
            coordinate = 0;
        }
        else {
            coordinate = 1;
        }
        // check if we are to sort by maximum genomic coordinate
        if(AlignHeap.HeapType.MAXHEAP == this.type) {
            coordinate *= -1;
        }


        // compare: contig, position, read offset, type, and base
        if(coordinate < 0 ||
                (coordinate == 0 &&
                 a.readOffset < b.readOffset) ||
                (coordinate == 0 &&
                 a.readOffset == b.readOffset &&
                 a.node.type < b.node.type) ||
                (coordinate == 0 &&
                 a.readOffset == b.readOffset &&
                 a.node.type == b.node.type &&
                 a.node.base < b.node.base) ||
                (coordinate == 0 &&
                 a.readOffset == b.readOffset &&
                 a.node.type == b.node.type &&
                 a.node.base == b.node.base &&
                 a.score < b.score)) {
            return -1;
                 }
        else if(null != b &&
                coordinate == 0 &&
                a.readOffset == b.readOffset &&
                a.node.type == b.node.type &&
                a.node.base == b.node.base &&
                a.score == b.score) {
            return 0;
                }
        else {
            return 1;
        }
    }
}
