/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class AlignHeap {
    public enum HeapType { MINHEAP, MAXHEAP }

    private HeapType type;
    private PriorityQueue<AlignHeapNode> queue;

    public AlignHeap(HeapType type)
    {
        this.queue = new PriorityQueue<AlignHeapNode>(0, new AlignHeapNodeComparator(type));
        this.type = type;
    }

    public void insert(AlignHeapNode heapNode) 
    {
        // TODO
    }

    public AlignHeapNode remove()
    {
        // TODO
        return null;
    }

    public AlignHeapNode peek()
    {
        // TODO
        return null;
    }
}
