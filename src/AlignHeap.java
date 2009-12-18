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
        this.queue.add(heapNode);
    }

    public AlignHeapNode remove()
    {
        return this.queue.poll();
    }

    public AlignHeapNode peek()
    {
        return this.queue.peek();
    }
}
