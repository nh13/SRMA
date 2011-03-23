/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class AlignHeap {
    public enum HeapType { MINHEAP, MAXHEAP }
    private final static int INITIAL_SIZE = 1024;

    private HeapType type;
    private PriorityQueue<AlignHeapNode> queue;

    public AlignHeap(HeapType type)
    {
        this.queue = new PriorityQueue<AlignHeapNode>(INITIAL_SIZE, new AlignHeapNodeComparator(type));
        this.type = type;
    }

    public void add(AlignHeapNode heapNode) 
    {
        this.queue.add(heapNode);
    }

    public AlignHeapNode poll()
    {
        return this.queue.poll();
    }

    public AlignHeapNode peek()
    {
        return this.queue.peek();
    }
    
    public int size()
    {
        return this.queue.size();
    }
}
