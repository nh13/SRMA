package srma;

import java.util.*;

// Notes: optimized for use in SRMA.java
public class ThreadPoolLinkedList<E> {

    private LinkedList<E> buffer; // internal use only

    public ThreadPoolLinkedList()
    {
        this.buffer = new LinkedList<E>();
    }

    // Add to buffer
    public void add(E e)
    {
        this.buffer.add(e);
    }

    public int size() 
    {
        return this.buffer.size();
    }

    public E getFirst()
    {
        if(0 == this.size()) {
            return null;
        }
        return this.buffer.getFirst();
    }

    public E getLast()
    {
        if(0 == this.size()) {
            return null;
        }
        return this.buffer.getLast();
    }

    public E removeFirst()
    {

        if(0 == this.size()) {
            return null;
        }
        return this.buffer.removeFirst();
    }

    public LinkedList<LinkedList<E>> getThreadLists(int numThreads)
    {
        int i = 0;
        int size = this.size();
        LinkedList<LinkedList<E>> threadLists = new LinkedList<LinkedList<E>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<E>());
        }

        i=0;
        while(0 != size) {
            threadLists.get(i).add(this.buffer.removeFirst());
            size--;
            i++;
            if(numThreads <= i) {
                i=0;
            }
        }

        return threadLists;
    }

    // Do not use this if we mean to flush
    public LinkedList<LinkedList<AlignRecord>> getAlignRecordThreadLists(int numThreads, int alignmentStartLowerBound)
    {
        int i = 0;
        int size = this.size();
        LinkedList<LinkedList<AlignRecord>> threadLists = new LinkedList<LinkedList<AlignRecord>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<AlignRecord>());
        }

        i=0;
        while(0 != size) {
            AlignRecord rec = (AlignRecord)this.buffer.getFirst();
            if(alignmentStartLowerBound <= rec.record.getAlignmentEnd()) {
                break;
            }
            threadLists.get(i).add((AlignRecord)this.buffer.removeFirst());
            size--;
            i++;
            if(numThreads <= i) {
                i=0;
            }
        }

        return threadLists;
    }

}
