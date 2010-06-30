package srma;

import java.util.*;
import net.sf.samtools.*;

// Notes: optimized for use in SRMA.java
public class ThreadPoolLinkedList {

    private LinkedList<AlignRecord> buffer; // internal use only

    public ThreadPoolLinkedList()
    {
        this.buffer = new LinkedList<AlignRecord>();
    }

    // Add to buffer
    public void add(AlignRecord e)
    {
        this.buffer.add(e);
    }

    public int size() 
    {
        return this.buffer.size();
    }

    public AlignRecord getFirst()
    {
        if(0 == this.size()) {
            return null;
        }
        return this.buffer.getFirst();
    }

    public AlignRecord getLast()
    {
        if(0 == this.size()) {
            return null;
        }
        return this.buffer.getLast();
    }

    public AlignRecord removeFirst()
    {

        if(0 == this.size()) {
            return null;
        }
        return this.buffer.removeFirst();
    }

    // Only returns items from the specified contig.  Stops when an entry is not found.
    public LinkedList<LinkedList<AlignRecord>> getThreadLists(int numThreads, int contig)
    {
        int i = 0;
        int size = this.size();
        LinkedList<LinkedList<AlignRecord>> threadLists = new LinkedList<LinkedList<AlignRecord>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<AlignRecord>());
        }

        i=0;
        while(0 != size) {
            AlignRecord rec = this.buffer.getFirst();
            if(rec.record.getReferenceIndex()+1 != contig) {
                break;
            }
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
    public LinkedList<LinkedList<AlignRecord>> getAlignRecordThreadLists(int numThreads, int contig, int alignmentStartUpperBound)
    {
        int i = 0;
        int size = this.size();
        LinkedList<LinkedList<AlignRecord>> threadLists = new LinkedList<LinkedList<AlignRecord>>();

        for(i=0;i<numThreads;i++) {
            threadLists.add(new LinkedList<AlignRecord>());
        }

        i=0;
        while(0 != size) {
            AlignRecord rec = this.buffer.getFirst();
            if(rec.record.getReferenceIndex()+1 != contig ||
                    alignmentStartUpperBound <= rec.record.getAlignmentEnd()) {
                break;
            }
            threadLists.get(i).add(this.buffer.removeFirst());
            size--;
            i++;
            if(numThreads <= i) {
                i=0;
            }
        }

        return threadLists;
    }
}
