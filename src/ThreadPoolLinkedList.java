package srma;

import java.util.*;

// Notes: optimized for use in SRMA.java
public class ThreadPoolLinkedList<E> {

    private int numThreads;
    private List<LinkedList<E>> lists;
    private int first;
    private int next;
    private int last;
    private int size;

    public ThreadPoolLinkedList(int numThreads) 
    {
        int i;
        this.numThreads = numThreads;
        this.lists = new ArrayList<LinkedList<E>>();
        for(i=0;i<this.numThreads;i++) {
            this.lists.add(new LinkedList<E>());
        }
        this.size = 0;
        this.first = 0;
        this.next = 0;
        this.last = 0;
    }

    public synchronized void add(E e)
    {
        this.lists.get(this.next).add(e);
        this.next++;
        if(this.numThreads <= this.next) {
            this.next = 0;
        }
        this.last++;
        if(this.numThreads <= this.last) {
            this.last = 0;
        }
        this.size++;
    }

    public int size() 
    {
        return this.size;
    }

    public E getFirst()
    {
        if(0 == this.size) {
            return null;
        }

        return this.lists.get(this.first).get(0);
    }

    public E getLast()
    {
        return this.lists.get(this.last).getLast();
    }

    public E removeFirst()
    {
        int prevFirst = this.first;

        if(0 == this.size) {
            return null;
        }
        this.size--;
        if(0 == this.size) {
            this.size = 0;
            this.first = 0;
            this.next = 0;
            this.last = 0;
        }
        else {
            this.first++;
            if(this.numThreads <= this.first) {
                this.first = 0;
            }
        }
        return this.lists.get(prevFirst).remove(0);
    }

    // This could cause some inconsistency if it is used to modify the underlying lists.  
    // Ignore this warning for performance.
    public ListIterator<E> listIterator(int threadID)
        throws Exception
    {
        if(threadID < 0 || this.numThreads <= threadID) {
            throw new Exception("threadID out of bounds");
        }
        return this.lists.get(threadID).listIterator();
    }

    // This could cause some inconsistency if it is used to modify the underlying lists.  
    // Ignore this warning for performance.
    public List<E> getListForThread(int threadID)
        throws Exception
    {
        if(threadID < 0 || this.numThreads <= threadID) {
            throw new Exception("threadID out of bounds");
        }
        return this.lists.get(threadID);
    }

    // TEST this
    public void rebalance()
        throws Exception
    {
        int i, size;

        size = this.lists.get(0).size();
        for(i=2;i<this.numThreads;i++) {
            int tmpSize = this.lists.get(i).size();
            if(tmpSize + 1 != size && tmpSize != size && tmpSize - 1 != size) {
                throw new Exception("rebalance will lead to inconsistency");
            }
        }

        if(0 == size) {
            this.size = 0;
            this.first = 0;
            this.next = 0;
            this.last = 0;
        }
        else {
            this.first = this.last = 0;
            this.size = 0;
            size = this.lists.get(0).size();
            for(i=2;i<this.numThreads;i++) {
                this.size += this.lists.get(i).size();
                if(size < this.lists.get(i).size()) {
                    this.first = this.last = i;
                }
                if(size == this.lists.get(i).size()) {
                    this.last = i;
                }
            }

            this.next = this.last + 1;
            if(this.numThreads <= this.next) {
                this.next = 0;
            }
        }
    }
}
