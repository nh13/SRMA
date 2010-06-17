package srma;

import java.util.*;

// Notes: optimized for use in SRMA.java
public class ThreadPoolLinkedList<E> {

    private int numThreads;
    private List<List<E>> lists;
    private int first;
    private int next;
    private int last;
    private int size;
    E lastE;

    public ThreadPoolLinkedList(int numThreads) 
    {
        int i;
        this.numThreads = numThreads;
        this.lists = new ArrayList<List<E>>();
        for(i=0;i<this.numThreads;i++) {
            this.lists.add(new LinkedList<E>());
        }
        this.size = 0;
        this.first = 0;
        this.next = 0;
        this.last = 0;
        this.lastE = null;
    }

    public void add(E e)
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
        this.lastE = e;
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
        return this.lastE;
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
            this.lastE = null;
        }
        else {
            this.first++;
            if(this.numThreads <= this.first) {
                this.first = 0;
            }
        }
        return this.lists.get(prevFirst).remove(0);
    }

    public ListIterator<E>listIterator(int threadID)
        throws Exception
    {
        if(threadID < 0 || this.numThreads <= threadID) {
            throw new Exception("threadID out of bounds");
        }
        return this.lists.get(threadID).listIterator();
    }
}
