/*
 * LICENSE to be determined
 */
package srma;

import java.io.*;
import java.util.*;
    
public class Node {

    public static final int MATCH       = 0;
    public static final int MISMATCH    = 1;
    public static final int INSERTION   = 2; // into the read
    public static final int DELETION    = 3; // from the read

    char base; // [acgtnACGTN]
    int type;
    int contig;
    int position; // one-based
    int offset; // for insertions
    int coverage;
    List<Node> next; // downstream nodes
    List<Integer> nextCov;
    List<Node> prev; // upstream nodes
    List<Integer> prevCov;

    public Node(char base, int type, int contig, int position, int offset, Node prev)
    {
        this.base = base;
        this.type = type;
        this.contig = contig;
        this.position = position;
        this.offset = offset;
        this.coverage = 1;
        this.next = new ArrayList<Node>();
        this.nextCov = new ArrayList<Integer>();
        this.prev = new ArrayList<Node>();
        this.prevCov = new ArrayList<Integer>();
        addToPrev(prev);
    }

    public void addToNext(Node node)
    {
        int indexOf = this.next.indexOf(node);
        if(indexOf < 0) {
            this.next.add(node);
            this.nextCov.add(1);
        }
        else {
            this.nextCov.add(indexOf, this.nextCov.remove(indexOf) + 1);
        }
    }

    public void addToPrev(Node node) 
    {
        int indexOf = this.prev.indexOf(node);
        if(indexOf < 0) {
            this.prev.add(node);
            this.prevCov.add(1);
        }
        else {
            this.prevCov.add(indexOf, this.prevCov.remove(indexOf) + 1);
        }
    }

    public void print(PrintStream out) 
    {
        ListIterator<Node> iter;
        ListIterator<Integer> iterCov;
        out.print("[" + (char)this.base + ":" + this.type + ":" + this.contig + ":" + 
                this.position + ":" + this.offset + ":" + this.coverage + "]");
        iter = this.next.listIterator();
        iterCov = this.nextCov.listIterator();

        while(iter.hasNext()) {
            Node next = iter.next();
            Integer nextCov = iterCov.next();
            out.print("\t" + (char)next.base + ":" + next.type + ":" + next.contig + 
                    ":" + next.position + ":" + next.offset + ":" + next.coverage + ":" + nextCov); 
        }
        out.println("");
    }

    public void print()
    {
        this.print(System.out);
    }

    public class NodeComparator implements Comparator<Node> {
        public int compare(Node o1, Node o2) 
        {
            if(o1.contig < o2.contig ||
                    (o1.contig == o2.contig && o1.position < o2.position) ||
                    (o1.contig == o2.contig && o1.position == o2.position && o1.offset < o2.offset) ||
                    (o1.contig == o2.contig && o1.position == o2.position && o1.offset == o2.offset && o1.type < o2.type) ||
                    (o1.contig == o2.contig && o1.position == o2.position && o1.offset == o2.offset && o1.type == o2.type && o1.base < o2.base)) {
                return -1;
            }
            else if(o1.contig == o2.contig && o1.position == o2.position && o1.offset == o2.offset && o1.type == o2.type && o1.base == o2.base) {
                return 0;
            }
            else {
                return 1;
            }
        }

        public boolean equals(Node o1, Node o2) 
        {
            return (0 == this.compare(o1, o2));
        }
    }
}
