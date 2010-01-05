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

    byte base; // [acgtnACGTN]
    int type;
    int contig;
    int position; // one-based
    int offset; // for insertions
    int coverage;
    List<Node> next; // downstream nodes
    List<Node> prev; // upstream nodes

    public Node(byte base, int type, int contig, int position, int offset, Node prev)
    {
        this.base = base;
        this.type = type;
        this.contig = contig;
        this.position = position;
        this.offset = offset;
        this.coverage = 1;
        this.next = new ArrayList<Node>();
        this.prev = new ArrayList<Node>();
        addToPrev(prev);
    }

    public void addToNext(Node node)
    {
        if(!this.next.contains(node)) {
            this.next.add(node);
        }
    }

    public void addToPrev(Node node) 
    {
        if(!this.prev.contains(node)) {
            this.prev.add(node);
        }
    }

    public void print(PrintStream out) 
    {
        ListIterator<Node> iter;
        out.print("[" + (char)this.base + ":" + this.type + ":" + this.contig + ":" + 
                this.position + ":" + this.offset + ":" + this.coverage + "]");
        iter = this.next.listIterator();
        while(iter.hasNext()) {
            Node next = iter.next();
            out.print("\t" + (char)next.base + ":" + next.type + ":" + next.contig + 
                    ":" + next.position + ":" + next.offset + ":" + next.coverage); 
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
