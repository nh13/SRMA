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
    List<NodeRecord> next; // downstram nodes
    List<NodeRecord> prev; // upstream nodes

    public Node(char base, int type, int contig, int position, Node prev, NodeComparator nodeComparator)
        throws Exception
    {
        this.base = base;
        this.type = type;
        this.contig = contig;
        this.position = position;
        this.offset = 0;
        this.coverage = 1;
        this.next = new LinkedList<NodeRecord>();
        this.prev = new LinkedList<NodeRecord>();
        if(null != prev) {
            if(Node.INSERTION == prev.type && Node.INSERTION == this.type) {
                this.offset = prev.offset + 1;
            }
        }
    }

    public void addToNext(Node node, NodeComparator nodeComparator)
        throws Exception
    {
        if(null == node) {
            throw new Exception("addToNext: node was null!");
        }
        
        this.addToList(this.next.listIterator(), node, nodeComparator);
    }


    public void addToPrev(Node node, NodeComparator nodeComparator)
        throws Exception
    {
        if(null == node) {
            throw new Exception("addToPrev: node was null!");
        }
        this.addToList(this.prev.listIterator(), node, nodeComparator);
    }
    
    private void addToList(ListIterator<NodeRecord> iter, Node node, NodeComparator nodeComparator) 
    {

        while(iter.hasNext()) {
            NodeRecord rec = iter.next();
            Node curN = rec.node;
            Integer curI = rec.coverage;
            int comparison = nodeComparator.compare(curN, node);

            if(0 == comparison) {
                iter.set(new NodeRecord(curN, curI+1));
                return;
            }
            else if(0 < comparison) {
                // move back one
                iter.previous();
                break;
            }
        }
        iter.add(new NodeRecord(node, 1));
    }

    public void checkList(ListIterator<NodeRecord> iter, NodeComparator nodeComparator)
        throws Exception
    {
        NodeRecord prev=null, cur=null;

        while(iter.hasNext()) {
            cur = iter.next();
            if(null != prev) {
                int comparison = nodeComparator.compare(prev.node, cur.node);
                if(0 < comparison) {
                    throw new Exception("OUT OF ORDER");
                }

                prev = cur;
            }
        }
    }

    public void print(PrintStream out) 
        throws Exception
    {
        ListIterator<NodeRecord> iter = null;

        out.print("[" + this.base + ":" + this.type + ":" + this.contig + ":" + 
                this.position + ":" + this.offset + ":" + this.coverage + "]");
        
        out.print("\tPREV[");
        iter = this.prev.listIterator();
        while(iter.hasNext()) {
            NodeRecord next = iter.next();
            Node nextNode = next.node;
            Integer nextCov = next.coverage;
            out.print("\t" + nextNode.base + ":" + nextNode.type + ":" + nextNode.contig + 
                    ":" + nextNode.position + ":" + nextNode.offset + ":" + nextNode.coverage + ":" + nextCov); 
        }
        out.print("]");
        out.print("\tNEXT[");
        iter = this.next.listIterator();
        while(iter.hasNext()) {
            NodeRecord next = iter.next();
            Node nextNode = next.node;
            Integer nextCov = next.coverage;
            out.print("\t" + nextNode.base + ":" + nextNode.type + ":" + nextNode.contig + 
                    ":" + nextNode.position + ":" + nextNode.offset + ":" + nextNode.coverage + ":" + nextCov); 
        }
        out.print("]");
        out.println("");
    }

    public void print()
        throws Exception
    {
        this.print(System.out);
    }

    public void destroy()
    {
        // Clear the previous list
        this.prev.clear();

        // Clear the next list
        this.next.clear();
    }

    public class NodeRecord {
        public Node node;
        public int coverage;

        public NodeRecord(Node node, int coverage) 
        {
            this.node = node;
            this.coverage = coverage;
        }
    }
}
