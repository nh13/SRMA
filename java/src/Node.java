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

    byte inUse = 0;
    char base; // [acgtnACGTN]
    int type;
    int contig;
    int position; // one-based
    int offset; // for insertions
    int coverage;
    List<NodeRecord> next; // downstram nodes
    List<NodeRecord> prev; // upstream nodes

    public Node(char base, int type, int contig, int position, Node prev, NodeRecordComparator comparator)
        throws Exception
    {
        this.inUse = 1;
        this.base = base;
        this.type = type;
        this.contig = contig;
        this.position = position;
        this.offset = 0;
        this.coverage = 1;
        this.next = new ArrayList<NodeRecord>();
        this.prev = new ArrayList<NodeRecord>();
        if(null != prev) {
            if(Node.INSERTION == prev.type && Node.INSERTION == this.type) {
                this.offset = prev.offset + 1;
            }
        }
    }

    public void addToNext(Node node, NodeRecordComparator comparator)
        throws Exception
    {
        NodeRecord rec = null;
        int index;
        if(null == node) {
            throw new Exception("addToNext: node was null!");
        }
        rec = new NodeRecord(node, 1);
        index = Collections.binarySearch(this.next, rec, comparator);
        if(0 <= index) { // exists
            this.next.get(index).coverage++;
        }
        else {
            // add it at the insertion point
            this.next.add(-index - 1, rec);
        }
    }


    public void addToPrev(Node node, NodeRecordComparator comparator)
        throws Exception
    {
        NodeRecord rec = null;
        int index;
        if(null == node) {
            throw new Exception("addToPrev: node was null!");
        }
        rec = new NodeRecord(node, 1);
        index = Collections.binarySearch(this.prev, rec, comparator);
        if(0 <= index) { // exists
            this.prev.get(index).coverage++;
        }
        else {
            // add it at the insertion point
            this.prev.add(-index - 1, rec);
        }
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

    private void removeFromList(List<NodeRecord> list, Node node, NodeRecordComparator comparator)
        throws Exception
    {
        int index = Collections.binarySearch(list, new NodeRecord(node, 1), comparator);
        if(0 <= index) { // exists
            list.remove(index);
        }
        else {
            throw new Exception("Could not remove node");
        }
    }

    public void removeFromNext(Node node, NodeRecordComparator comparator)
        throws Exception
    {
        this.removeFromList(this.next, node, comparator);
    }

    public void removeFromPrev(Node node, NodeRecordComparator comparator)
        throws Exception
    {
        this.removeFromList(this.prev, node, comparator);
    }

    public void removeLinks(NodeRecordComparator comparator)
        throws Exception
    {
        ListIterator<NodeRecord> iter = null;

        // remove links in prev
        iter = this.prev.listIterator();
        while(iter.hasNext()) {
            NodeRecord nodeRecord = iter.next();
            // sever the link
            nodeRecord.node.removeFromNext(this, comparator);
            // remove the node
            iter.remove();
        }
        // remove links in next
        iter = this.next.listIterator();
        while(iter.hasNext()) {
            NodeRecord nodeRecord = iter.next();
            // sever the link
            nodeRecord.node.removeFromPrev(this, comparator);
            // remove the node
            iter.remove();
        }
    }

    public void destroy()
    {
        // Clear the previous list
        this.prev.clear();

        // Clear the next list
        this.next.clear();
    }
}
