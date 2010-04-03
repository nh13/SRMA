/*
 * LICENSE to be determined
 */
package srma;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;

import java.util.*;
import java.io.*;
import srma.Node;

public class Graph {
    int contig; // one based
    int position_start; // one based
    int position_end; // one based
    ArrayList<PriorityQueue<Node>> nodes; // zero based
    ArrayList<Integer> coverage; // does not count insertions with offset > 0
    SAMFileHeader header;

    public Graph(SAMFileHeader header)
    {
        this.header = header;
        this.contig = 1; 
        this.position_start = 1; 
        this.position_end = 1;
        this.nodes = new ArrayList<PriorityQueue<Node>>(); 
        this.coverage = new ArrayList<Integer>();
    }

    // Returns start/end node in the alignment graph with respect to strand
    public Node addSAMRecord(SAMRecord record, ReferenceSequence sequence) throws Exception
    {
        Alignment alignment;
        PriorityQueue<Node> nodeQueue = null;
        int i, ref_i, offset, node_type;
        Node prev=null, cur=null, ret=null;
        boolean strand = false;

        // Get the alignment
        alignment = new Alignment(record, sequence);
        strand = record.getReadNegativeStrandFlag(); 

        // Reset if there are no nodes
        if(0 == this.nodes.size()) {
            this.position_start = record.getAlignmentStart();
            if(Alignment.GAP == alignment.reference[0]) { // insertion
                this.position_start--;
            }
            // TOD0: could be insertions then deletions at the start, which will cause errors, not implemented yet
            this.position_end = this.position_start;
            this.contig = record.getReferenceIndex() + 1;
        }

        // HERE
        /*
        System.err.println(record.toString()); 
        alignment.print(System.err);
        System.err.println("HERE: " + record.getAlignmentStart() + "-" + record.getAlignmentEnd());
        */

        /* Reminders:
           i - index from 0 to 'alignment.length' 
           ref_i - index within 'alignment.reference'
           */
            
        for(i=0,ref_i=-1;i<alignment.length;i++,prev=cur) 
        { // go through the alignment

            // Skip over a deletion
            while(Alignment.GAP == alignment.read[i]) { 
                i++;
                ref_i++;
            }

            // Get the node type
            offset = 0;
            if(alignment.read[i] == alignment.reference[i]) { // match
                node_type = Node.MATCH;
            }
            else if(alignment.reference[i] == Alignment.GAP) { // insertion
                node_type = Node.INSERTION; 
                if(null != prev && Node.INSERTION == prev.type) {
                    offset = prev.offset + 1;
                }
            }
            else { // mismatch
                node_type = Node.MISMATCH;
            }
            if(null == prev || Node.INSERTION != prev.type) { // previous was an insertion, already on the position
                ref_i++;
            }

            // Create the node
            cur = this.addNode(new Node((char)alignment.read[i], 
                        node_type,
                        record.getReferenceIndex() + 1,
                        record.getAlignmentStart() + ref_i,
                        offset,
                        prev),
                    prev);

            // save return node
            if(null == prev && !strand) { // first node and forward strand
                ret = cur;
            }
        }

        if(strand) { // negative strand
            ret = cur;
        }

        return ret;
    }

    /* 
     * Adds the node to the graph.  Merges if the graph already
     * contains a similar node.
     * */
    private Node addNode(Node node, Node prev)
        throws Exception
    {

        Node curNode = null;
        int i;

        // Check if such a node exists
        // - if such a node exists, return it
        // - else insert it

        curNode = this.contains(node);
        if(null == curNode) { // new node, "how exciting!"
            if(node.contig != this.contig) { // destroy
                this.destroy();
                this.contig = node.contig;
                this.position_start = this.position_end = node.position;
            }
            // Add new queues if necessary
            for(i=this.position_end;i<=node.position;i++) {
                this.nodes.add(new PriorityQueue<Node>(1, new NodeComparator()));
                this.coverage.add(new Integer(0));
            }
            // Get the proper queue and add
            this.nodes.get(node.position - this.position_start).add(node);
            if(node.offset == 0) { // do not include insertions that extend an insertion
                this.coverage.set(node.position - this.position_start, node.coverage + this.coverage.get(node.position - this.position_start)); // set coverage
            }
            this.position_end = node.position;
            curNode = node;
        }
        else { // already contains
            curNode.coverage++; 
        }
        // Update edges
        if(null != prev) {
            curNode.addToPrev(prev);
            prev.addToNext(curNode);
        }

        return curNode;
    }

    /* 
     * Returns the Node in the graph if already exists,
     * null otherwise
     * */
    private Node contains(Node node)
    {
        PriorityQueue<Node> nodeQueue = null;
        Iterator<Node> nodeQueueIter = null;
        NodeComparator nodeComparator = null;
        Node curNode = null;

        // See if there are any nodes at this position
        try {
            nodeQueue = this.nodes.get(node.position - this.position_start);
        } catch (IndexOutOfBoundsException e) {
            return null;
        }

        // Go through all nodes at this position etc.
        nodeQueueIter = nodeQueue.iterator();
        nodeComparator = new NodeComparator();
        while(nodeQueueIter.hasNext()) {
            curNode = nodeQueueIter.next();
            if(nodeComparator.equals(curNode, node)) {
                return curNode;
            }
        }

        return null;
    }

    public int getPriorityQueueIndexAtPosition(int position)
    {
        PriorityQueue<Node> nodeQueue = null;

        if(position < this.position_start || this.position_end < position) {
            return 0;
        }

        nodeQueue = this.nodes.get(position - this.position_start);
        if(0 < nodeQueue.size()) {
            return position;
        }

        return 0;
    }

    public int getPriorityQueueIndexAtPositionOrGreater(int position)
        throws Exception
    {
        PriorityQueue<Node> nodeQueue = null;

        if(position < this.position_start) {
            position = this.position_start;
        }

        while(position <= this.position_end) {
            nodeQueue = this.nodes.get(position - this.position_start);
            if(0 < nodeQueue.size()) {
                return position;
            }
            position++;
        }

        return 0;
    }

    public int getPriorityQueueIndexAtPositionOrBefore(int position)
        throws Exception
    {
        PriorityQueue<Node> nodeQueue = null;

        if(this.position_end < position) {
            position = this.position_end;
        }

        while(this.position_start <= position) {
            nodeQueue = this.nodes.get(position - this.position_start);
            if(0 < nodeQueue.size()) {
                return position;
            }
            position--;
        }

        return 0;
    }

    public PriorityQueue<Node> getPriorityQueue(int position)
    {
        try {
            return this.nodes.get(position - this.position_start);
        } catch (IndexOutOfBoundsException e) {
            return null;
        }
    }

    public int getCoverage(int position)
    {
        try {
            return this.coverage.get(position - this.position_start);
        } catch (IndexOutOfBoundsException e) {
            return 0;
        }
    }

    public void prune(int referenceIndex, int alignmentStart, int offset)
        throws Exception
    {
        PriorityQueue<Node> queue = null;

        if(this.contig != referenceIndex+1) {
            throw new Exception("Pruning expects the same contig");
        }

        while(this.position_start < alignmentStart - offset) {
            // remove nodes from the queue
            queue = this.nodes.get(0); 
            if(null != queue) {
                while(null != queue.peek()) {
                    queue.poll().destroy();
                }
            }
            // destroy the first node in the queue
            queue = null;
            this.nodes.remove(0); 
            this.coverage.remove(0);
            this.position_start++;
        }
        // update position_start further
        while(this.position_start < alignmentStart) {
            if(0 < this.position_start) {
                if(null != this.nodes.get(0)) {
                    break;
                }
                this.nodes.remove(0); // remove empty
                this.coverage.remove(0);
            }
            this.position_start++;
        }
    }

    private void destroy()
    {
        int i;
        PriorityQueue<Node> queue;

        for(i=0;i<this.nodes.size();i++) {
            queue = this.nodes.get(i);
            while(null != queue.peek()) {
                queue.poll().destroy();
            }
        }

        this.contig = 1;
        this.position_start = this.position_end = 1;
        this.nodes.clear();
        this.coverage.clear();
    }

    public void print()
        throws Exception
    {
        this.print(System.out);
    }

    public void print(PrintStream out)
        throws Exception
    {
        int i;
        PriorityQueue<Node> queue;
        Iterator<Node> iter;

        out.println((1+contig)+":"+position_start+"-"+position_end);

        for(i=0;i<this.nodes.size();i++) {
            queue = this.nodes.get(i);
            iter = queue.iterator();
            while(iter.hasNext()) {
                iter.next().print(out);
            }
        }
    }
}
