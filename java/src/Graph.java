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
    List<PriorityQueue<Node>> nodes; // zero based
    List<Integer> coverage; // does not count insertions with offset > 0
    NodeComparator nodeComparator; 
    private boolean isEmpty;

    public Graph()
    {
        this.contig = 1; 
        this.position_start = 1; 
        this.position_end = 1;
        this.nodes = new ArrayList<PriorityQueue<Node>>(); 
        this.coverage = new ArrayList<Integer>();
        this.nodeComparator = new NodeComparator();
        // Add two initial dummy elements
        this.nodes.add(new PriorityQueue<Node>(1, this.nodeComparator));
        this.coverage.add(new Integer(0));
        this.isEmpty = true;
    }

    // Returns start/end node in the alignment graph with respect to strand
    public Node addSAMRecord(SAMRecord record, ReferenceSequence sequence) throws Exception
    {
        Alignment alignment;
        PriorityQueue<Node> nodeQueue = null;
        int i, ref_i, offset, node_type, alignment_start, alignment_reference_index;
        Node prev=null, cur=null, ret=null;
        boolean strand = false;
        
        // DEBUG
        //String readName = record.getReadName();

        alignment_start = record.getAlignmentStart();
        alignment_reference_index = record.getReferenceIndex();

        if(alignment_reference_index != sequence.getContigIndex()) {
            throw new Exception("SAMRecord contig does not match the current reference sequence contig");
        }

        // Get the alignment
        alignment = new Alignment(record, sequence);
        strand = record.getReadNegativeStrandFlag(); 

        synchronized (this) {
            //System.err.println("HERE " + alignment_start + " " + this.position_start + ":" + this.position_end);
            if(alignment_start < this.position_start) {
                for(i=alignment_start;i<this.position_start;i++) {
                    this.nodes.add(0, new PriorityQueue<Node>(1, this.nodeComparator));
                    this.coverage.add(0, new Integer(0));
                }
                this.position_start = alignment_start;
            }

            // Reset if there are no nodes
            if(this.isEmpty) {
                this.position_start = alignment_start;
                if(Alignment.GAP == alignment.reference[0]) { // insertion
                    this.position_start--;
                }
                // TODO: could be insertions then deletions at the start, which will cause errors, not implemented yet
                this.position_end = this.position_start;
                this.contig = alignment_reference_index + 1;
                this.nodes.clear();
                this.coverage.clear();
                this.nodes.add(new PriorityQueue<Node>(1, this.nodeComparator));
                this.coverage.add(new Integer(0));
                this.isEmpty = false;
            }
        }

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
            if(alignment.read[i] == alignment.reference[i]) { // match
                node_type = Node.MATCH;
            }
            else if(alignment.reference[i] == Alignment.GAP) { // insertion
                node_type = Node.INSERTION; 
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
                        alignment_reference_index + 1,
                        alignment_start + ref_i,
                        prev,
                        this.nodeComparator),
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
    private synchronized Node addNode(Node node, Node prev)
        throws Exception
    {
        Node curNode = null;
        int i;

        // Check if such a node exists
        // - if such a node exists, return it
        // - else insert it
        curNode = this.contains(node);
        if(null == curNode) { // new node, "how exciting!"
            if(node.contig != this.contig) { // same contig 
                throw new Exception("NOT IMPLEMENTED");
            }
            // Add new queues if necessary
            for(i=this.position_end;i<node.position;i++) {
                this.nodes.add(new PriorityQueue<Node>(1, this.nodeComparator));
                this.coverage.add(new Integer(0));
            }
            // Get the proper queue and add
            this.nodes.get(node.position - this.position_start).add(node);

            // do not include insertions that extend an insertion
            if(Node.INSERTION != node.type || 0 == node.offset) {
                this.coverage.set(node.position - this.position_start, node.coverage + this.coverage.get(node.position - this.position_start)); // set coverage
            }
            if(this.position_end < node.position) {
                this.position_end = node.position;
            }
            curNode = node;
            this.isEmpty = false;
        }
        else { // already contains
            curNode.coverage++; 
            // do not include insertions that extend an insertion
            if(Node.INSERTION != curNode.type || 0 != curNode.offset) {
                // increment coverage
                this.coverage.set(curNode.position - this.position_start, 1 + this.coverage.get(curNode.position - this.position_start)); 
            }
        }
        // Update edges
        if(null != prev) {
            curNode.addToPrev(prev, this.nodeComparator);
            prev.addToNext(curNode, this.nodeComparator);
        }

        return curNode;
    }

    /* 
     * Returns the Node in the graph if already exists,
     * null otherwise
     * */
    private Node contains(Node node)
        throws Exception
    {
        PriorityQueue<Node> nodeQueue = null;
        Iterator<Node> nodeQueueIter = null;
        Node curNode = null;

        // See if there are any nodes at this position
        if(node.position - this.position_start < 0 || this.nodes.size() <= node.position - this.position_start) {
            return null;
        }
        nodeQueue = this.nodes.get(node.position - this.position_start);

        // Go through all nodes at this position etc.
        nodeQueueIter = nodeQueue.iterator();
        while(nodeQueueIter.hasNext()) {
            curNode = nodeQueueIter.next();
            if(this.nodeComparator.equals(curNode, node)) {
                return curNode;
            }
        }

        return null;
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

    public synchronized void prune(int referenceIndex, int alignmentStart, int offset)
        throws Exception
    {
        boolean shouldClear = false;
        
        if(this.contig != referenceIndex+1) {
            shouldClear = true;
        }
        else {
            if(this.position_start < alignmentStart - offset) { // unreachable nodes at the start
                if(this.position_end < alignmentStart - offset) { // all are unreachable
                    shouldClear = true;
                }
                else {
                    this.nodes = this.nodes.subList(alignmentStart - offset - this.position_start, this.nodes.size());
                    this.coverage = this.coverage.subList(alignmentStart - offset - this.position_start, this.coverage.size());
                    this.position_start = alignmentStart - offset;
                }
            }
        }
        if(shouldClear) {
            this.nodes.clear();
            this.coverage.clear();
            this.contig = referenceIndex + 1;
            this.position_start = this.position_end = alignmentStart;
            this.nodes.add(new PriorityQueue<Node>(1, this.nodeComparator));
            this.coverage.add(new Integer(0));
            this.isEmpty = true;
        }
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

    // Debugging function
    public void check()
        throws Exception
    {
        if(0 < this.nodes.size()) {
            int i;
            for(i=this.position_start;i<=this.position_end;i++) {
                PriorityQueue<Node> q = this.nodes.get(i - this.position_start);
                Iterator<Node> iter = q.iterator();
                while(iter.hasNext()) {
                    Node n = iter.next();
                    if(i != n.position) {
                        System.err.println("i="+i+"\tn.position="+n.position);
                        throw new Exception("Inconsistent graph");
                    }
                }
            }
        }
    }

    // Debugging function
    public void printDebug()
        throws Exception
    {
        int i;
        System.err.println(this.contig + ":" + this.position_start + "-" + this.position_end + " " + this.nodes.size() + " " + this.coverage.size());

        for(i=0;i<this.coverage.size();i++) {
            Node prev = null;
            PriorityQueue<Node> q1 = this.nodes.get(i);
            // copy queue
            PriorityQueue<Node> q2 = new PriorityQueue<Node>(1, new NodeComparator());
            Iterator<Node> iter = q1.iterator();

            while(iter.hasNext()) {
                q2.add(iter.next());
            }

            System.err.println((i+1)+" "+this.coverage.get(i)+" ");
            while(0 != q2.size()) {
                Node n = q2.poll();
                n.print(System.err);
                n.checkList(n.prev.listIterator(), this.nodeComparator);
                n.checkList(n.next.listIterator(), this.nodeComparator);
                if(null != prev) {
                    int c = this.nodeComparator.compare(prev, n);
                    if(0 < c) {
                        throw new Exception("OUT OF ORDER");
                    }
                    //System.err.println("comparison="+c);
                }
                prev = n;
            }
        }
    }
}
