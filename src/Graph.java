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
    int position_end;
    Vector<PriorityQueue<Node>> nodes;
    SAMFileHeader header;
    List<ReferenceSequence> sequences;

    public Graph(SAMFileHeader header, List<ReferenceSequence> sequences)
    {
        this.header = header;
        this.sequences = sequences;
        this.contig = 1; 
        this.position_start = 1; 
        this.position_end = 1;
        this.nodes = new Vector<PriorityQueue<Node>>(); 
    }

    public void addSAMRecord(SAMRecord record) throws Exception
    {
        Alignment alignment;

        // Get the alignment
        alignment = new Alignment(record, sequences);

        // HERE
        /*
           System.err.println(record.toString()); 
           System.err.println("refr:" + new String(alignment.reference));
           System.err.println("read:" + new String(alignment.read));
           */

        int i, j, ref_i, offset, node_type;
        Node prev=null, cur=null;

        /* Reminders:
           i - index from 0 to 'alignment.length' 
           ref_i - index within 'alignment.reference'
           */

        for(i=j=0,ref_i=-1;
                i<alignment.length;
                i++,prev=cur) 
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
                ref_i++;
            }
            else if(alignment.reference[i] == Alignment.GAP) { // insertion
                node_type = Node.INSERTION; 
                j=0; // not advancing the reference
                if(null != prev && Node.INSERTION == prev.type) {
                    offset = prev.offset + 1;
                }
            }
            else { // mismatch
                node_type = Node.MISMATCH;
                j=1;
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
        }
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
        if(null == curNode) { // new node, how exciting
            if(node.contig != this.contig) {
                this.destroy();
                this.contig = node.contig;
                this.position_start = this.position_end = node.position;
            }
            // Add new queues if necessary
            for(i=this.position_end;i<=node.position;i++) {
                this.nodes.add(new PriorityQueue<Node>(1, new NodeComparator()));
            }
            // Get the proper queue and add
            this.nodes.get(node.position - this.position_start).add(node);
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
        } catch (ArrayIndexOutOfBoundsException e) {
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

    public PriorityQueue<Node> getPriorityQueueAtPositionOrGreater(int position)
        throws Exception
    {
        PriorityQueue<Node> nodeQueue = null;

        if(position < this.position_start) {
            position = this.position_start;
        }

        while(position <= this.position_end) {
            nodeQueue = this.nodes.get(position - this.position_start);
            if(0 < nodeQueue.size()) {
                return nodeQueue;
            }
            position++;
        }

        throw new GraphException(GraphException.OTHER, "Could not find an adequate node to start re-alignment.");
    }

    public PriorityQueue<Node> getPriorityQueueAtPositionOrBefore(int position)
        throws Exception
    {
        PriorityQueue<Node> nodeQueue = null;
        int t = position; // HERE 

        if(this.position_end < position) {
            position = this.position_end;
        }

        while(this.position_start <= position) {
            nodeQueue = this.nodes.get(position - this.position_start);
            if(0 < nodeQueue.size()) {
                return nodeQueue;
            }
            position--;
        }

        return null;
    }

    public void prune(int start) 
        throws GraphException
    {
        // TODO
        throw new GraphException(GraphException.NOT_IMPLEMENTED);
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

    public class GraphException extends Exception
    {
        int type;
        public static final int NOT_IMPLEMENTED = 0;
        public static final int OTHER = 1;

        public GraphException(int type)
        {
            super("Not implemented");
            this.type = type;
        }

        public GraphException(int type, String message) {
            super(message);
            this.type = type;
        }
    }
}
