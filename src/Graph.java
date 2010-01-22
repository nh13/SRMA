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

        /*
           System.err.println(record.toString()); 
           System.err.println("refr:" + new String(alignment.reference));
           System.err.println("read:" + new String(alignment.read));
           */

        int i, j, ref_i, read_i, node_type;
        Node prev=null, cur=null;

        /* Reminders:
           i - index from 0 to 'alignment.length' 
           ref_i - index within 'alignment.reference'
           read_i = index within 'alignment.read'
           */

        for(i=ref_i=read_i=0;
                i<alignment.length;
                i++,prev=cur) 
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

            // Create the node
            cur = this.addNode(new Node((char)alignment.read[i], 
                        node_type,
                        record.getReferenceIndex(),
                        record.getAlignmentStart() + ref_i,
                        0,
                        prev),
                    prev);
        }
    }

    /* 
     * Adds the node to the graph.  Merges if the graph already
     * contains a similar node.
     * */
    private Node addNode(Node node, Node prev)
    {

        Node curNode = null;
        PriorityQueue<Node> nodeQueue = null;

        // Check if such a node exists
        // - if such a node exists, return it
        // - else insert it

        curNode = this.contains(node);
        if(null == curNode) { // new node, how exciting
            if(node.contig != this.contig) {
                this.destroy();
            }
            try {
                nodeQueue = this.nodes.get(node.position - this.position_start);
            } catch (ArrayIndexOutOfBoundsException e) {
                nodeQueue = new PriorityQueue<Node>(10, new NodeComparator());
                this.nodes.add(node.position - this.position_start, nodeQueue);
            }
            nodeQueue.add(node);
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

        while(null == nodeQueue && position <= this.position_end) {
            try {
                nodeQueue = this.nodes.get(position - this.position_start);
                return nodeQueue;
            } catch (ArrayIndexOutOfBoundsException e) {
                // ignore
            }
            position++;
        }

        throw new GraphException(GraphException.OTHER, "Could not find an adequate node to start re-alignment.");
    }
    
    public PriorityQueue<Node> getPriorityQueueAtPositionOrBefore(int position)
        throws Exception
    {
        PriorityQueue<Node> nodeQueue = null;

        while(null == nodeQueue && this.position_start <= position) {
            try {
                nodeQueue = this.nodes.get(position - this.position_start);
                return nodeQueue;
            } catch (ArrayIndexOutOfBoundsException e) {
                // ignore
            }
            position--;
        }

        throw new GraphException(GraphException.OTHER, "Could not find an adequate node to start re-alignment.");
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
    {
        this.print(System.out);
    }

    public void print(PrintStream out)
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
