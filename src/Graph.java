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
    Vector<Node> referenceNodes;
    SAMFileHeader header;
    List<ReferenceSequence> sequences;

    public Graph(SAMFileHeader header, List<ReferenceSequence> sequences)
    {
        this.header = header;
        this.sequences = sequences;
        this.contig = 0; 
        this.position_start = -1; 
        this.position_end = -1;
        referenceNodes = new Vector<Node>(); 
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

        int i, j, ref_i, seq_i;
        Node prev=null, cur=null;

        if(this.position_start <= record.getAlignmentStart() && record.getAlignmentStart() <= this.position_end && alignment.reference[0] != Alignment.GAP) {
            prev = referenceNodes.get(record.getAlignmentStart() - this.position_start);
        }

        for(i=ref_i=seq_i=0;i<alignment.length;i++) { // go through the alignment
            // TODO: lower/upper case?
            if(alignment.read[i] == alignment.reference[i]) { // match
                cur = insertMatch((char)alignment.read[i],
                        record.getReferenceIndex(),
                        record.getAlignmentStart() + ref_i,
                        prev);
                ref_i++; seq_i++;
            }
            else if(alignment.reference[i] == Alignment.GAP) { // insertion
                assert Alignment.GAP != alignment.read[i];
                // begins with an insertion
                if(null == prev) { 
                    assert i == 0;
                    // See if such an insertion exists

                    // Get insertion range
                    int start=0, end=0, found=0;

                    for(end=0;end<alignment.length;end++) {
                        if(alignment.reference[end] != Alignment.GAP) {
                            break;
                        }
                    }

                    try {

                        Queue<Node> nodeQueue = new LinkedList<Node>();
                        Queue<Integer> intQueue = new LinkedList<Integer>(); 
                        // Get Node and try to see if there is an incoming insertion etc...
                        nodeQueue.add(this.referenceNodes.get(record.getAlignmentStart() - position_start));
                        intQueue.add(end);

                        while(0 < nodeQueue.size() && 0 == found) {
                            // Get element off the queue
                            Node nodeCur = nodeQueue.poll();
                            int endCur = intQueue.poll();

                            // Go through all previous
                            ListIterator<Node> nodeIter = nodeCur.prev.listIterator();
                            while(nodeIter.hasNext() && 0 == found) {
                                Node node = nodeIter.next();
                                // only accept insertions that have the potential to account for the above alignment
                                if(node.type == Node.INSERTION && node.offset <= endCur - start + 1 && alignment.read[endCur] == node.base) {
                                    nodeQueue.add(node);
                                    intQueue.add(endCur-1);
                                }
                                else if(endCur == start - 1) {
                                    prev = node;
                                    found = 1;
                                }
                            }
                        }

                    } catch (ArrayIndexOutOfBoundsException e) {
                        // unsuccessful retrieval of the reference nodes -> ignore
                    }
                    if(0 == found) {
                        // Must insert from the first reference base (ugh)
                        throw new GraphException(GraphException.NOT_IMPLEMENTED);
                    }
                }
                // check if the insertion exists
                for(cur=null,j=0;null != prev && j < prev.next.size();j++) {
                    Node tmpNode = prev.next.get(j);
                    if(Node.INSERTION == tmpNode.type &&
                            alignment.read[i] == tmpNode.base) {
                        cur = tmpNode;
                        cur.coverage++;
                        break; // Found one
                            }
                }
                if(null == cur) { // No such insertion
                    cur = new Node((char)alignment.read[i],
                            Node.INSERTION,
                            prev.contig,
                            prev.position,
                            prev.offset+1,
                            prev);
                }
                seq_i++;
            }
            else if(alignment.read[i] == Alignment.GAP) { // deletion
                assert Alignment.GAP != alignment.read[i];
                cur = prev;
                // Skip over deletion.  Note: impossible to end with a deletion.
                while(Alignment.GAP == alignment.read[i]) {
                    assert i < alignment.length;
                    assert Alignment.GAP != alignment.reference[i];
                    i++;
                    ref_i++;
                }
                i--;
            }
            else { // mismatch
                // TODO check mismatch does not exist 
                cur = new Node((char)alignment.read[i],
                        Node.MISMATCH,
                        record.getReferenceIndex(),
                        record.getAlignmentStart() + ref_i,
                        0,
                        prev);
                if(prev != null) {
                    ListIterator<Node> iter = prev.next.listIterator();
                    NodeComparator comp = new NodeComparator();
                    while(iter.hasNext()) {
                        Node node = iter.next();
                        if(0 == comp.compare(node, cur)) {
                            cur = node; 
                            break;
                        }
                    }
                }
                ref_i++; seq_i++;
            }
            if(null != prev 
                    && prev != cur) {
                prev.addToNext(cur);
                cur.addToPrev(prev);
                fixReverseInsertion(cur);
                    }
            prev = cur;
        }
    }

    public Node getReferenceNode(int position) 
        throws Exception
    {
        if(position < this.position_start) {
            position = this.position_start;
        }
        else if(this.position_end < position) {
            throw new Exception("Out of range");
        }
        return this.referenceNodes.get(position - this.position_start);
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
        Node node;
        for(i=0;i<this.referenceNodes.size();i++) {
            node = this.referenceNodes.get(i);
            destroyNode(node);
        }

        this.contig = 0;
        this.position_start = this.position_end = -1;
        this.referenceNodes.clear();
    }

    private void destroyNode(Node node)
    {
        ListIterator iter;
        Node n;

        // Clear the previous list
        iter = node.prev.listIterator();
        node.prev.clear();

        // Clear the next list
        iter = node.next.listIterator();
        while(iter.hasNext()) {
            n = (Node)iter.next();
            if(Node.MATCH != n.type) {
                destroyNode(n);
            }
        }
        node.next.clear();
    }

    private Node insertMatch(char base, int contig, int position, Node prev)
    {
        Node cur;

        if(contig != this.contig || this.position_end < position) { // Not within the range
            cur = new Node(base, Node.MATCH, contig, position, 0, prev);

            if(contig != this.contig) {
                this.destroy();
            }

            assert this.position_start <= position;

            if(this.referenceNodes.isEmpty()) {
                this.position_start = position;
                if(0 == this.referenceNodes.size()) {
                    this.referenceNodes.setSize(1);
                }
                this.referenceNodes.add(0, cur);
            }
            else {
                if(this.referenceNodes.size() < position - this.position_start + 1) {
                    this.referenceNodes.setSize(position - this.position_start + 1);
                }
                this.referenceNodes.add(position - this.position_start, cur);
            }
            this.contig = contig;
            this.position_end = position;
        }
        else {
            cur = this.referenceNodes.get(position - this.position_start);
            if(null == cur) {
                cur = new Node(base, Node.MATCH, contig, position, 0, prev);
                if(this.referenceNodes.size() < position - this.position_start + 1) {
                    this.referenceNodes.setSize(position - this.position_start + 1);
                }
                referenceNodes.set(position - this.position_start, cur);
            }
            cur.coverage++;
        }
        return cur;
    }

    private void fixReverseInsertion(Node node)
    {
        // TODO
        // Ignore for now since this is computationally expensive
    }

    public void print()
    {
        this.print(System.out);
    }

    public void print(PrintStream out)
    {
        out.println((1+contig)+":"+position_start+"-"+position_end);

        PriorityQueue<Node> queue;
        ListIterator<Node> iter;

        queue = new PriorityQueue<Node>(10, new NodeComparator());
        Node node = (Node)referenceNodes.get(0);
        if(null != node) {
            queue.add(node);
        }

        // BFS
        while(null != queue.peek()) {
            node = queue.poll();
            node.print(out);
            iter = node.next.listIterator();
            while(iter.hasNext()) {
                Node next = iter.next();
                if(!queue.contains(next)) {
                    queue.add(next);
                }
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
    }
}
