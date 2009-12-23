/*
 * LICENSE to be determined
 */
package srma;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.*;

import java.util.*;
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

        // System.out.println("refr:" + new String(alignment.reference));
        // System.out.println("read:" + new String(alignment.read));

        int i, j, ref_i, seq_i;
        Node prev=null, cur=null;

        if(this.position_start <= record.getAlignmentStart() && record.getAlignmentStart() <= this.position_end) {
            prev = referenceNodes.get(record.getAlignmentStart() - this.position_start);
        }

        for(i=ref_i=seq_i=0;i<alignment.length;i++) { // go through the alignment
            // TODO: lower/upper case?
            if(alignment.read[i] == alignment.reference[i]) { // match
                cur = insertMatch(alignment.read[i],
                        record.getReferenceIndex(),
                        record.getAlignmentStart() + ref_i,
                        prev);
                ref_i++; seq_i++;
            }
            else if(alignment.reference[i] == Alignment.GAP) { // insertion
                assert Alignment.GAP != alignment.read[i];
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
                    if(null == prev) {
                        // This could throw an exception if the vector does not have the prev base
                        prev = this.referenceNodes.get(record.getAlignmentStart() - position_start); 
                    }
                    cur = new Node(alignment.read[i],
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
                cur = new Node(alignment.read[i],
                            Node.MISMATCH,
                            record.getReferenceIndex(),
                            record.getAlignmentStart() + ref_i,
                            0,
                            prev);
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

    private Node insertMatch(byte base, int contig, int position, Node prev)
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
        System.out.println((1+contig)+":"+position_start+"-"+position_end);

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
            node.print();
            iter = node.next.listIterator();
            while(iter.hasNext()) {
                Node next = iter.next();
                if(!queue.contains(next)) {
                    queue.add(next);
                }
            }
        }
    }
}
