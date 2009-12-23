/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class Align {

    public Align(Graph graph, SAMRecord rec, int offset)
        throws Exception
    {
        if(rec.getReadNegativeStrandFlag()) {
            // Reverse
            throw new Exception("Error.  Reverse alignments not implemented.");
        }
        else {
            // Forward
            AlignForwardStrand(graph, rec, offset);
        }
    }

    private void AlignForwardStrand(Graph graph, SAMRecord rec, int offset)
        throws Exception
    {
        Node startNode=null;
        Node nextNode=null;
        AlignHeapNode curAlignHeapNode=null;
        AlignHeapNode nextAlignHeapNode=null;
        AlignHeapNode bestAlignHeapNode=null;
        AlignHeap heap=null;
        byte readBases[]=null;
        byte readQualities[]=null;
        AlignHeapNode.Space space=AlignHeapNode.Space.NTSPACE;
        ListIterator<Node> iter=null;
        AlignHeapNodeComparator comp;
        int alignmentStart = rec.getAlignmentStart();

        assert AlignHeapNode.Space.COLORSPACE != space;
        // TODO:
        // - reverse strand
        // - infer COLORSPACE
        // - get first node based on ???
        // - get colors and color qualities for COLORSPACE, must normalize colors to adaptor too
        // - recover alignment and print

        readBases = rec.getReadBases();
        if(readBases.length <= 0) {
            throw new Exception("Error.  The current alignment has no bases.");
        }
        readQualities = rec.getBaseQualities();
        if(readQualities.length <= 0) {
            throw new Exception("Error.  The current alignment has no qualities.");
        }

        heap = new AlignHeap(AlignHeap.HeapType.MINHEAP); // should be set based on strand
        comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MINHEAP); // should be set based on strand

        // Add first node - should be set based on strand
        heap.add(new AlignHeapNode(null, 
                    graph.getReferenceNode(alignmentStart - offset),
                    readBases[0], 
                    readQualities[0], 
                    space));

        while(null != heap.peek()) {
            curAlignHeapNode = heap.poll();

            /*
            System.out.println("size:" + heap.size());
            curAlignHeapNode.node.print();
            */

            // Remove all non-insertions with the same contig/pos/read offset/type/base and lower score 
            nextAlignHeapNode = heap.peek();
            while(Node.INSERTION != curAlignHeapNode.node.type && null != nextAlignHeapNode && 0 == comp.compare(curAlignHeapNode, nextAlignHeapNode)) {
                if(curAlignHeapNode.score < nextAlignHeapNode.score) {
                    // Update current node
                    curAlignHeapNode = heap.poll();
                }
                else {
                    // Ignore next node
                    heap.poll();
                }
                nextAlignHeapNode = heap.peek();
            }
            nextAlignHeapNode=null;
            // Check if the alignment is complete
            if(curAlignHeapNode.readOffset == readBases.length - 1) {
                // Complete, store if has the best alignment.
                // TODO: must guarantee left-most alignment based on strand
                if(null == bestAlignHeapNode || bestAlignHeapNode.score < curAlignHeapNode.score) {
                    bestAlignHeapNode = curAlignHeapNode;
                }
            }
            else {
                // Based on strand etc.
                // Assuming only forward for now
                
                // Go to all "next" nodes
                iter = curAlignHeapNode.node.next.listIterator();
                while(iter.hasNext()) {
                    nextNode = iter.next();
                    heap.add(new AlignHeapNode(curAlignHeapNode, 
                                nextNode, 
                                readBases[curAlignHeapNode.readOffset+1], 
                                readQualities[curAlignHeapNode.readOffset+1], 
                                space));
                    // Add a start node
                    // TODO should be conditioned on strand
                    if(nextNode.position <= alignmentStart + offset) {
                        heap.add(new AlignHeapNode(null, 
                                    nextNode, 
                                    readBases[0], 
                                    readQualities[0], 
                                    space));
                    }
                }
                iter=null;
            }
        }
        
        // Recover alignment
        this.updateSAM(rec, bestAlignHeapNode);
    }

    private void updateSAM(SAMRecord rec, AlignHeapNode bestAlignHeapNode)
        throws Exception
    {
        AlignHeapNode curAlignHeapNode;

        int alignmentStart = 0;
        byte read[]=null;
        int readIndex=-1;
        byte readBases[]=null;
        
        // To generate a new CIGAR
        List<CigarElement> cigarElements=null;
        CigarOperator prevCigarOperator=null, curCigarOperator=null;
        int prevCigarOperatorLength=0;
        
        // TODO 
        // setAttributes
        // setInferredInsertSize (invalidates paired end reads)
        // setMappingQuality (?)
        // setFlag
        // setReadBases (for color space)
        //
        //
        // color errors
        // do not use readBases in color space

        readBases = rec.getReadBases();
        readIndex = bestAlignHeapNode.alignmentLength-1;
        read = new byte[readIndex+1];
        cigarElements = new LinkedList<CigarElement>();
        alignmentStart=bestAlignHeapNode.startPosition;

        // TODO
        assert null != bestAlignHeapNode;
        curAlignHeapNode = bestAlignHeapNode;
        while(null != curAlignHeapNode) {
            // Do stuff
            //curAlignHeapNode.node.print(System.out);
            switch(curAlignHeapNode.node.type) {
                case Node.MATCH: /* Fall through */
                case Node.MISMATCH:
                    if(null == curAlignHeapNode.prev || curAlignHeapNode.prev.node.position == curAlignHeapNode.node.position-1) { 
                        //System.out.println("M/MM");
                        read[readIndex]=readBases[curAlignHeapNode.readOffset];
                        curCigarOperator = CigarOperator.MATCH_OR_MISMATCH;
                        readIndex--;
                    }
                    else {
                        //System.out.println("DEL");
                        curCigarOperator = CigarOperator.DELETION;
                    }
                    break;
                case Node.INSERTION:
                    //System.out.println("INS");
                    read[readIndex]=readBases[curAlignHeapNode.readOffset];
                    readIndex--;
                    curCigarOperator = CigarOperator.INSERTION;
                    break;
                default:
                    throw new Exception("Unknown Type");
            }
            if(null == prevCigarOperator) {
                prevCigarOperator = curCigarOperator;
                prevCigarOperatorLength=1;
            }
            else if(0 < prevCigarOperatorLength && prevCigarOperator != curCigarOperator) {
                cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                prevCigarOperator = curCigarOperator;
                prevCigarOperatorLength=1;
            }
            else {
                prevCigarOperatorLength++;
            }

            // Update
            curAlignHeapNode = curAlignHeapNode.prev;
        }
        if(0 < prevCigarOperatorLength) {
            cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
        }

        // Update SAM record
        rec.setCigar(new Cigar(cigarElements));
        System.out.println("rec.getAlignmentStart()="+rec.getAlignmentStart());
        System.out.println("alignmentStart="+alignmentStart);
        rec.setAlignmentStart(alignmentStart);
        rec.setReadBases(read);
        //rec.setAttributes(null);

        assert null != curAlignHeapNode; 
    }
}
