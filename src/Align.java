/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class Align {

    public static SAMRecord Align(Graph graph, SAMRecord rec, int offset, int coverage)
        throws Exception
    {
        if(rec.getReadNegativeStrandFlag()) {
            // Reverse
            throw new Exception("Error.  Reverse alignments not implemented.");
        }
        else {
            // Forward
            return AlignForwardStrand(graph, rec, offset, coverage);
        }
    }

    private static SAMRecord AlignForwardStrand(Graph graph, SAMRecord rec, int offset, int coverage)
        throws Exception
    {
        int i;
        PriorityQueue<Node> startNodeQueue = null;
        Iterator<Node> startNodeQueueIter = null;
        Node startNode;
        Node nextNode=null;
        AlignHeapNode curAlignHeapNode=null;
        AlignHeapNode nextAlignHeapNode=null;
        AlignHeapNode bestAlignHeapNode=null;
        AlignHeap heap=null;
        String read=null;
        String qualities=null;
        SRMAUtil.Space space=SRMAUtil.Space.NTSPACE;
        ListIterator<Node> iter=null;
        ListIterator<Integer> iterCov=null;
        AlignHeapNodeComparator comp;
        int alignmentStart = rec.getAlignmentStart();

        assert SRMAUtil.Space.COLORSPACE != space;

        // TODO:
        // - reverse strand
        // - get first node based on ???
        // - get colors and color qualities for COLORSPACE, must normalize colors to adaptor too
        // - recover alignment and print
        // - remove/fix paired end reads?

        // HERE
        System.err.println("ALI:"+rec.toString());

        read = (String)rec.getAttribute("CS");
        if(null == read) {
            space = SRMAUtil.Space.NTSPACE;
            read = new String(rec.getReadBases());
            if(read.length() <= 0) {
                throw new Exception("Error.  The current alignment has no bases.");
            }
            qualities = new String(rec.getBaseQualities());
            if(qualities.length() <= 0) {
                throw new Exception("Error.  The current alignment has no qualities.");
            }
        }
        else {
            read = SRMAUtil.normalizeColorSpaceRead(read);
            space = SRMAUtil.Space.COLORSPACE;
            qualities = (String)rec.getAttribute("CQ");
            if(null == qualities) {
                throw new Exception("Error.  The current color space alignment has no color qualities.");
            }
        }

        heap = new AlignHeap(AlignHeap.HeapType.MINHEAP); // should be set based on strand
        comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MINHEAP); // should be set based on strand

        for(i=alignmentStart-offset;i<=alignmentStart+offset;) {
            startNodeQueue = graph.getPriorityQueueAtPositionOrGreater(i);
            startNodeQueueIter = startNodeQueue.iterator();
            while(startNodeQueueIter.hasNext()) {
                startNode = startNodeQueueIter.next();
                heap.add(new AlignHeapNode(null, 
                            startNode,
                            startNode.coverage,
                            read.charAt(0),
                            qualities.charAt(0),
                            space));
                i=startNode.position+1;
            }
        }

        while(null != heap.peek()) {
            curAlignHeapNode = heap.poll();

            // HERE
            /*
               System.err.println("size:" + heap.size() + "\talignmentStart:" + alignmentStart + "\toffset:" + offset + "\treadOffset:" + curAlignHeapNode.readOffset);
               curAlignHeapNode.node.print(System.err);
               */

            // Remove all non-insertions with the same contig/pos/read-offset/type/base and lower score 
            nextAlignHeapNode = heap.peek();
            while(Node.INSERTION != curAlignHeapNode.node.type 
                    && null != nextAlignHeapNode 
                    && 0 == comp.compare(curAlignHeapNode, nextAlignHeapNode)) 
            {
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
            if(curAlignHeapNode.readOffset == read.length() - 1) {
                // Complete, store if has the best alignment.
                // TODO: must guarantee left-most alignment based on strand
                
                // HERE
                System.err.print(curAlignHeapNode.coverageSum + ":" + curAlignHeapNode.score + ":");
                curAlignHeapNode.node.print(System.err);
           
                if(null == bestAlignHeapNode 
                        || bestAlignHeapNode.score < curAlignHeapNode.score 
                        || (bestAlignHeapNode.score == curAlignHeapNode.score 
                            && bestAlignHeapNode.coverageSum < curAlignHeapNode.coverageSum)) 
                {
                    bestAlignHeapNode = curAlignHeapNode;
                }
            }
            else {
                // Based on strand etc.
                // Assuming only forward for now

                // Go to all "next" nodes
                iter = curAlignHeapNode.node.next.listIterator();
                iterCov = curAlignHeapNode.node.nextCov.listIterator();
                while(iter.hasNext()) {
                    nextNode = iter.next();
                    int nextCoverage = iterCov.next();
                    if(coverage <= nextCoverage) {
                        heap.add(new AlignHeapNode(curAlignHeapNode, 
                                    nextNode, 
                                    nextCoverage,
                                    read.charAt(curAlignHeapNode.readOffset+1), 
                                    qualities.charAt(curAlignHeapNode.readOffset+1), 
                                    space));
                    }
                }
                iter=null;
            }
        }

        // Recover alignment
        Align.updateSAM(rec, bestAlignHeapNode, space, read);

        return rec;
    }

    private static void updateSAM(SAMRecord rec, AlignHeapNode bestAlignHeapNode, SRMAUtil.Space space, String read)
        throws Exception
    {
        AlignHeapNode curAlignHeapNode;

        int alignmentStart = 0;
        int readIndex=-1;
        byte decodedBases[];
        byte colorErrors[];
        int i;
        int length=0;// HERE
        int length2=0;// HERE

        if(null == bestAlignHeapNode) {
            return;
        }

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

        decodedBases = new byte[read.length()];
        readIndex = read.length()-1;
        cigarElements = new LinkedList<CigarElement>();
        alignmentStart=bestAlignHeapNode.startPosition;

        // TODO
        assert null != bestAlignHeapNode;
        curAlignHeapNode = bestAlignHeapNode;
        while(null != curAlignHeapNode) {
            length2++;
            // Do stuff
            switch(curAlignHeapNode.node.type) {
                case Node.MATCH: /* Fall through */
                case Node.MISMATCH:
                    if(null != prevCigarOperator && null != curAlignHeapNode.prev && curAlignHeapNode.prev.node.position != curAlignHeapNode.node.position-1) { 
                        //System.out.println("DEL");
                        length+=prevCigarOperatorLength;
                        //if(null == prevCigarOperator) { throw new Exception("HERE 111"); }
                        cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                        cigarElements.add(0, new CigarElement(curAlignHeapNode.node.position - curAlignHeapNode.prev.node.position - 1, CigarOperator.DELETION));
                        // Start mismatch
                        prevCigarOperator = null;
                    }
                    else {
                        //System.out.println("M/MM");
                    }
                    if(space == SRMAUtil.Space.COLORSPACE) {
                        decodedBases[readIndex]  = (byte)curAlignHeapNode.node.base;
                        readIndex--;
                    }
                    curCigarOperator = CigarOperator.MATCH_OR_MISMATCH;
                    break;
                case Node.INSERTION:
                    //System.out.println("INS");
                    if(space == SRMAUtil.Space.COLORSPACE) {
                        decodedBases[readIndex]  = (byte)curAlignHeapNode.node.base;
                        readIndex--;
                    }
                    curCigarOperator = CigarOperator.INSERTION;
                    break;
                default:
                    throw new Exception("Unknown Type");
            }
            if(null == prevCigarOperator) {
                // first cigar operator
                prevCigarOperator = curCigarOperator;
                prevCigarOperatorLength=1;
            }
            else if(0 < prevCigarOperatorLength && prevCigarOperator != curCigarOperator) {
                // new cigar operator
                length += prevCigarOperatorLength;
                if(null == prevCigarOperator) { throw new Exception("HERE 222"); }
                cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                prevCigarOperator = curCigarOperator;
                prevCigarOperatorLength=1;
            }
            else {
                // same cigar operator
                prevCigarOperatorLength++;
            }

            // Update
            curAlignHeapNode = curAlignHeapNode.prev;
        }
        if(0 < prevCigarOperatorLength) {
            if(null == prevCigarOperator) { throw new Exception("HERE 333:" + prevCigarOperatorLength); }
            length += prevCigarOperatorLength;
            cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
        }

        // Get color error string
        if(space == SRMAUtil.Space.COLORSPACE) {
            colorErrors = new byte[read.length()];
            char prevBase = SRMAUtil.COLORSPACE_ADAPTOR;
            for(i=0;i<read.length();i++) {
                if(SRMAUtil.colorSpaceNextBase(prevBase, read.charAt(i)) == decodedBases[i]) {
                    colorErrors[i] = Alignment.GAP;
                }
                else {
                    colorErrors[i] = (byte)read.charAt(i);
                }
            }
        }
        else {
            decodedBases = new byte[read.length()];
            for(i=0;i<read.length();i++) {
                decodedBases[i] = (byte)read.charAt(i);
            }
        }

        // Update SAM record
        rec.setCigar(new Cigar(cigarElements));
        //System.out.println("rec.getAlignmentStart()="+rec.getAlignmentStart());
        //System.out.println("alignmentStart="+alignmentStart);
        rec.setAlignmentStart(alignmentStart);
        rec.setReadBases(decodedBases);
        // Clear attributes
        // TODO: use the new picard clear attributes function
        List<SAMRecord.SAMTagAndValue> atts = rec.getAttributes();
        ListIterator<SAMRecord.SAMTagAndValue> attsIter = atts.listIterator();
        while(attsIter.hasNext()) {
            SAMRecord.SAMTagAndValue att = attsIter.next();
            rec.setAttribute(att.tag, null);
        }
        // Set new attributes
        rec.setAttribute("AS", bestAlignHeapNode.score);
        rec.setAttribute("XC", bestAlignHeapNode.coverageSum);
        // set the XE attribute for colorError string

        // HERE
        if(length != rec.getReadBases().length) {
            System.out.println("length="+length+"\tlength2="+length2+"\trec.getReadBases().length="+rec.getReadBases().length);
            System.out.println("CIGAR="+rec.getCigarString());
            System.out.println("readOffset="+bestAlignHeapNode.readOffset);
            throw new Exception("CIGAR length and sequence length were inconsistent");
        }
        assert null != curAlignHeapNode; 
    }
}
