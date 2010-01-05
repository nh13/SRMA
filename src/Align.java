/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;
import srma.*;

public class Align {

    public static void Align(Graph graph, SAMRecord rec, int offset, int coverage)
        throws Exception
    {
        if(rec.getReadNegativeStrandFlag()) {
            // Reverse
            throw new Exception("Error.  Reverse alignments not implemented.");
        }
        else {
            // Forward
            AlignForwardStrand(graph, rec, offset, coverage);
        }
    }

    private static void AlignForwardStrand(Graph graph, SAMRecord rec, int offset, int coverage)
        throws Exception
    {
        Node startNode=null;
        Node nextNode=null;
        AlignHeapNode curAlignHeapNode=null;
        AlignHeapNode nextAlignHeapNode=null;
        AlignHeapNode bestAlignHeapNode=null;
        AlignHeap heap=null;
        byte read[]=null;
        byte qualities[]=null;
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

        read = (byte[])rec.getAttribute("CS");
        if(null == read) {
            space = SRMAUtil.Space.NTSPACE;
            read = rec.getReadBases();
            if(read.length <= 0) {
                throw new Exception("Error.  The current alignment has no bases.");
            }
            qualities = rec.getBaseQualities();
            if(qualities.length <= 0) {
                throw new Exception("Error.  The current alignment has no qualities.");
            }
        }
        else {
            SRMAUtil.normalizeColorSpaceRead(read);
            space = SRMAUtil.Space.COLORSPACE;
            qualities = (byte[])rec.getAttribute("CQ");
            if(null == qualities) {
                throw new Exception("Error.  The current color space alignment has no color qualities.");
            }
        }

        heap = new AlignHeap(AlignHeap.HeapType.MINHEAP); // should be set based on strand
        comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MINHEAP); // should be set based on strand

        // Add first node - should be set based on strand
        heap.add(new AlignHeapNode(null, 
                    graph.getReferenceNode(alignmentStart - offset),
                    read[0], 
                    qualities[0], 
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
            if(curAlignHeapNode.readOffset == read.length - 1) {
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
                iterCov = curAlignHeapNode.node.nextCov.listIterator();
                while(iter.hasNext()) {
                    nextNode = iter.next();
                    if(coverage <= iterCov.next()) {
                        heap.add(new AlignHeapNode(curAlignHeapNode, 
                                    nextNode, 
                                    read[curAlignHeapNode.readOffset+1], 
                                    qualities[curAlignHeapNode.readOffset+1], 
                                    space));
                        // Add a start node
                        // TODO should be conditioned on strand
                        if(nextNode.position <= alignmentStart + offset) {
                            heap.add(new AlignHeapNode(null, 
                                        nextNode, 
                                        read[0], 
                                        qualities[0], 
                                        space));
                        }
                    }
                }
                iter=null;
            }
        }

        // Recover alignment
        Align.updateSAM(rec, bestAlignHeapNode, space, read);
    }

    private static void updateSAM(SAMRecord rec, AlignHeapNode bestAlignHeapNode, SRMAUtil.Space space, byte read[])
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

        decodedBases = new byte[bestAlignHeapNode.alignmentLength];
        readIndex = bestAlignHeapNode.alignmentLength-1;
        cigarElements = new LinkedList<CigarElement>();
        alignmentStart=bestAlignHeapNode.startPosition;

        // TODO
        assert null != bestAlignHeapNode;
        curAlignHeapNode = bestAlignHeapNode;
        while(null != curAlignHeapNode) {
            length2++;
            // Do stuff
            //curAlignHeapNode.node.print(System.out);
            switch(curAlignHeapNode.node.type) {
                case Node.MATCH: /* Fall through */
                case Node.MISMATCH:
                    if(null != curAlignHeapNode.prev && curAlignHeapNode.prev.node.position != curAlignHeapNode.node.position-1) { 
                        //System.out.println("DEL");
                        length+=prevCigarOperatorLength;
                        cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                        cigarElements.add(0, new CigarElement(curAlignHeapNode.node.position - curAlignHeapNode.prev.node.position - 1, CigarOperator.DELETION));
                        // Start mismatch
                        prevCigarOperator = null;
                    }
                    else {
                        //System.out.println("M/MM");
                    }
                    if(space == SRMAUtil.Space.COLORSPACE) {
                        decodedBases[readIndex]  = curAlignHeapNode.node.base;
                        readIndex--;
                    }
                    curCigarOperator = CigarOperator.MATCH_OR_MISMATCH;
                    break;
                case Node.INSERTION:
                    //System.out.println("INS");
                    if(space == SRMAUtil.Space.COLORSPACE) {
                        decodedBases[readIndex]  = curAlignHeapNode.node.base;
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
            length += prevCigarOperatorLength;
            cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
        }

        // Get color error string
        if(space == SRMAUtil.Space.COLORSPACE) {
            colorErrors = new byte[bestAlignHeapNode.alignmentLength];
            byte prevBase = SRMAUtil.COLORSPACE_ADAPTOR;
            for(i=0;i<read.length;i++) {
                if(SRMAUtil.colorSpaceNextBase(prevBase, read[i]) == decodedBases[i]) {
                    colorErrors[i] = Alignment.GAP;
                }
                else {
                    colorErrors[i] = read[i];
                }
            }
        }

        // Update SAM record
        rec.setCigar(new Cigar(cigarElements));
        //System.out.println("rec.getAlignmentStart()="+rec.getAlignmentStart());
        //System.out.println("alignmentStart="+alignmentStart);
        rec.setAlignmentStart(alignmentStart);
        if(space == SRMAUtil.Space.COLORSPACE) {
            rec.setReadBases(decodedBases);
        }
        else {
            rec.setReadBases(read);
        }
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
