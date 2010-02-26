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

        int i;
        PriorityQueue<Node> startNodeQueue = null;
        Iterator<Node> startNodeQueueIter = null;
        Node startNode;
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
        int alignmentStart = -1;
        int numStartNodesAdded = 0;
        boolean strand = rec.getReadNegativeStrandFlag(); // false -> forward, true -> reverse

        // Debugging stuff
        String readName = rec.getReadName();

        assert SRMAUtil.Space.COLORSPACE != space;

        // TODO:
        // - remove/fix paired end reads?

        read = (String)rec.getAttribute("CS");
        if(null == read) {
            // Use base space
            space = SRMAUtil.Space.NTSPACE;
            if(strand) { // reverse
                byte tmp[] = rec.getReadBases();
                SAMRecordUtil.reverseArray(tmp);
                read = new String(tmp);
                SAMRecordUtil.reverseArray(tmp); // reverse back!
                tmp = rec.getBaseQualities();
                SAMRecordUtil.reverseArray(tmp);
                qualities = new String(tmp);
                SAMRecordUtil.reverseArray(tmp); // reverse back!
            }
            else { // forward
                read = new String(rec.getReadBases());
                qualities = new String(rec.getBaseQualities());
            }
            if(read.length() <= 0) {
                throw new Exception("Error.  The current alignment has no bases.");
            }
            if(qualities.length() <= 0) {
                throw new Exception("Error.  The current alignment has no qualities.");
            }
        }
        else {
            // assumes CS and CQ are always in sequencing order
            read = SRMAUtil.normalizeColorSpaceRead(read);
            space = SRMAUtil.Space.COLORSPACE;
            qualities = (String)rec.getAttribute("CQ");
            if(null == qualities) {
                throw new Exception("Error.  The current color space alignment has no color qualities.");
            }
        }

        heap = new AlignHeap(AlignHeap.HeapType.MINHEAP); // should be set based on strand
        if(strand) { // reverse
            comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MAXHEAP); 
        }
        else { // forward
            comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MINHEAP); 
        }

        // Add start nodes
        if(strand) { // reverse
            alignmentStart = rec.getAlignmentEnd();
            for(i=alignmentStart+offset;alignmentStart-offset<=i;) {
                startNodeQueue = graph.getPriorityQueueAtPositionOrBefore(i);
                if(null != startNodeQueue) {
                    startNodeQueueIter = startNodeQueue.iterator();
                    int prev_i = i;
                    while(startNodeQueueIter.hasNext()) {
                        startNode = startNodeQueueIter.next();
                        heap.add(new AlignHeapNode(null, 
                                    startNode,
                                    startNode.coverage,
                                    read.charAt(0),
                                    qualities.charAt(0),
                                    space));
                        i=startNode.position-1;
                        numStartNodesAdded++;
                    }
                }
                else {
                    i--;
                }
            }
        }
        else {
            alignmentStart = rec.getAlignmentStart();
            for(i=alignmentStart-offset;i<=alignmentStart+offset;) {
                startNodeQueue = graph.getPriorityQueueAtPositionOrGreater(i);
                if(null != startNodeQueue) {
                    startNodeQueueIter = startNodeQueue.iterator();
                    int prev_i = i;
                    while(startNodeQueueIter.hasNext()) {
                        startNode = startNodeQueueIter.next();
                        heap.add(new AlignHeapNode(null, 
                                    startNode,
                                    startNode.coverage,
                                    read.charAt(0),
                                    qualities.charAt(0),
                                    space));
                        i=startNode.position+1;
                        numStartNodesAdded++;
                    }
                }
                else {
                    i++;
                }
            }
        }
        if(numStartNodesAdded == 0) {
            throw new Exception("Did not add any start nodes!");
        }

        while(null != heap.peek()) {
            curAlignHeapNode = heap.poll();

            // HERE
                /*
                    //System.err.println("size:" + heap.size() + "\talignmentStart:" + alignmentStart + "\toffset:" + offset + "\treadOffset:" + curAlignHeapNode.readOffset);
                    System.err.print("size:" + heap.size() + ":" + curAlignHeapNode.readOffset + ":" + curAlignHeapNode.score + ":" + curAlignHeapNode.coverageSum + ":" + curAlignHeapNode.startPosition + ":");
                    curAlignHeapNode.node.print(System.err);
                */

            // Remove all non-insertions with the same contig/pos/read-offset/type/base and lower score 
            nextAlignHeapNode = heap.peek();
            while(Node.INSERTION != curAlignHeapNode.node.type 
                    && null != nextAlignHeapNode 
                    && 0 == comp.compare(curAlignHeapNode, nextAlignHeapNode)) 
            {
                if(curAlignHeapNode.score < nextAlignHeapNode.score ||
                        (curAlignHeapNode.score == nextAlignHeapNode.score && 
                         curAlignHeapNode.coverageSum < nextAlignHeapNode.coverageSum)) {
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
                /*
                    System.err.print(curAlignHeapNode.coverageSum + ":" + curAlignHeapNode.score + ":");
                    System.err.print(curAlignHeapNode.startPosition + ":");
                    curAlignHeapNode.node.print(System.err);
                    */

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

                if(strand) { // reverse
                    // Go to all the "prev" nodes
                    iter = curAlignHeapNode.node.prev.listIterator();
                    iterCov = curAlignHeapNode.node.prevCov.listIterator();
                }
                else { // forward
                    // Go to all "next" nodes
                    iter = curAlignHeapNode.node.next.listIterator();
                    iterCov = curAlignHeapNode.node.nextCov.listIterator();
                }
                while(iter.hasNext()) {
                    Node nextNode = iter.next();
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
        Align.updateSAM(rec, bestAlignHeapNode, space, read, strand);

        return rec;
    }

    private static void updateSAM(SAMRecord rec, AlignHeapNode bestAlignHeapNode, SRMAUtil.Space space, String read, boolean strand)
        throws Exception
    {
        AlignHeapNode curAlignHeapNode=null;
        AlignHeapNode prevAlignHeapNode=null;

        int alignmentStart = 0;
        int readIndex=-1;
        byte readBases[];
        byte colorErrors[];
        int i;
        String readColors=null, readColorQualities=null;

        if(null == bestAlignHeapNode) {
            System.err.println("\nNo alignments!");
            return;
        }

        // To generate a new CIGAR
        List<CigarElement> cigarElements=null;
        CigarOperator prevCigarOperator=null, curCigarOperator=null;
        int prevCigarOperatorLength=0;

        // TODO 
        // setInferredInsertSize (invalidates paired end reads)
        // setMappingQuality (?)
        // setFlag
        // update base qualities for color space reads 

        readBases = new byte[read.length()];
        if(strand) {
            readIndex=0;
        }
        else {
            readIndex = read.length()-1;
        }
        cigarElements = new LinkedList<CigarElement>();
        if(strand) {
            alignmentStart=bestAlignHeapNode.node.position;
        }
        else {
            alignmentStart=bestAlignHeapNode.startPosition;
        }

        assert null != bestAlignHeapNode;
        curAlignHeapNode = bestAlignHeapNode;
        while(null != curAlignHeapNode) {
            // Get the current cigar operator
            if(null != prevAlignHeapNode && CigarOperator.DELETION != prevCigarOperator && 1 < Math.abs(curAlignHeapNode.node.position - prevAlignHeapNode.node.position)) {
                //System.out.println("DEL");
                curCigarOperator = CigarOperator.DELETION;
            }
            else {
                switch(curAlignHeapNode.node.type) {
                    case Node.MATCH: // Fall through
                    case Node.MISMATCH:
                        curCigarOperator = CigarOperator.MATCH_OR_MISMATCH;
                        break;
                    case Node.INSERTION:
                        //System.out.println("INS");
                        curCigarOperator = CigarOperator.INSERTION;
                        break;
                    default:
                        throw new Exception("Unknown node type");
                }
                if(space == SRMAUtil.Space.COLORSPACE) {
                    readBases[readIndex]  = (byte)curAlignHeapNode.node.base;
                    if(strand) {
                        readIndex++;
                    }
                    else {
                        readIndex--;
                    }
                }
            }
            if(prevCigarOperator != curCigarOperator) {
                // different cigar operator

                // add the previous cigar operator
                if(null != prevCigarOperator) {
                    if(strand) { // reverse
                        // append 
                        cigarElements.add(new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                    }
                    else {
                        // prepend
                        cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                    }
                }

                // update prevCigarOperator
                prevCigarOperator = curCigarOperator;
                if(curCigarOperator == CigarOperator.DELETION) {
                    // length of deletion
                    prevCigarOperatorLength = Math.abs(curAlignHeapNode.node.position - prevAlignHeapNode.node.position) - 1;
                }
                else {
                    prevCigarOperatorLength=1;
                }
            }
            else {
                // same cigar operator
                prevCigarOperatorLength++;
            }

            // Update
            if(CigarOperator.DELETION != curCigarOperator) {
                prevAlignHeapNode = curAlignHeapNode;
                curAlignHeapNode = curAlignHeapNode.prev;
            }
        }
        if(0 < prevCigarOperatorLength) {
            if(null == prevCigarOperator || CigarOperator.DELETION == prevCigarOperator) { 
                throw new Exception("Ended with a null cigar operator or a deletion cigar operator");
            }
            if(strand) { // reverse
                // append 
                cigarElements.add(new CigarElement(prevCigarOperatorLength, prevCigarOperator));
            }
            else {
                // prepend
                cigarElements.add(0, new CigarElement(prevCigarOperatorLength, prevCigarOperator));
                // adjust if we start with an insertion
                if(CigarOperator.INSERTION == prevCigarOperator) {
                    alignmentStart++;
                }
            }
        }

        // Get color error string
        if(space == SRMAUtil.Space.COLORSPACE) {
            colorErrors = new byte[read.length()];
            char prevBase = SRMAUtil.COLORSPACE_ADAPTOR;
            for(i=0;i<read.length();i++) {
                if(SRMAUtil.colorSpaceNextBase(prevBase, read.charAt(i)) == readBases[i]) {
                    colorErrors[i] = Alignment.GAP;
                }
                else {
                    colorErrors[i] = (byte)read.charAt(i);
                }
            }
        }
        else {
            readBases = new byte[read.length()];
            if(strand) {
                for(i=0;i<read.length();i++) {
                    readBases[i] = (byte)read.charAt(read.length() - i - 1);
                }
            }
            else {
                for(i=0;i<read.length();i++) {
                    readBases[i] = (byte)read.charAt(i);
                }
            }
        }

        // Update SAM record
        rec.setCigar(new Cigar(cigarElements));
        //System.out.println("rec.getAlignmentStart()="+rec.getAlignmentStart());
        //System.out.println("alignmentStart="+alignmentStart);
        rec.setAlignmentStart(alignmentStart);
        if(strand) { // reverse
            // reverse read bases and qualities so it is on the + strand
            //SAMRecordUtil.reverseArray(readBases);
            // TODO: qualities
        }
        rec.setReadBases(readBases);
        if(null != rec.getAttribute("CS")) {
            readColors = (String)rec.getAttribute("CS");
            readColorQualities = (String)rec.getAttribute("CQ");
        }
        // Clear attributes
        rec.clearAttributes();
        // Set new attributes
        if(null != readColors) {
            rec.setAttribute("CS", readColors);
            rec.setAttribute("CQ", readColorQualities);
        }
        rec.setAttribute("AS", bestAlignHeapNode.score);
        rec.setAttribute("XC", bestAlignHeapNode.coverageSum);
        // set the XE attribute for colorError string
    }
}
