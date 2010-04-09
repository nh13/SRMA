/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import net.sf.samtools.*;
import net.sf.picard.reference.*;
import srma.*;

public class Align {

    private static final int MAX_HEAP_SIZE = 8192;

    public static SAMRecord align(Graph graph, SAMRecord rec, Node recNode, 
            ReferenceSequence sequence, 
            int offset, 
            double minimumAlleleFrequency,
            int minimumAlleleCoverage)
        throws Exception
    {

        int i;
        AlignHeapNode curAlignHeapNode = null;
        AlignHeapNode nextAlignHeapNode = null;
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
        //String readName = rec.getReadName();

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

        // Remove mate pair information
        Align.removeMateInfo(rec);

        // Bound by original alignment if possible
        bestAlignHeapNode = Align.boundWithOriginalAlignment(rec, 
                graph,
                recNode, 
                strand, 
                read, 
                qualities, 
                space, 
                sequence, 
                minimumAlleleFrequency,
                minimumAlleleCoverage);

        // HERE 
        /*
           if(null != bestAlignHeapNode) {
           System.err.println("FOUND BEST");
           }
           else {
           System.err.println("NOT FOUND (BEST)");
           }
        //Align.updateSAM(rec, bestAlignHeapNode, space, read, strand);
        //return rec;
        */

        if(strand) { // reverse
            comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MAXHEAP); 
            heap = new AlignHeap(AlignHeap.HeapType.MAXHEAP); 
        }
        else { // forward
            comp = new AlignHeapNodeComparator(AlignHeap.HeapType.MINHEAP); 
            heap = new AlignHeap(AlignHeap.HeapType.MINHEAP); 
        }

        // Add start nodes
        if(strand) { // reverse
            alignmentStart = rec.getAlignmentEnd();
            for(i=alignmentStart+offset;alignmentStart-offset<=i;i--) {
                int position = graph.getPriorityQueueIndexAtPositionOrBefore(i);
                PriorityQueue<Node> startNodeQueue = graph.getPriorityQueue(position);
                if(0 != position && null != startNodeQueue) {
                    Iterator<Node> startNodeQueueIter = startNodeQueue.iterator();
                    while(startNodeQueueIter.hasNext()) {
                        Node startNode = startNodeQueueIter.next();
                        if(passFilters(graph,
                                    startNode,
                                    minimumAlleleFrequency,
                                    minimumAlleleCoverage)) {
                            heap.add(new AlignHeapNode(null, 
                                        startNode,
                                        startNode.coverage,
                                        read.charAt(0),
                                        qualities.charAt(0),
                                        space));
                        }
                        if(startNode.position < i) {
                            i = startNode.position;
                        }
                        numStartNodesAdded++;
                    }
                }
            }
        }
        else {
            alignmentStart = rec.getAlignmentStart();
            for(i=alignmentStart-offset;i<=alignmentStart+offset;i++) {
                int position = graph.getPriorityQueueIndexAtPositionOrGreater(i);
                PriorityQueue<Node> startNodeQueue = graph.getPriorityQueue(position);
                if(0 != position && null != startNodeQueue) {
                    Iterator<Node> startNodeQueueIter = startNodeQueue.iterator();
                    while(startNodeQueueIter.hasNext()) {
                        Node startNode = startNodeQueueIter.next();
                        if(passFilters(graph,
                                    startNode,
                                    minimumAlleleFrequency,
                                    minimumAlleleCoverage)) {
                            heap.add(new AlignHeapNode(null, 
                                        startNode,
                                        startNode.coverage,
                                        read.charAt(0),
                                        qualities.charAt(0),
                                        space));
                        }
                        if(i < startNode.position) {
                            i = startNode.position;
                        }
                        numStartNodesAdded++;
                    }
                }
            }
        }
        if(numStartNodesAdded == 0) {
            throw new Exception("Did not add any start nodes!");
        }

        // Get first node off the heap
        curAlignHeapNode = heap.poll();

        while(null != curAlignHeapNode) {

            if(Align.MAX_HEAP_SIZE <= heap.size()) {
                // too many to consider
                return rec;
            }

            // HERE
            //System.err.println("strand:" + strand + "\tsize:" + heap.size() + "\talignmentStart:" + alignmentStart + "\toffset:" + offset + "\treadOffset:" + curAlignHeapNode.readOffset);
            //System.err.print("size:" + heap.size() + ":" + curAlignHeapNode.readOffset + ":" + curAlignHeapNode.score + ":" + curAlignHeapNode.alleleCoverageSum + ":" + curAlignHeapNode.startPosition + "\t");
            //curAlignHeapNode.node.print(System.err);
            //System.err.print("\rposition:" + curAlignHeapNode.node.position + "\treadOffset:" + curAlignHeapNode.readOffset);

            // Remove all non-insertions with the same contig/pos/read-offset/type/base and lower score 
            nextAlignHeapNode = heap.peek();
            while(Node.INSERTION != curAlignHeapNode.node.type 
                    && null != nextAlignHeapNode 
                    && 0 == comp.compare(curAlignHeapNode, nextAlignHeapNode)) 
            {
                if(curAlignHeapNode.score < nextAlignHeapNode.score ||
                        (curAlignHeapNode.score == nextAlignHeapNode.score && 
                         curAlignHeapNode.alleleCoverageSum < nextAlignHeapNode.alleleCoverageSum)) {
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

                // HERE
                //System.err.print(curAlignHeapNode.alleleCoverageSum + ":" + curAlignHeapNode.score + ":");
                //System.err.print(curAlignHeapNode.startPosition + ":");
                //curAlignHeapNode.node.print(System.err);

                if(null == bestAlignHeapNode 
                        || bestAlignHeapNode.score < curAlignHeapNode.score 
                        || (bestAlignHeapNode.score == curAlignHeapNode.score 
                            && bestAlignHeapNode.alleleCoverageSum < curAlignHeapNode.alleleCoverageSum)) 
                {
                    bestAlignHeapNode = curAlignHeapNode;
                }
            }
            else {
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
                    if(passFilters(graph,
                                nextNode,
                                nextCoverage,
                                minimumAlleleFrequency,
                                minimumAlleleCoverage)) {
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
            // Get next node
            curAlignHeapNode = heap.poll();
        }

        // Recover alignment
        Align.updateSAM(rec, bestAlignHeapNode, space, read, strand);
        
        return rec;
    }

    private static void removeMateInfo(SAMRecord rec)
    {
        if(rec.getReadPairedFlag() || rec.getFirstOfPairFlag() || rec.getSecondOfPairFlag()) { // paired end
            // Remove all information of its mate
            
            // flag
            rec.setProperPairFlag(false); // not paired any more
            rec.setMateUnmappedFlag(false);
            rec.setMateNegativeStrandFlag(false);
            rec.setFirstOfPairFlag(false);
            rec.setFirstOfSecondFlag(false);

            // entries
            rec.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec.setMateAlignmentStart(0);
            rec.setInferredInsertSize(0);

             // TODO: remove tags and values that are mate pair inclined.
        }
    }

    private static boolean passFilters(Graph graph,
            Node node,
            int nodeCoverage,
            double minimumAlleleFrequency,
            int minimumAlleleCoverage) 
    {
        // must pass either minimum allele frequency or minimum allele coverage
        // - the former does not penalize low coverage regions, while the later does 
        int coverage = graph.getCoverage(node.position);
        if(minimumAlleleFrequency <= nodeCoverage / ((double)coverage) || minimumAlleleCoverage <= nodeCoverage) {
            return true;
        }
        else {
            return false;
        }
    }

    private static boolean passFilters(Graph graph,
            Node node,
            double minimumAlleleFrequency,
            int minimumAlleleCoverage) 
    {
        return passFilters(graph, node, node.coverage, minimumAlleleFrequency, minimumAlleleCoverage);
    }


    private static AlignHeapNode boundWithOriginalAlignment(SAMRecord rec, 
            Graph graph,
            Node recNode, 
            boolean strand, 
            String read, 
            String qualities, 
            SRMAUtil.Space space,
            ReferenceSequence sequence, 
            double minimumAlleleFrequency,
            int minimumAlleleCoverage)
        throws Exception
    {
        Alignment alignment = null;
        AlignHeapNode curAlignHeapNode = null;
        ListIterator<Node> iter=null;
        ListIterator<Integer> iterCov=null;
        AlignHeap heap = null;

        // Cannot bound
        if(!passFilters(graph,
                    recNode,
                    minimumAlleleFrequency,
                    minimumAlleleCoverage)) {
            return null;
        }

        // Get original alignment
        alignment = new Alignment(rec, sequence);

        // Initialize heap
        if(strand) { // reverse
            heap = new AlignHeap(AlignHeap.HeapType.MAXHEAP); 
        }
        else { // forward
            heap = new AlignHeap(AlignHeap.HeapType.MINHEAP); 
        }
        heap.add(new AlignHeapNode(null, 
                    recNode,
                    recNode.coverage,
                    read.charAt(0),
                    qualities.charAt(0),
                    space));

        // HERE
        /*
           System.err.println("\n" + rec.getAlignmentStart() + ":" + rec.getAlignmentEnd());
           System.err.println("reverse:" + strand);
           recNode.print(System.err);
           alignment.print(System.err); 
           */

        curAlignHeapNode = heap.poll();
        while(null != curAlignHeapNode) {

            // Check if alignment was found
            if(curAlignHeapNode.readOffset == read.length() - 1) { // Found
                return curAlignHeapNode;
            }

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

            // We should find original alignment
            if(iter.hasNext()) {

                // Get the expected next position in the alignment
                int nextReadOffset = (strand) ? (read.length() - 1 - curAlignHeapNode.readOffset - 1) : (curAlignHeapNode.readOffset + 1);
                int nextPosition = alignment.positions[nextReadOffset];
                int nextIndex = alignment.positionsIndex[nextReadOffset];
                char nextBase = (char)alignment.read[nextIndex];
                int nextType = -1;

                /*
                   System.err.println("HERE BOUND: " + nextPosition + ":" + nextIndex + ":" + (char)alignment.reference[nextIndex] + ":" + (char)alignment.read[nextIndex]);
                   */

                if(Alignment.GAP == alignment.reference[nextIndex]) {
                    nextType = Node.INSERTION;
                }
                else if(alignment.read[nextIndex] == alignment.reference[nextIndex]) {
                    nextType = Node.MATCH;
                }
                else {
                    // HERE
                    if(Alignment.GAP == alignment.read[nextIndex]) {
                        throw new Exception("Alignment error");
                    }
                    nextType = Node.MISMATCH;
                }

                // HERE
                if(nextBase == Alignment.GAP) {
                    throw new Error("Alignment error");
                }

                // HERE
                /*
                   System.err.println("nextReadOffset: " + nextReadOffset
                   + "\tnextPosition: " + nextPosition
                   + "\tnextBase: " + nextBase
                   + "\tnextType: " + nextType);
                   */

                while(iter.hasNext()) {
                    Node nextNode = iter.next();
                    int nextCoverage = iterCov.next();

                    // HERE
                    /*
                       System.err.println("Found nextNode.position: " + nextNode.position
                       + "\tnextNode.type: " + nextNode.type
                       + "\tnextNode.offset: " + nextNode.offset
                       + "\tnextNode.base: " + nextNode.base
                       + "\tnextNode.position: " + nextNode.position);
                       */

                    if(nextNode.position == nextPosition && nextNode.type == nextType && nextNode.base == nextBase) { // bases match
                        if(passFilters(graph,
                                    nextNode,
                                    nextCoverage,
                                    minimumAlleleFrequency,
                                    minimumAlleleCoverage)) {
                            heap.add(new AlignHeapNode(curAlignHeapNode, 
                                        nextNode, 
                                        nextCoverage,
                                        read.charAt(curAlignHeapNode.readOffset+1), 
                                        qualities.charAt(curAlignHeapNode.readOffset+1), 
                                        space));
                        }
                        else {
                            return null;
                        }
                    }
                }
            }
            iter=null;
            iterCov=null;

            // Get next
            curAlignHeapNode = heap.poll();
        }

        System.err.println("");
        alignment.print(System.err);
        throw new Exception("Control reached unexpected point");
    }

    private static void updateSAM(SAMRecord rec, AlignHeapNode bestAlignHeapNode, SRMAUtil.Space space, String read, boolean strand)
        throws Exception
    {
        AlignHeapNode curAlignHeapNode=null;
        AlignHeapNode prevAlignHeapNode=null;

        int alignmentStart = 0;
        int readIndex=-1;
        byte readBases[] = null;
        byte colorErrors[] = null;
        int i;
        String readColors=null, readColorQualities=null;

        // Debugging stuff
        //String readName = rec.getReadName();

        if(null == bestAlignHeapNode) {
            // Do not modify the alignment
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

        // Adjust position on the (-) strand if we end with an insertion
        if(strand && Node.INSERTION == bestAlignHeapNode.node.type) { // reverse
            alignmentStart++;
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
        rec.setAttribute("XC", bestAlignHeapNode.alleleCoverageSum);
        // set the XE attribute for colorError string
        //rec.setAttribute("CE", colorErrors);
    }
}
