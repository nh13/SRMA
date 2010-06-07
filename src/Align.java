/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import net.sf.samtools.*;
import net.sf.picard.reference.*;
import srma.*;

public class Align {

    private static final int CORRECT_BASE_QUALITY_PENALTY = 20; // TODO: should be a parameter to SRMA

    public static SAMRecord align(Graph graph, SAMRecord rec, Node recNode, 
            ReferenceSequence sequence, 
            SAMProgramRecord programRecord,
            int offset, 
            AlleleCoverageCutoffs alleleCoverageCutoffs,
            boolean correctBases,
            boolean useSequenceQualities,
            int MAXIMUM_TOTAL_COVERAGE,
            int MAX_HEAP_SIZE)
        throws Exception
    {

        int i;
        AlignHeapNode curAlignHeapNode = null;
        AlignHeapNode nextAlignHeapNode = null;
        AlignHeapNode bestAlignHeapNode=null;
        AlignHeap heap=null;
        String read=null; // could be cs
        String readBases = null; // always nt
        String qualities=null; // could be cq
        SRMAUtil.Space space=SRMAUtil.Space.NTSPACE;
        ListIterator<Node> iter=null;
        ListIterator<Integer> iterCov=null;
        AlignHeapNodeComparator comp=null;
        int alignmentStart = -1;
        int numStartNodesAdded = 0;
        boolean strand = rec.getReadNegativeStrandFlag(); // false -> forward, true -> reverse

        // Debugging stuff
        //String readName = rec.getReadName();

        assert SRMAUtil.Space.COLORSPACE != space;

        // Get space
        read = (String)rec.getAttribute("CS");
        if(null == read) {
            // Use base space
            space = SRMAUtil.Space.NTSPACE;
        }
        else {
            // assumes CS and CQ are always in sequencing order
            space = SRMAUtil.Space.COLORSPACE;
        }

        // Get read and qualities
        if(space == SRMAUtil.Space.NTSPACE) {
            byte tmpRead[] = rec.getReadString().getBytes();
            byte tmpQualities[] = rec.getBaseQualityString().getBytes();
            // Reverse once
            if(strand) { // reverse
                SAMRecordUtil.reverseArray(tmpRead);
                SAMRecordUtil.reverseArray(tmpQualities);
            }
            read = new String(tmpRead);
            readBases = new String(tmpRead);
            qualities = new String(tmpQualities);
            if(strand) { // reverse
                SAMRecordUtil.reverseArray(tmpRead);
                SAMRecordUtil.reverseArray(tmpQualities);
            }
        }
        else {
            byte tmpRead[] = rec.getReadString().getBytes();
            // Reverse once
            if(strand) { // reverse
                SAMRecordUtil.reverseArray(tmpRead);
            }
            readBases = new String(tmpRead);
            if(strand) { // reverse
                SAMRecordUtil.reverseArray(tmpRead);
            }
            read = SRMAUtil.normalizeColorSpaceRead(read);
            qualities = (String)rec.getAttribute("CQ");
            // Some aligners include a quality value for the adapter.  A quality value
            // IMHO should not be given for an unobserved (assumed) peice of data.  Trim
            // the first quality in this case
            if(qualities.length() == 1 + read.length()) { // trim the f
                qualities = qualities.substring(1);
            }
        }
        // Reverse back
        if(readBases.length() <= 0) {
            throw new Exception("Error.  The current alignment has no bases.");
        }
        if(read.length() <= 0) {
            throw new Exception("Error.  The current alignment has no bases.");
        }
        if(qualities.length() <= 0) {
            throw new Exception("Error.  The current alignment has no qualities.");
        }
        if(readBases.length() != read.length()) {
            if(space == SRMAUtil.Space.COLORSPACE) {
                throw new Exception("Error.  The current alignment's read bases length does not match the length of the colors in the CS tag.");
            }
            else {
                throw new Exception("Error.  Internal error: readBases.length() != read.length()");
            }
        }

        // Remove mate pair information
        Align.removeMateInfo(rec);

        comp = new AlignHeapNodeComparator((strand) ? AlignHeap.HeapType.MAXHEAP : AlignHeap.HeapType.MINHEAP);

        // Bound by original alignment if possible
        bestAlignHeapNode = Align.boundWithOriginalAlignment(rec, 
                graph,
                recNode, 
                comp,
                strand, 
                read,
                qualities,
                readBases,
                space, 
                sequence, 
                alleleCoverageCutoffs,
                useSequenceQualities,
                MAXIMUM_TOTAL_COVERAGE,
                MAX_HEAP_SIZE);

        // HERE 
        /*
        //System.err.println("readName="+rec.getReadName());
        if(null != bestAlignHeapNode) {
        System.err.println("\nFOUND BEST:" + rec.toString());
        }
        else {
        System.err.println("\nNOT FOUND (BEST)" + rec.toString());
        }
        //Align.updateSAM(rec, programRecord, bestAlignHeapNode, space, read, qualities, strand, correctBases);
        //return rec;
        */

        heap = new AlignHeap((strand) ? AlignHeap.HeapType.MAXHEAP : AlignHeap.HeapType.MINHEAP);

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
                        int f = passFilters(graph,
                                startNode,
                                alleleCoverageCutoffs,
                                MAXIMUM_TOTAL_COVERAGE);
                        if(0 == f) {
                            heap.add(new AlignHeapNode(null, 
                                        startNode,
                                        startNode.coverage,
                                        read.charAt(0),
                                        qualities.charAt(0),
                                        useSequenceQualities,
                                        space));
                        }
                        else if(f < 0) {
                            return rec;
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
                        int f = passFilters(graph,
                                startNode,
                                alleleCoverageCutoffs,
                                MAXIMUM_TOTAL_COVERAGE);
                        if(0 == f) {
                            heap.add(new AlignHeapNode(null, 
                                        startNode,
                                        startNode.coverage,
                                        read.charAt(0),
                                        qualities.charAt(0),
                                        useSequenceQualities,
                                        space));
                        }
                        else if(f < 0) {
                            return rec;
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

            if(MAX_HEAP_SIZE <= heap.size()) {
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
                // All read bases examined, store if has the best alignment.

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
            else if(null != bestAlignHeapNode && curAlignHeapNode.score < bestAlignHeapNode.score) {
                // ignore, under the assumption that scores can only become more negative.
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
                    int f = passFilters(graph,
                            nextNode,
                            nextCoverage,
                            alleleCoverageCutoffs,
                            MAXIMUM_TOTAL_COVERAGE);
                    if(0 == f) {
                        heap.add(new AlignHeapNode(curAlignHeapNode, 
                                    nextNode, 
                                    nextCoverage,
                                    read.charAt(curAlignHeapNode.readOffset+1), 
                                    qualities.charAt(curAlignHeapNode.readOffset+1), 
                                    useSequenceQualities,
                                    space));
                    }
                    else if(f < 0) {
                        return rec;
                    }
                }
                iter=null;
            }
            // Get next node
            curAlignHeapNode = heap.poll();
        }

        // Recover alignment
        Align.updateSAM(rec, programRecord, bestAlignHeapNode, space, read, qualities, strand, correctBases);

        return rec;
    }

    private static void removeMateInfo(SAMRecord rec)
    {
        if(rec.getReadPairedFlag()) {
            // Remove all information of its mate

            // flag
            rec.setProperPairFlag(false); // not paired any more
            rec.setMateUnmappedFlag(false);
            rec.setMateNegativeStrandFlag(false);
            //rec.setFirstOfPairFlag(false);
            //rec.setSecondOfPairFlag(false);

            // entries
            rec.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
            rec.setMateAlignmentStart(0);
            rec.setInferredInsertSize(0);

            // TODO: remove tags and values that are mate pair inclined.
        }
    }

    private static AlignHeapNode boundWithOriginalAlignment(SAMRecord rec, 
            Graph graph,
            Node recNode, 
            AlignHeapNodeComparator comp,
            boolean strand, 
            String read, // could be cs 
            String qualities, // could be cq
            String readBases, // always nt
            SRMAUtil.Space space,
            ReferenceSequence sequence, 
            AlleleCoverageCutoffs alleleCoverageCutoffs,
            boolean useSequenceQualities,
            int MAXIMUM_TOTAL_COVERAGE,
            int MAX_HEAP_SIZE) 
        throws Exception
    {
        AlignHeapNode curAlignHeapNode = null;
        AlignHeapNode nextAlignHeapNode = null;
        AlignHeapNode bestAlignHeapNode=null;
        ListIterator<Node> iter=null;
        ListIterator<Integer> iterCov=null;
        AlignHeap heap = null;

        // Cannot bound
        if(0 != passFilters(graph,
                    recNode,
                    alleleCoverageCutoffs,
                    MAXIMUM_TOTAL_COVERAGE)) {
            return null;
        }

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
                    useSequenceQualities,
                    space));

        curAlignHeapNode = heap.poll();
        while(null != curAlignHeapNode) {
            if(MAX_HEAP_SIZE <= heap.size()) {
                // too many to consider
                return null;
            }
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

            if(curAlignHeapNode.readOffset == readBases.length() - 1) { // found, keep beset
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


                // Get the expected next position in the alignment
                while(iter.hasNext()) {
                    Node nextNode = iter.next();
                    int nextCoverage = iterCov.next();

                    // Base should match alignment
                    if(nextNode.base == readBases.charAt(curAlignHeapNode.readOffset+1)) {
                        int f = passFilters(graph, 
                                nextNode, 
                                nextCoverage, 
                                alleleCoverageCutoffs,
                                MAXIMUM_TOTAL_COVERAGE);
                        if(0 == f) {
                            heap.add(new AlignHeapNode(curAlignHeapNode, 
                                        nextNode, 
                                        nextCoverage,
                                        read.charAt(curAlignHeapNode.readOffset+1), 
                                        qualities.charAt(curAlignHeapNode.readOffset+1), 
                                        useSequenceQualities,
                                        space));
                        }
                        else if(f < 0) {
                            return null;
                        }
                    }
                }
                iter=null;
                iterCov=null;
            }

            // Get next
            curAlignHeapNode = heap.poll();
        }

        return bestAlignHeapNode;
    }

    private static void updateSAM(SAMRecord rec, SAMProgramRecord programRecord, AlignHeapNode bestAlignHeapNode, SRMAUtil.Space space, String read, String qualities, boolean strand, boolean correctBases)
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

        // Get color space attributes
        if(null != rec.getAttribute("CS")) {
            readColors = (String)rec.getAttribute("CS");
            readColorQualities = (String)rec.getAttribute("CQ");
        }
        // Clear attributes
        rec.clearAttributes();

        readBases = new byte[read.length()];
        if(strand) {
            readIndex=0;
        }
        else {
            readIndex = read.length()-1;
        }
        cigarElements = new LinkedList<CigarElement>();
        if(strand) { // reverse strand is the current position
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
                if(space == SRMAUtil.Space.COLORSPACE || correctBases) {
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
            }
        }

        if(space == SRMAUtil.Space.COLORSPACE) { // color space, read bases already inferred
            // Get color error string
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
        else if(correctBases) { // bases were corrected
            byte newBaseQualities[] = new byte[read.length()];
            if(strand) {
                for(i=0;i<read.length();i++) {
                    if(readBases[i] == (byte)read.charAt(read.length() - i - 1)) {
                        newBaseQualities[i] = (byte)(qualities.charAt(read.length() - i - 1) - 33);
                    }
                    else {
                        // TODO: how much to down-weight ?
                        newBaseQualities[i] = (byte)(SRMAUtil.QUAL2CHAR(SRMAUtil.CHAR2QUAL(qualities.charAt(read.length() - i - 1)) - CORRECT_BASE_QUALITY_PENALTY) - 33);
                        if(newBaseQualities[i] <= 0) {
                            newBaseQualities[i]=1;
                        }
                    }
                }
            }
            else {
                for(i=0;i<read.length();i++) {
                    if(readBases[i] == (byte)read.charAt(i)) {
                        newBaseQualities[i] = (byte)(qualities.charAt(i) - 33);
                    }
                    else {
                        // TODO: how much to down-weight ?
                        newBaseQualities[i] = (byte)(SRMAUtil.QUAL2CHAR(SRMAUtil.CHAR2QUAL(qualities.charAt(i)) - CORRECT_BASE_QUALITY_PENALTY) - 33);
                        if(newBaseQualities[i] <= 0) {
                            newBaseQualities[i]=1;
                        }
                    }
                }
            }
            rec.setBaseQualities(newBaseQualities);
            rec.setAttribute("XO", read);
            rec.setAttribute("XQ", qualities);
        }
        else { // bases not corrected 
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
        rec.setAlignmentStart(alignmentStart);
        rec.setReadBases(readBases);
        // Set new attributes
        if(null != readColors) {
            rec.setAttribute("CS", readColors);
            rec.setAttribute("CQ", readColorQualities);
        }
        rec.setAttribute("AS", bestAlignHeapNode.score);
        rec.setAttribute("XC", bestAlignHeapNode.alleleCoverageSum);
        rec.setAttribute("PG", programRecord.getId());
        // set the XE attribute for colorError string
        //rec.setAttribute("CE", colorErrors);
    }

    /*
     * -1 if the alignment process should be aborted 
     *  0 if the alignment should continue 
     *  1 if the alignment should not be considered any further
     * */
    private static int passFilters(Graph graph,
            Node node,
            int toNodeCoverage,
            AlleleCoverageCutoffs alleleCoverageCutoffs,
            int MAXIMUM_TOTAL_COVERAGE) 
    {
        int totalCoverage = graph.getCoverage(node.position);
        if(MAXIMUM_TOTAL_COVERAGE < totalCoverage) {
            return -1;
        }
        else if(alleleCoverageCutoffs.getQ(totalCoverage) <= toNodeCoverage) {
            // HERE
            //System.err.println("TRUE totalCoverage="+totalCoverage+"\ttoNodeCoverage="+toNodeCoverage+"\tcutoff="+alleleCoverageCutoffs.getQ(totalCoverage));
            return 0;
        }
        else {
            // HERE
            //System.err.println("FALSE totalCoverage="+totalCoverage+"\ttoNodeCoverage="+toNodeCoverage+"\tcutoff="+alleleCoverageCutoffs.getQ(totalCoverage));
            return 1;
        }
    }

    private static int passFilters(Graph graph,
            Node node,
            AlleleCoverageCutoffs alleleCoverageCutoffs,
            int MAXIMUM_TOTAL_COVERAGE) 
    {
        return passFilters(graph, node, node.coverage, alleleCoverageCutoffs, MAXIMUM_TOTAL_COVERAGE);
    }
}
