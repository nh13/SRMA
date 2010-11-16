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

    private static final List<String> saveTags =
        Arrays.asList("RG", "LB", "PU", "PG", "CS", "CQ");

    public static void align(Graph graph, SAMRecord rec, Node recNode, 
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
        ListIterator<Node.NodeRecord> iter=null;
        AlignHeapNodeComparator comp=null;
        int alignmentStart = -1;
        int numStartNodesAdded = 0;

        AlignData data = null;

        boolean strand = rec.getReadNegativeStrandFlag(); // false -> forward, true -> reverse
        data = new AlignData(rec);

        space = data.space;
        readBases = data.readBases;
        if(space == SRMAUtil.Space.NTSPACE) {
            read = data.readBases;
            qualities = data.readBaseQualities;
        }
        else {
            read = data.readColors;
            qualities = data.readColorQualities;
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
                            return;
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
                            return;
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
                return;
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

            // Check if the alignment is complete
            if(curAlignHeapNode.readOffset == read.length() - 1) {
                // All read bases examined, store if has the best alignment.

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
                }
                else { // forward
                    // Go to all "next" nodes
                    iter = curAlignHeapNode.node.next.listIterator();
                }
                while(iter.hasNext()) {
                    Node.NodeRecord next = iter.next();
                    int f = passFilters(graph,
                            next.node,
                            next.coverage,
                            alleleCoverageCutoffs,
                            MAXIMUM_TOTAL_COVERAGE);
                    if(0 == f) {
                        heap.add(new AlignHeapNode(curAlignHeapNode, 
                                    next.node,
                                    next.coverage,
                                    read.charAt(curAlignHeapNode.readOffset+1), 
                                    qualities.charAt(curAlignHeapNode.readOffset+1), 
                                    useSequenceQualities,
                                    space));
                    }
                    else if(f < 0) {
                        return;
                    }
                }
                iter=null;
            }
            // Get next node
            curAlignHeapNode = heap.poll();
        }

        // Recover alignment
        Align.updateSAM(rec, programRecord, bestAlignHeapNode, read, qualities, data, correctBases);
    }

    private static void removeMateInfo(SAMRecord rec)
    {
        if(rec.getReadPairedFlag()) {
            // Remove all information of its mate

            // flag
            rec.setProperPairFlag(false); // not paired any more
            rec.setMateUnmappedFlag(false);
            rec.setMateNegativeStrandFlag(false);

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
        AlignHeapNode bestAlignHeapNode = null;
        ListIterator<Node.NodeRecord> iter=null;
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

        // Add start nodes
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
                         curAlignHeapNode.alleleCoverageSum < nextAlignHeapNode.alleleCoverageSum)) 
                {
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
                }
                else { // forward
                    // Go to all "next" nodes
                    iter = curAlignHeapNode.node.next.listIterator();
                }

                // Get the expected next position in the alignment
                while(iter.hasNext()) {
                    Node.NodeRecord next = iter.next();

                    // Base should match alignment
                    if(next.node.base == readBases.charAt(curAlignHeapNode.readOffset+1)) {
                        int f = passFilters(graph, 
                                next.node,
                                next.coverage,
                                alleleCoverageCutoffs,
                                MAXIMUM_TOTAL_COVERAGE);
                        if(0 == f) {
                            heap.add(new AlignHeapNode(curAlignHeapNode, 
                                        next.node,
                                        next.coverage,
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
            }

            // Get next
            curAlignHeapNode = heap.poll();
        }

        return bestAlignHeapNode;
    }

    private static byte getColorQuality(byte e1, byte e2, byte q1, byte q2)
    {
        int val;
        if(e1 == (byte)Alignment.GAP && e2 == (byte)Alignment.GAP) {
            val = q1 + q2 + 10;
        }
        else if(e1 == (byte)Alignment.GAP) {
            val = q1 - q2;
        }
        else if(e2 == (byte)Alignment.GAP) {
            val = q2 - q1;
        }
        else {
            val = 1;
        }
        if(val <= 0) {
            val = 1;
        }
        else if(63 < val) {
            val = 63;
        }
        return (byte)val;
    }

    private static void updateSAM(SAMRecord rec, SAMProgramRecord programRecord, AlignHeapNode bestAlignHeapNode, 
            String read, String qualities, AlignData data, boolean correctBases)
        throws Exception
    {
        AlignHeapNode curAlignHeapNode=null;
        AlignHeapNode prevAlignHeapNode=null;

        int alignmentStart = 0;
        int readIndex=-1;
        byte readBases[] = null;
        byte baseQualities[] = null;
        byte colorErrors[] = null;
        int i;
        List<String> optFieldTags = new LinkedList<String>();
        List<Object> optFieldValues = new LinkedList<Object>();
        Object attr;
        SRMAUtil.Space space;
        boolean strand;
        List<CigarElement> tmpCigarElements=null;

        if(null == bestAlignHeapNode) {
            // Do not modify the alignment
            return;
        }

        strand = data.strand;
        space = data.space;

        // To generate a new CIGAR
        List<CigarElement> cigarElements=null;
        CigarOperator prevCigarOperator=null, curCigarOperator=null;
        int prevCigarOperatorLength=0;

        // TODO 
        // setInferredInsertSize (invalidates paired end reads)
        // setMappingQuality (?)
        // setFlag
        // update base qualities for color space reads 

        // clear attributes, but save some
        Align.clearAttributes(rec, optFieldTags, optFieldValues);

        readBases = new byte[read.length()];
        baseQualities = new byte[qualities.length()];
        for(i=0;i<qualities.length();i++) {
            // Must subtract 33 for PHRED scaling
            baseQualities[i] = (byte)(qualities.charAt(i) - 33);
        }

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
            if(strand) { // reverse
                for(i=0;i<read.length();i++) {
                    char nextBase = SRMAUtil.colorSpaceNextBase(prevBase, read.charAt(i));
                    if(nextBase == SRMAUtil.getCompliment((char)readBases[read.length()-i-1])) {
                        colorErrors[i] = (byte)Alignment.GAP;
                    }
                    else {
                        colorErrors[i] = (byte)read.charAt(i);
                    }
                    if(0 < i) {
                        // qualities are assumed to be always in the same direction as the color errors
                        baseQualities[read.length()-i] = getColorQuality(colorErrors[i-1],
                                colorErrors[i],
                                (byte)(qualities.charAt(i-1) - 33),
                                (byte)(qualities.charAt(i) - 33));
                    }
                    prevBase = SRMAUtil.getCompliment((char)readBases[read.length()-i-1]);
                }
                // last color
                baseQualities[0] = (byte)(qualities.charAt(read.length()-1)-33);
            }
            else {
                for(i=0;i<read.length();i++) {
                    char nextBase = SRMAUtil.colorSpaceNextBase(prevBase, read.charAt(i));
                    if(nextBase == readBases[i]) {
                        colorErrors[i] = (byte)Alignment.GAP;
                    }
                    else {
                        colorErrors[i] = (byte)read.charAt(i);
                    }
                    if(0 < i) {
                        baseQualities[i-1] = getColorQuality(colorErrors[i-1],
                                colorErrors[i],
                                (byte)(qualities.charAt(i-1) - 33),
                                (byte)(qualities.charAt(i) - 33));
                    }
                    prevBase = (char)readBases[i];
                }
                // last color
                baseQualities[read.length()-1] = (byte)(qualities.charAt(read.length()-1)-33);
            }
        }
        else if(correctBases) { // bases were corrected
            if(strand) {
                for(i=0;i<read.length();i++) {
                    if(readBases[i] == (byte)read.charAt(read.length() - i - 1)) {
                        baseQualities[i] = (byte)(qualities.charAt(read.length() - i - 1) - 33);
                    }
                    else {
                        // TODO: how much to down-weight ?
                        baseQualities[i] = (byte)(SRMAUtil.QUAL2CHAR(SRMAUtil.CHAR2QUAL(qualities.charAt(read.length() - i - 1)) - CORRECT_BASE_QUALITY_PENALTY) - 33);
                        if(baseQualities[i] <= 0) {
                            baseQualities[i]=1;
                        }
                    }
                }
            }
            else {
                for(i=0;i<read.length();i++) {
                    if(readBases[i] == (byte)read.charAt(i)) {
                        baseQualities[i] = (byte)(qualities.charAt(i) - 33);
                    }
                    else {
                        // TODO: how much to down-weight ?
                        baseQualities[i] = (byte)(SRMAUtil.QUAL2CHAR(SRMAUtil.CHAR2QUAL(qualities.charAt(i)) - CORRECT_BASE_QUALITY_PENALTY) - 33);
                        if(baseQualities[i] <= 0) {
                            baseQualities[i]=1;
                        }
                    }
                }
            }
            rec.setAttribute("XO", read);
            rec.setAttribute("XQ", qualities);
        }
        else { // bases not corrected 
            readBases = new byte[read.length()];
            baseQualities = new byte[qualities.length()]; // qualities.length() == read.length()
            if(strand) { // reverse
                for(i=0;i<read.length();i++) {
                    readBases[i] = (byte)read.charAt(read.length() - i - 1);
                    baseQualities[i] = (byte)(qualities.charAt(read.length() - i -1) - 33);
                }
            }
            else {
                for(i=0;i<read.length();i++) {
                    readBases[i] = (byte)read.charAt(i);
                    baseQualities[i] = (byte)(qualities.charAt(i) - 33);
                }
            }
        }

        // Add in soft-clipping and hard-clipping
        tmpCigarElements = rec.getCigar().getCigarElements();
        if(data.clippedStart) {
            CigarElement e1 = tmpCigarElements.get(0); // first
            int l = e1.getLength();
            CigarOperator op = e1.getOperator();

            // add cigar operator
            cigarElements.add(0, new CigarElement(l, op));

            if(CigarOperator.S == op) { // soft-clipping, add bases & quality
                byte tmpBases[] = new byte[readBases.length + l];
                // copy over aligned bases
                System.arraycopy(readBases, 0, tmpBases, l, readBases.length);
                readBases = tmpBases;
                // copy over soft-clipped bases
                System.arraycopy(rec.getReadBases(), 0, readBases, 0, l);

                byte tmpQualities[] = new byte[baseQualities.length + l];
                // copy over aligned qualities
                System.arraycopy(baseQualities, 0, tmpQualities, l, baseQualities.length);
                baseQualities = tmpQualities;
                // copy over soft-clipped qualities
                System.arraycopy(rec.getBaseQualities(), 0, baseQualities, 0, l);
            }
        }
        if(data.clippedEnd) {
            CigarElement e2 = tmpCigarElements.get(tmpCigarElements.size()-1); // last 
            int l = e2.getLength();
            CigarOperator op = e2.getOperator();

            // add cigar operator
            cigarElements.add(new CigarElement(l, op));

            if(CigarOperator.S == op) { // soft-clipping, add bases & quality
                byte tmpBases[] = new byte[readBases.length + l];
                // copy over aligned bases 
                System.arraycopy(readBases, 0, tmpBases, 0, readBases.length);
                // copy over soft-clipped bases
                System.arraycopy(rec.getReadBases(), readBases.length, tmpBases, readBases.length, l);
                readBases = tmpBases;

                byte tmpQualities[] = new byte[baseQualities.length + l];
                // copy over aligned qualities
                System.arraycopy(baseQualities, 0, tmpQualities, 0, baseQualities.length);
                // copy over soft-clipped qualities
                System.arraycopy(rec.getBaseQualities(), baseQualities.length, tmpQualities, baseQualities.length, l);
                baseQualities = tmpQualities;
            }
        }

        // Update SAM record
        rec.setCigar(new Cigar(cigarElements));
        rec.setAlignmentStart(alignmentStart);
        rec.setReadBases(readBases);
        rec.setBaseQualities(baseQualities);
        // Reset saved attributes
        Align.resetAttributes(rec, optFieldTags, optFieldValues);
        // Set new attributes
        if(space == SRMAUtil.Space.COLORSPACE) { 
            // set the XE attribute for colorError string
            rec.setAttribute("XE", new String(colorErrors));
        }
        rec.setAttribute("AS", bestAlignHeapNode.score);
        rec.setAttribute("XC", bestAlignHeapNode.alleleCoverageSum);
        rec.setAttribute("PG", programRecord.getId());
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
            return 0;
        }
        else {
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

    private static void clearAttributes(SAMRecord rec, List<String> optFieldTags, List<Object> optFieldValues)
    {
        ListIterator<String> iter = saveTags.listIterator();

        while(iter.hasNext()) {
            String tag = iter.next();
            Object attr = rec.getAttribute(tag);
            if(null != attr) {
                optFieldTags.add(tag);
                optFieldValues.add(attr);
            }
        }
        rec.clearAttributes();
    }

    private static void resetAttributes(SAMRecord rec, List<String> optFieldTags, List<Object> optFieldValues)
    {
        ListIterator<String> iterTags = optFieldTags.listIterator();
        ListIterator<Object> iterValues = optFieldValues.listIterator();

        while(iterTags.hasNext()) {
            rec.setAttribute(iterTags.next(), iterValues.next());
        }
    }

    private static class AlignData {
        // General info
        public SRMAUtil.Space space;
        public boolean strand;

        // Base space
        public String readBases = null;
        public String readBaseQualities=null; 

        // Color space
        public String readColors = null; // color space only, no adapter
        public String readColorQualities  = null; // color space only

        // Alignment info
        public boolean clippedStart = false;
        public boolean clippedEnd= false;

        public AlignData(SAMRecord rec)
            throws Exception
        {
            String read = null;

            this.strand = rec.getReadNegativeStrandFlag();

            // Get space
            read = (String)rec.getAttribute("CS");
            if(null == read) {
                // Use base space
                this.space = SRMAUtil.Space.NTSPACE;
            }
            else {
                // assumes CS and CQ are always in sequencing order
                this.space = SRMAUtil.Space.COLORSPACE;
            }

            // Note: the reads on the negative strand will not be complimented, only reversed
            // Get read and qualities
            if(this.space == SRMAUtil.Space.NTSPACE) {
                byte tmpRead[] = rec.getReadString().getBytes();
                byte tmpQualities[] = rec.getBaseQualityString().getBytes();
                // Reverse once
                if(this.strand) { // reverse
                    SAMRecordUtil.reverseArray(tmpRead);
                    SAMRecordUtil.reverseArray(tmpQualities);
                }
                this.readBases = new String(tmpRead);
                this.readBaseQualities = new String(tmpQualities);
                // Reverse again
                if(this.strand) { // reverse
                    SAMRecordUtil.reverseArray(tmpRead);
                    SAMRecordUtil.reverseArray(tmpQualities);
                }
                this.readColors = null;
                this.readColorQualities = null;
                if(this.readBaseQualities.length() <= 0) {
                    throw new Exception("Error.  The current alignment has no base qualities.");
                }
            }
            else {
                byte tmpRead[] = rec.getReadString().getBytes();
                // Reverse once
                if(this.strand) { // reverse
                    SAMRecordUtil.reverseArray(tmpRead);
                }
                this.readBases = new String(tmpRead);
                this.readBaseQualities = null;
                // Reverse again
                if(strand) { // reverse
                    SAMRecordUtil.reverseArray(tmpRead);
                }
                this.readColors = SRMAUtil.normalizeColorSpaceRead((String)rec.getAttribute("CS"));
                this.readColorQualities = (String)rec.getAttribute("CQ");
                // Some aligners include a quality value for the adapter.  A quality value
                // IMHO should not be given for an unobserved (assumed) peice of data.  Trim
                // the first quality in this case
                if(this.readColorQualities.length() > this.readColors.length()) { // trim the first quality
                    this.readColorQualities = readColorQualities.substring(1);
                }
                if(this.readColors.length() <= 0) {
                    throw new Exception("Error.  The current alignment has no colors.");
                }
                if(this.readColorQualities.length() <= 0) {
                    throw new Exception("Error.  The current alignment has no color qualities.");
                }
            }
            if(this.readBases.length() <= 0) {
                throw new Exception("Error.  The current alignment has no bases.");
            }

            // Deal with soft-clipping and hard-cliping
            {
                List<CigarElement> cigarElements = null;

                cigarElements = rec.getCigar().getCigarElements();
                CigarElement e1 = cigarElements.get(0); // first
                CigarElement e2 = cigarElements.get(cigarElements.size()-1); // last 

                // get the bases that are clipped
                if(CigarOperator.S == e1.getOperator()) {
                    int l = e1.getLength();
                    this.clippedStart = true;
                    if(this.strand) { // reverse
                        this.readBases = readBases.substring(0, this.readBases.length() - l);
                        if(this.space == SRMAUtil.Space.NTSPACE) {
                            this.readBaseQualities = this.readBaseQualities.substring(0, this.readBaseQualities.length() - l);
                        }
                        else { // adapter is not trimmed
                            this.readColors = this.readColors.substring(0, this.readColors.length() - l); 
                            this.readColorQualities = this.readColorQualities.substring(0, this.readColorQualities.length() - l); 
                        }
                    }
                    else { 
                        this.readBases = this.readBases.substring(l);
                        if(this.space == SRMAUtil.Space.NTSPACE) {
                            this.readBaseQualities = this.readBaseQualities.substring(l);
                        }
                        else { // adapter is trimmed
                            this.readColors =  this.readColors.substring(l); // skip the first l colors 
                            // the adapter is now the previously clipped base
                            this.readColors = Character.toString(SRMAUtil.getCompliment(this.readBases.charAt(l-1))) + this.readColors; 
                            this.readColors = SRMAUtil.normalizeColorSpaceRead(this.readColors); // re-normalize
                            this.readColorQualities = this.readColorQualities.substring(l); // skip the first l qualities
                        }
                    }
                }
                else if(CigarOperator.H == e1.getOperator()) {
                    this.clippedStart = true;
                    if(this.space == SRMAUtil.Space.COLORSPACE) {
                        int l = e1.getLength();
                        // the adapter is the last non-clipped base, inferred by the first non-clipped color
                        if(this.strand) { // reverse
                            // not trimed
                            this.readColors = this.readColors.substring(0, this.readColors.length() - l); 
                            this.readColorQualities = this.readColorQualities.substring(0, this.readColorQualities.length() - l); 
                        }
                        else {
                            this.readColors = this.readColors.substring(l); // skip the first l colors 
                            this.readColors = SRMAUtil.colorSpaceNextBase(this.readBases.charAt(0), this.readColors.charAt(0)) + this.readColors;
                            this.readColors = SRMAUtil.normalizeColorSpaceRead(this.readColors); // re-normalize
                            this.readColorQualities = this.readColorQualities.substring(l); // skip the first l qualities
                        }
                    }
                }

                if(CigarOperator.S == e2.getOperator()) {
                    int l = e2.getLength();
                    this.clippedEnd = true;
                    if(this.strand) { // reverse
                        this.readBases = this.readBases.substring(l);
                        if(this.space == SRMAUtil.Space.NTSPACE) {
                            this.readBaseQualities = this.readBaseQualities.substring(l);
                        }
                        else { // adapter is trimmed
                            this.readColors = this.readColors.substring(l); // skip the first l colors 
                            // the adapter is now the previously clipped base
                            this.readColors = Character.toString(SRMAUtil.getCompliment(this.readBases.charAt(l+1))) + this.readColors; 
                            this.readColors = SRMAUtil.normalizeColorSpaceRead(this.readColors); // re-normalize
                            this.readColorQualities = this.readColorQualities.substring(l); // skip the first
                        }
                    }
                    else {
                        this.readBases = this.readBases.substring(0, readBases.length() - l);
                        if(this.space == SRMAUtil.Space.NTSPACE) {
                            this.readBaseQualities = this.readBaseQualities.substring(0, this.readBaseQualities.length() - l);
                        }
                        else { // adapter is not trimmed
                            this.readColors = this.readColors.substring(0, readColors.length() - l);
                            this.readColorQualities = this.readColorQualities.substring(0, this.readColorQualities.length() - l);
                        }
                    }
                }
                else if(CigarOperator.H == e2.getOperator()) {
                    this.clippedEnd = true;
                    if(space == SRMAUtil.Space.COLORSPACE) {
                        int l = e2.getLength();
                        // the adapter is the last non-clipped base, inferred by the first non-clipped color
                        if(this.strand) { // reverse
                            this.readColors = this.readColors.substring(l); // skip the first l colors 
                            char b = SRMAUtil.getCompliment(this.readBases.charAt(0));
                            this.readColors = SRMAUtil.colorSpaceNextBase(b, this.readColors.charAt(0)) + this.readColors;
                            this.readColors = SRMAUtil.normalizeColorSpaceRead(this.readColors); // re-normalize
                            this.readColorQualities = this.readColorQualities.substring(l); // skip the first l qualities
                        }
                        else {
                            // not trimed
                            this.readColors = this.readColors.substring(0, readColors.length() - l);
                            this.readColorQualities = this.readColorQualities.substring(0, this.readColorQualities.length() - l);
                        }
                    }
                }
            }
        }
    }
}
