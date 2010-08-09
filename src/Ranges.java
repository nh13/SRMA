package srma;

import java.io.*;
import java.util.*;
import net.sf.picard.reference.*;
import net.sf.samtools.*;

// TODO:
// - need to check references are in order
// - better error messages
// - state how many ranges were found

public class Ranges {

    private LinkedList<Range> ranges = null;

    public Ranges(File file, SAMSequenceDictionary referenceDictionary, int offset)
    {
        BufferedReader br = null;
        String line = null;
        int i, lineNumber = 1;
        Map<String, Integer> hm = new HashMap<String, Integer>();

        try {
            // open
            br = new BufferedReader(new FileReader(file));

            // init
            this.ranges = new LinkedList<Range>();
            for(i=0;i<referenceDictionary.size();i++) {
                hm.put(referenceDictionary.getSequence(i).getSequenceName(), new Integer(i));
            }

            // read the file
            while(null != (line = br.readLine())) {
                this.addRange(line, lineNumber, hm, referenceDictionary, offset);
                lineNumber++;
            }

            // close
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public Ranges(File file, SAMSequenceDictionary referenceDictionary)
    {
        this(file, referenceDictionary, 0);
    }

    public Ranges(String range, SAMSequenceDictionary referenceDictionary, int offset)
        throws Exception
    {
        Map<String, Integer> hm = new HashMap<String, Integer>();
        int i, colon, dash, referenceIndex;
        int startPosition, endPosition;

        // init
        this.ranges = new LinkedList<Range>();
        for(i=0;i<referenceDictionary.size();i++) {
            hm.put(referenceDictionary.getSequence(i).getSequenceName(), new Integer(i));
        }

        // get delimiters
        colon = range.indexOf(':');
        dash = range.indexOf('-', colon+1);

        // Only chromosome name
        if(colon < 0) { // contig only
            referenceIndex = (int)hm.get(range);
            startPosition = 1;
            endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
        }
        else {
            if(colon <= 0 || (0 < dash && (dash - colon) <= 1) || (0 < dash && (range.length() - dash) <= 1)) {
                throw new Exception("RANGE was improperly specified.");
            }   
                
            String chrName = range.substring(0, colon);
            if(null == hm.get(chrName)) {
                throw new Exception("Could not find reference name [" + chrName + "] in RANGE");
            }
            referenceIndex = (int)hm.get(chrName);

            if(dash <= 0) {
                startPosition = Integer.parseInt(range.substring(colon+1));
                endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
            }
            else {
                startPosition = Integer.parseInt(range.substring(colon+1, dash));
                endPosition = Integer.parseInt(range.substring(dash+1));
            }
        }

        if(startPosition <= 0 || referenceDictionary.getSequence(referenceIndex).getSequenceLength() < startPosition) {
            throw new Exception("startPosition was out of bounds in RANGE");
        }
        else if(endPosition <= 0 || referenceDictionary.getSequence(referenceIndex).getSequenceLength() < endPosition) {
            throw new Exception("endPosition was out of bounds in RANGE");
        }
        else if(endPosition < startPosition) {
            throw new Exception("endPosition < startPosition in RANGE");
        }

        startPosition -= offset;
        if(startPosition <= 0) {
            startPosition = 1;
        }
        endPosition += offset;
        if(referenceDictionary.getSequence(referenceIndex).getSequenceLength() < endPosition) {
            endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
        }

        this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
    }

    public Ranges(String range, SAMSequenceDictionary referenceDictionary)
        throws Exception
    {
        this(range, referenceDictionary, 0);
    }

    private void addRange(String line, int lineNumber, Map<String, Integer> m, SAMSequenceDictionary referenceDictionary, int offset)
        throws Exception
    {
        StringTokenizer st = new StringTokenizer(line);

        int referenceIndex=-1;
        int startPosition=-1, endPosition=-1;
        int i=0;

        while(st.hasMoreTokens()) {
            if(0 == i) {
                Integer tmpInteger = -1;
                String chrName = null;
                try {
                    chrName = new String(st.nextToken());
                    tmpInteger = (int)m.get(chrName);
                } catch(java.lang.NullPointerException e) {
                    throw new Exception("Could not find reference name ["+chrName+"] in RANGES on line " + lineNumber);
                }
                if(null == tmpInteger) {
                    throw new Exception("Could not find reference name ["+chrName+"] in RANGES on line " + lineNumber);
                }
                referenceIndex = (int)tmpInteger;
            }
            else if(1 == i) {
                startPosition = Integer.parseInt(new String(st.nextToken()));
                if(startPosition <= 0 || referenceDictionary.getSequence(referenceIndex).getSequenceLength() < startPosition) {
                    throw new Exception("start position was out of bounds in RANGES on line " + lineNumber);
                }
                // add offset
                startPosition -= offset;
                if(startPosition <= 0) {
                    startPosition = 1;
                }
            }
            else if(2 == i) {
                endPosition = Integer.parseInt(new String(st.nextToken()));
                if(endPosition <= 0 || referenceDictionary.getSequence(referenceIndex).getSequenceLength() < endPosition) {
                    throw new Exception("end position was out of bounds in RANGES on line " + lineNumber);
                }
                // add offset
                endPosition += offset;
                if(referenceDictionary.getSequence(referenceIndex).getSequenceLength() < endPosition) {
                    endPosition = referenceDictionary.getSequence(referenceIndex).getSequenceLength();
                }
            }
            else {
                new Exception("Too many entries in RANGES on line " + lineNumber);
            }
            i++;
        }
        if(3 != i) {
            new Exception("Too few entries in RANGES on line " + lineNumber);
        }
        if(endPosition < startPosition) {
            throw new Exception("End position < start position in RANGES on line " + lineNumber);
        }

        if(null != this.ranges.peek()) {
            Range last = this.ranges.getLast();
            if(referenceIndex < last.referenceIndex) {
                throw new Exception("Ranges must be in sorted order in RANGES (line numbers " + (lineNumber-1) + "-" + lineNumber);
            }
            else if(referenceIndex == last.referenceIndex) { // sam reference
                if(startPosition < last.startPosition) {
                    throw new Exception("Ranges must be in sorted order in RANGES (line numbers " + (lineNumber-1) + "-" + lineNumber);
                }
                else if(last.endPosition + 1 < startPosition) {
                    // just add
                    this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
                }
                else if(last.endPosition < endPosition) {
                    // merge over-lapping
                    last.endPosition = endPosition;
                }
            }
            else {
                // just add
                this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
            }
        }
        else {
            // just add
            this.ranges.add(new Range(referenceIndex, startPosition, endPosition));
        }
    }

    public Iterator<Range> iterator() 
    {
        return this.ranges.iterator();
    }

    public Range get(int i)
    {
        return this.ranges.get(i);
    }

    public int size()
    {
        return this.ranges.size();
    }
}
