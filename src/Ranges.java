package srma;

import java.io.*;
import java.util.*;
import net.sf.picard.reference.*;

// TODO:
// - need to check references are in order
// - better error messages
// - state how many ranges were found

public class Ranges {

    private LinkedList<Range> ranges = null;

    public Ranges(File file, List<ReferenceSequence> referenceSequences)
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
            for(i=0;i<referenceSequences.size();i++) {
                hm.put(referenceSequences.get(i).getName(), new Integer(i));
            }

            // read the file
            while(null != (line = br.readLine())) {
                this.addRange(line, lineNumber, hm);
                lineNumber++;
            }

            // close
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void addRange(String line, int lineNumber, Map<String, Integer> m)
        throws Exception
    {
        StringTokenizer st = new StringTokenizer(line);

        int referenceIndex=-1;
        int startPosition=-1, endPosition=-1;
        int i=0;

        while(st.hasMoreTokens()) {
            if(0 == i) {
                Integer tmpInteger = (int)m.get(new String(st.nextToken()));
                if(null == tmpInteger) {
                    throw new Exception("Could not find reference name on line " + lineNumber);
                }
                referenceIndex = (int)tmpInteger;
            }
            else if(1 == i) {
                startPosition = Integer.parseInt(new String(st.nextToken()));
            }
            else if(2 == i) {
                endPosition = Integer.parseInt(new String(st.nextToken()));
            }
            else {
                new Exception("Too many entries in ranges on line " + lineNumber);
            }
            i++;
        }
        if(3 != i) {
            new Exception("Too few entries in ranges on line " + lineNumber);
        }
        if(endPosition < startPosition) {
            throw new Exception("End position < start position on line " + lineNumber);
        }

        if(null != this.ranges.peek()) {
            Range last = this.ranges.getLast();
            if(referenceIndex < last.referenceIndex) {
                throw new Exception("Ranges must be in sorted order (line numbers " + (lineNumber-1) + "-" + lineNumber);
            }
            else if(referenceIndex == last.referenceIndex) { // sam reference
                if(startPosition < last.startPosition) {
                    throw new Exception("Ranges must be in sorted order (line numbers " + (lineNumber-1) + "-" + lineNumber);
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
