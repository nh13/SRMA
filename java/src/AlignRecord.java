package srma;

import java.util.*;
import net.sf.samtools.*;

public class AlignRecord {
    public SAMRecord record;
    public Node node;
    public int fileIndex; // to which input/output file does this belong?

    public AlignRecord(SAMRecord record, Node node, int fileIndex)
    {
        this.record = record;
        this.node = node;
        this.fileIndex = fileIndex;
    }   

    public void setNode(Node node)
    {
        this.node = node;
    }
}  
