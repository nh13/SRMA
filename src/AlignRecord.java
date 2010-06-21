package srma;

import java.util.*;
import net.sf.samtools.*;

public class AlignRecord {
    public SAMRecord record;
    public Node node;

    public AlignRecord(SAMRecord record, Node node)
    {
        this.record = record;
        this.node = node;
    }   
}  
