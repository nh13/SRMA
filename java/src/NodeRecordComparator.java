/*
 *  * LICENSE to be determined
 *   */
package srma;

import java.io.*;
import java.util.*;

public class NodeRecordComparator implements Comparator<NodeRecord> 
{
    private NodeComparator nodeComparator;

    public NodeRecordComparator()
    {
        this.nodeComparator = new NodeComparator();
    }

    public int compare(NodeRecord o1, NodeRecord o2)
    {
        return this.nodeComparator.compare(o1.node, o2.node);
    }

    public boolean equals(NodeRecord o1, NodeRecord o2)
    {
        return (0 == this.compare(o1, o2));
    }
}
