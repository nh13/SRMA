/*
 *  * LICENSE to be determined
 *   */
package srma;

import java.io.*;
import java.util.*;

public class NodeComparator implements Comparator<Node> 
{

    // TODO: optimize this
    static public int compareNode(Node o1, Node o2)
    {
        if(o1.contig == o2.contig 
                && o1.position == o2.position 
                && o1.offset == o2.offset 
                && o1.type == o2.type 
                && o1.base == o2.base) 
        {
            return 0;
        }
        else if(o1.contig < o2.contig 
                || (o1.contig == o2.contig 
                    && o1.position < o2.position) 
                || (o1.contig == o2.contig 
                    && o1.position == o2.position 
                    && o1.offset < o2.offset) 
                || (o1.contig == o2.contig 
                    && o1.position == o2.position 
                    && o1.offset == o2.offset 
                    && o1.type < o2.type) 
                || (o1.contig == o2.contig 
                    && o1.position == o2.position 
                    && o1.offset == o2.offset 
                    && o1.type == o2.type 
                    && o1.base < o2.base)) 
        {
            return -1;
        }
        else {
            return 1;
        }   
    }
    
    public int compare(Node o1, Node o2)
    {
        if(o1.contig == o2.contig 
                && o1.position == o2.position 
                && o1.offset == o2.offset 
                && o1.type == o2.type 
                && o1.base == o2.base) 
        {
            return 0;
        }
        else if(o1.contig < o2.contig 
                || (o1.contig == o2.contig 
                    && o1.position < o2.position) 
                || (o1.contig == o2.contig 
                    && o1.position == o2.position 
                    && o1.offset < o2.offset) 
                || (o1.contig == o2.contig 
                    && o1.position == o2.position 
                    && o1.offset == o2.offset 
                    && o1.type < o2.type) 
                || (o1.contig == o2.contig 
                    && o1.position == o2.position 
                    && o1.offset == o2.offset 
                    && o1.type == o2.type 
                    && o1.base < o2.base)) 
        {
            return -1;
        }
        else {
            return 1;
        }   
    }
    
    public boolean equals(Node o1, Node o2)
    {
        return (0 == this.compare(o1, o2));
    }
}
