/*
 * LICENSE to be determined
 */
package srma;

import java.io.*;
import java.util.*;

public class NodeRecord {
    public Node node;
    public int coverage;

    public NodeRecord(Node node, int coverage) 
    {
        this.node = node;
        this.coverage = coverage;
    }
}
