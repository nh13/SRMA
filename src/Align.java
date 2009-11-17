/*
 * LICENSE to be determined
 */
package srma;

import net.sf.samtools.*;
import java.util.*;

public class Align {

	public Align(Graph graph, SAMRecord rec)
	{
		if(rec.getReadNegativeStrandFlag()) {
			// Reverse
			System.err.println("Not implemented\n");
		}
		else {
			// Forward
			AlignForwardStrand(graph, rec);
		}
	}

	private void AlignForwardStrand(Graph graph, SAMRecord rec)
	{
		Node curNode=null;

		// TODO implement this
		//curNode = graph.getStartNode(rec.getAlignmentStart());
	}

}
