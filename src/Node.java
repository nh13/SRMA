/*
 * LICENSE to be determined
 */
package srma;

import java.io.*;
import java.util.*;
    
public class Node {

    public enum NodeType { MATCH, MISMATCH, INSERTION, DELETION }

    byte base; // [acgtnACGTN]
    NodeType type;
    int contig;
    int position;
    List<Node> next; // downstream nodes
    List<Node> prev; // upstream nodes

    public Node(byte base, NodeType type, int contig, int position, Node node)
    {
        this.base = base;
        this.type = type;
        this.contig = contig;
        this.position = position;
        this.next = new ArrayList<Node>();
        this.prev = new ArrayList<Node>();
    }

    public void addToNext(Node node)
    {
        if(!this.next.contains(node)) {
            this.next.add(node);
        }
    }

    public void addToPrev(Node node) 
    {
        if(!this.prev.contains(node)) {
            this.prev.add(node);
        }
    }

    public void print(PrintStream out) 
    {
        out.println("" + this.base + "\t" + this.type + "\t" + this.contig + "\t" + 
                this.position + "\t" + this.prev.size() + "\t" + this.next.size());
    }

    public void print()
    {
        this.print(System.out);
    }
}
