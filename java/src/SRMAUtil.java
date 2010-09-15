/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import srma.*;

public class SRMAUtil {
    // These should be in a util class
    public enum Space { NTSPACE, COLORSPACE }
    public static final char COLORSPACE_ADAPTOR = 'A';
    public static final char DNA[] = {'A', 'C', 'G', 'T', 'N'};
    public static final char COLORS[] = {'0', '1', '2', '3', '4'};

    public static int CHAR2QUAL(char qual)
    {
        int q = (((int)qual) - 33);
        if(q < 0) {
            return 1;
        }
        else if(255 < q) {
            return 255;
        }
        else {
            return q;
        }
    }

    public static char QUAL2CHAR(int qual)
    {
        if(93 < qual) {
            return (char)(126);
        }
        else if(qual < 0) {
            return (char)(33);
        }
        return (char)(qual + 33);
    }

    public static char colorSpaceNextBase(char base, char color) throws Exception
    {
        int start=0, by=0, result=0;

        switch(base) {
            case 'A':
            case 'a':
                start=0; by=1; break;
            case 'C':
            case 'c':
                start=1; by=-1; break;
            case 'G':
            case 'g':
                start=2; by=1; break;
            case 'T':
            case 't':
                start=3; by=-1; break;
            case 'N':
            case 'n':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base: " + base);
        }

        switch(color) {
            case '0':
                result = start; break;
            case '1':
                result = start + by; break;
            case '2':
                result = start + 2*by; break;
            case '3':
                result = start + 3*by; break;
            case '.':
            case '4':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the color");
        }

        if(result < 0) {
            return DNA[4 - ( (-1*result) % 4)];
        }
        else {
            return DNA[result % 4];
        }
    }

    public static char colorSpaceEncode(char b1, char b2) throws Exception
    {
        int start=0, by=0, result=0;

        switch(b1) {
            case 'A':
            case 'a':
                start=0; by=1; break;
            case 'C':
            case 'c':
                start=1; by=-1; break;
            case 'G':
            case 'g':
                start=2; by=1; break;
            case 'T':
            case 't':
                start=3; by=-1; break;
            case 'N':
            case 'n':
                return '4'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        switch(b2) {
            case 'A':
            case 'a':
                result = start; break;
            case 'C':
            case 'c':
                result = start + by; break;
            case 'G':
            case 'g':
                result = start + 2*by; break;
            case 'T':
            case 't':
                result = start + 3*by; break;
            case 'N':
            case 'n':
                return '4'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        if(result < 0) {
            return COLORS[4 - ( (-1*result) % 4)];
        }
        else {
            return COLORS[result % 4];
        }
    }

    public static String normalizeColorSpaceRead(String read)
        throws Exception
    {
        // Remove the adaptor
        String ret = new String(read.substring(1));;

        if(read.length() < 2) {
            throw new Exception("Read was too short");
        }
        else if(read.charAt(0) != COLORSPACE_ADAPTOR) {
            ret = replaceCharAt(ret, 0, colorSpaceEncode(COLORSPACE_ADAPTOR, 
                        colorSpaceNextBase(read.charAt(0), read.charAt(1))));
        }
        return ret;
    }

    public static String replaceCharAt(String s, int pos, char c) {
        return s.substring(0,pos) + c + s.substring(pos+1);
    }

    public static char getCompliment(char base) {
        switch(base) {
            case 'A':
                return 'T';
            case 'a':
                return 't';
            case 'C':
                return 'G';
            case 'c':
                return 'g';
            case 'G':
                return 'C';
            case 'g':
                return 'c';
            case 'T':
                return 'A';
            case 't':
                return 'a';
            default:
                return base;
        }
    }
}
