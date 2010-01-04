/*
 * LICENSE to be determined
 */
package srma;

import java.util.*;
import srma.*;

public class SRMAUtil {
    // These should be in a util class
    public enum Space { NTSPACE, COLORSPACE }
    public static final byte COLORSPACE_ADAPTOR = 'A';
    public static final byte DNA[] = {'A', 'C', 'G', 'T', 'N'};
    public static final byte COLORS[] = {'0', '1', '2', '3', '4'};

    public static int CHAR2QUAL(byte qual)
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

    public static byte colorSpaceNextBase(byte base, byte color) throws Exception
    {
        int start=0, by=0, result=0;

        switch(base) {
            case 'A':
                start=0; by=1; break;
            case 'C':
                start=1; by=-1; break;
            case 'G':
                start=2; by=1; break;
            case 'T':
                start=3; by=-1; break;
            case 'N':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base");
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
            case '4':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        if(result < 0) {
            return DNA[4 - ( (-1*result) % 4)];
        }
        else {
            return DNA[result % 4];
        }
    }
    
    public static byte colorSpaceEncode(byte b1, byte b2) throws Exception
    {
        int start=0, by=0, result=0;

        switch(b1) {
            case 'A':
                start=0; by=1; break;
            case 'C':
                start=1; by=-1; break;
            case 'G':
                start=2; by=1; break;
            case 'T':
                start=3; by=-1; break;
            case 'N':
                return 'N'; 
            default:
                throw new Exception("Error: could not understand the base");
        }

        switch(b2) {
            case 'A':
                result = start; break;
            case 'C':
                result = start + by; break;
            case 'G':
                result = start + 2*by; break;
            case 'T':
                result = start + 3*by; break;
            case 'N':
                return 'N'; 
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

    public static void normalizeColorSpaceRead(byte read[])
        throws Exception
    {
        if(read.length < 2) {
            throw new Exception("Read was too short");
        }
        else if(read[0] != COLORSPACE_ADAPTOR) {
            byte base = colorSpaceNextBase(read[0], read[1]);
            read[0] = COLORSPACE_ADAPTOR;
            read[1] = colorSpaceEncode(read[0], base);
        }
        // Remove adaptor
        int i;
        byte tmp[] = new byte[read.length-1];
        for(i=0;i<read.length-1;i++) {
            tmp[i] = read[i+1];
        }
        read = tmp;
    }
}
