package srma;

import java.io.*;
import java.util.*;
import java.lang.Math;

public class AlleleCoverageCutoffs {
    Vector<Integer> minQ = null;
    int maxCoverage = 0;

    public AlleleCoverageCutoffs(int MINIMUM_ALLELE_COVERAGE, double MINIMUM_ALLELE_PROBABILITY, boolean QUIET) 
    {

        Vector<Vector<Integer>> binomialCoefficients = new Vector<Vector<Integer>>();
        binomialCoefficients.add(new Vector<Integer>()); // dummy for 0

        this.minQ = new Vector<Integer>();
        this.minQ.add(0); // dummy for 0

        this.maxCoverage = 0;
        int curCoverage = 1;
        int lastQ = 0;
        if(!QUIET) {
            System.err.println("Allele coverage cutoffs:");
        }
        while(lastQ < MINIMUM_ALLELE_COVERAGE) {
            // Dynamic programming for binomial co-efficients
            int k;
            Vector<Integer> v = new Vector<Integer>();
            for(k=0;k<=curCoverage;k++) {
                if(0 == k || curCoverage == k) {
                    v.add(1);
                }
                else {
                    v.add(binomialCoefficients.get(curCoverage-1).get(k-1) + binomialCoefficients.get(curCoverage-1).get(k));
                }
            }
            binomialCoefficients.add(v);

            // Get binomial quantile that satisfies probability
            double prob = 0.0;
            double p_term_2 = Math.pow(0.5, curCoverage); // p^i (1-p)^{n-i} => always the same same since p=0.5 is assumed
            int curQ = -1; // will always be incremented 
            do {
                curQ++;
                prob += p_term_2 * binomialCoefficients.get(curCoverage).get(curQ);
            }
            while(prob < MINIMUM_ALLELE_PROBABILITY);
            //System.err.println("cov="+curCoverage+" curQ="+curQ+" prob="+prob);
            lastQ = curQ;
            if(MINIMUM_ALLELE_COVERAGE < curQ) {
                this.minQ.add(MINIMUM_ALLELE_COVERAGE); 
                if(!QUIET) {
                    System.err.println("coverage: "+curCoverage+"\tminimum allele coverage: "+MINIMUM_ALLELE_COVERAGE);
                }
            }
            else {
                this.minQ.add(curQ); 
                if(!QUIET) {
                    System.err.println("coverage: "+curCoverage+"\tminimum allele coverage: "+curQ);
                }
            }
            this.maxCoverage = curCoverage;
            curCoverage++;
        }
        if(!QUIET) {
            System.err.println("coverage: >"+this.maxCoverage+"\tminimum allele coverage: "+MINIMUM_ALLELE_COVERAGE);
        }
    }

    public int getQ(int coverage)
    {
        if(coverage < 0) {
            return 0;
        }
        else if(this.maxCoverage < coverage) {
            return (int)this.minQ.get(this.maxCoverage);
        }
        else {
            return (int)this.minQ.get(coverage);
        }
    }
}
