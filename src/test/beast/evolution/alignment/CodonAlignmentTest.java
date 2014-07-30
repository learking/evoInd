package test.beast.evolution.alignment;

import org.junit.Test;

import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import junit.framework.TestCase;

public class CodonAlignmentTest extends TestCase {
    static public Alignment getAlignment() throws Exception {
        Sequence human = new Sequence("human", "AAAACCCCGGGGAAA");
        Sequence chimp = new Sequence("chimp", "ACGTACGTACGTACG");

        UserDataType dataType = new UserDataType();
        dataType.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");
        
        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", chimp,
        		"userDataType", dataType
        );
        return data;
    }

    @Test
    public void testRangeFiltered() throws Exception {
    	//construct codon alignment
        Alignment data = getAlignment();
        assertEquals(5, data.getSiteCount());
        assertEquals(4, data.getPatternCount());
        
        //try to get one site from our original alignment
        FilteredAlignment data2 = new FilteredAlignment();
        //Note that "2" refers to second codon !!!!! sequence position starts at "1" instead of "0"!!!!
        data2.initByName("data", data, "filter", "2");
        assertEquals(1, data2.getSiteCount());
        assertEquals(1, data2.getPatternCount());
        
        //iTaxon = 1 is chimp
        System.out.println(alignmentToString(data2, 1));
    }
    
    String alignmentToString(Alignment data, int iTaxon) throws Exception {
        int[] nStates = new int[data.getSiteCount()];
        for (int i = 0; i < data.getSiteCount(); i++) {
            int iPattern = data.getPatternIndex(i);
            int[] sitePattern = data.getPattern(iPattern);
            nStates[i] = sitePattern[iTaxon];
        }
        return data.getDataType().state2string(nStates);
    }
}
