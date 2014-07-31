package test.beast.evolution.likelihood;

import org.junit.Test;

import test.beast.BEASTTestCase;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN98;
import beast.evolution.tree.Tree;
import junit.framework.TestCase;

/**
 * This test aims to check whether YN98 is compatible with the rest of BEAST2 program (before spending time to 
 * compile my plugin, construct .xml input file and test using real-world sequence alignment)
 * *
 */

public class YN98TreeLikelihoodTest extends TestCase {
	
    protected TreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new TreeLikelihood();
    }
	
    static public Alignment getCodonAlignment() throws Exception {
        Sequence human = new Sequence("human", "AGAAATATGTCTGATAAAAGAGTTACTTTG");
        Sequence chimp = new Sequence("chimp", "AGAAATATGTCTGATAAAAGAATTACTTTG");
        Sequence bonobo = new Sequence("bonobo", "AGAAATATGTCTGATAAAAGAATTACTTTG");
        Sequence gorilla = new Sequence("gorilla", "AGAAATATGTCTGATAAAAGAGTTACTTTG");
        Sequence orangutan = new Sequence("orangutan", "AGAAATATGTCTGACAAAAGAGTTACTTTG");
        Sequence siamang = new Sequence("siamang", "AGAAATACGTCTGACGAAAGAGTTACTTTG");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", chimp, "sequence", bonobo, "sequence", gorilla, "sequence", orangutan, "sequence", siamang,
                "userDataType", codon);

        return data;
    }
    
    static public Alignment getCodonAlignment3Sites() throws Exception {
        Sequence human = new Sequence("human", "AGAAATATG");
        Sequence chimp = new Sequence("chimp", "AGAAATATG");
        Sequence bonobo = new Sequence("bonobo", "AGAAATATG");
        Sequence gorilla = new Sequence("gorilla", "AGAAATATG");
        Sequence orangutan = new Sequence("orangutan", "AGAAATATG");
        Sequence siamang = new Sequence("siamang", "AGAAATACG");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", chimp, "sequence", bonobo, "sequence", gorilla, "sequence", orangutan, "sequence", siamang,
                "userDataType", codon);

        return data;
    }
    
    @Test
    public void testYN98Likelihood() throws Exception {
    	
        Alignment data = getCodonAlignment();
        Tree tree = BEASTTestCase.getTree(data);

        RealParameter f = new RealParameter(new Double[]{0.2, 0.2, 0.3, 0.3});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);

        YN98 yn98 = new YN98();
        yn98.initByName("kappa", "1.5", "omega", "0.9", "nucleoFrequencies", nucleoFrequencies);

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("substModel", yn98);

        TreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);

        /*
        double fLogP = 0;
        fLogP = likelihood.calculateLogP();
        assertEquals(fLogP, -1789.912401996943, BEASTTestCase.PRECISION);

        likelihood.initByName("useAmbiguities", true, "data", data, "tree", tree, "siteModel", siteModel);
        fLogP = likelihood.calculateLogP();
        assertEquals(fLogP, -1789.912401996943, BEASTTestCase.PRECISION);
        */
    }
    
    @Test
    public void testYN98EachSiteLikelihood() throws Exception {
    	
        Alignment data = getCodonAlignment3Sites();
        Tree tree = BEASTTestCase.getTree(data);
        
        //use FilteredAlignment to create a sub data for each site
        FilteredAlignment site1 = new FilteredAlignment();
        //Note that "2" refers to second codon !!!!! sequence position starts at "1" instead of "0"!!!!
        site1.initByName("data", data, "filter", "1");
        FilteredAlignment site2 = new FilteredAlignment();
        site2.initByName("data", data, "filter", "2");
        FilteredAlignment site3 = new FilteredAlignment();
        site3.initByName("data", data, "filter", "3");
        
        RealParameter f = new RealParameter(new Double[]{0.2, 0.2, 0.3, 0.3});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);

        YN98 site1_yn98 = new YN98();
        YN98 site2_yn98 = new YN98();
        YN98 site3_yn98 = new YN98();
        site1_yn98.initByName("kappa", "1.5", "omega", "0.9", "nucleoFrequencies", nucleoFrequencies);
        site2_yn98.initByName("kappa", "2.5", "omega", "1.5", "nucleoFrequencies", nucleoFrequencies);
        site3_yn98.initByName("kappa", "5", "omega", "3.4", "nucleoFrequencies", nucleoFrequencies);
        
        SiteModel site1_siteModel = new SiteModel();
        site1_siteModel.initByName("substModel", site1_yn98);
        SiteModel site2_siteModel = new SiteModel();
        site2_siteModel.initByName("substModel", site2_yn98);
        SiteModel site3_siteModel = new SiteModel();
        site3_siteModel.initByName("substModel", site3_yn98);
        
        TreeLikelihood site1_likelihood = newTreeLikelihood();
        site1_likelihood.initByName("data", site1, "tree", tree, "siteModel", site1_siteModel);
        TreeLikelihood site2_likelihood = newTreeLikelihood();
        site2_likelihood.initByName("data", site2, "tree", tree, "siteModel", site2_siteModel);
        TreeLikelihood site3_likelihood = newTreeLikelihood();
        site3_likelihood.initByName("data", site3, "tree", tree, "siteModel", site3_siteModel);
    
    }
}
