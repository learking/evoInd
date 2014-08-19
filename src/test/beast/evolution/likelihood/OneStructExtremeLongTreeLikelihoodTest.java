package test.beast.evolution.likelihood;

import org.junit.Test;

import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.datatype.UserDataType;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.OneStruct;
import beast.evolution.tree.Tree;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class OneStructExtremeLongTreeLikelihoodTest extends TestCase {
    protected TreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new TreeLikelihood();
    }
    
    static public Tree getTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "(human:10000,chimp:10000);",
                "IsLabelledNewick", true);
        return tree;
    }

    static public Alignment getCodonAlignment3Sites() throws Exception {
        Sequence human = new Sequence("human", "GTTATGAAC");
        Sequence chimp = new Sequence("chimp", "GTAATTAAA");

        UserDataType codon = new UserDataType();
        codon.initByName("states", 61, "codelength", 3, "codeMap", "AAA=0, AAC=1, AAG=2, AAT=3, ACA=4, ACC=5, ACG=6, ACT=7, AGA=8, AGC=9, AGG=10, AGT=11, ATA=12, ATC=13, ATG=14, ATT=15, CAA=16, CAC=17, CAG=18, CAT=19, CCA=20, CCC=21, CCG=22, CCT=23, CGA=24, CGC=25, CGG=26, CGT=27, CTA=28, CTC=29, CTG=30, CTT=31, GAA=32, GAC=33, GAG=34, GAT=35, GCA=36, GCC=37, GCG=38, GCT=39, GGA=40, GGC=41, GGG=42, GGT=43, GTA=44, GTC=45, GTG=46, GTT=47, TAC=48, TAT=49, TCA=50, TCC=51, TCG=52, TCT=53, TGC=54, TGG=55, TGT=56, TTA=57, TTC=58, TTG=59, TTT=60");

        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", chimp, "userDataType", codon);

        return data;
    }
 
    @Test
    public void testOneStructExtremeLongTreeLikelihood() throws Exception {
        Alignment data = getCodonAlignment3Sites();
        Tree tree = getTree(data);
        
        //use FilteredAlignment to create a sub data for each site
        FilteredAlignment site1 = new FilteredAlignment();
        //Note that "2" refers to second codon !!!!! sequence position starts at "1" instead of "0"!!!!
        site1.initByName("data", data, "filter", "1");
        FilteredAlignment site2 = new FilteredAlignment();
        site2.initByName("data", data, "filter", "2");
        FilteredAlignment site3 = new FilteredAlignment();
        site3.initByName("data", data, "filter", "3");
        
        RealParameter f = new RealParameter(new Double[]{0.2, 0.1, 0.3, 0.4});
        Frequencies nucleoFrequencies = new Frequencies();
        nucleoFrequencies.initByName("frequencies", f, "estimate", false);
        
        RealParameter codonProb1 = new RealParameter(new Double[]{
        		9.726911e-02,3.710371e-02,1.013584e-01,3.560677e-02,8.560782e-03,1.038485e-02,3.445223e-03,7.038088e-03,1.549992e-02,1.317710e-02,
        		1.490684e-02,8.523435e-03,9.480222e-05,2.709963e-04,1.067433e-03,2.339801e-04,3.958545e-02,7.734925e-03,4.124965e-02,7.422862e-03,
        		2.429324e-02,2.643607e-02,7.080957e-03,2.436986e-02,1.015356e-02,1.494012e-02,1.499864e-02,6.872558e-03,1.504320e-04,4.103283e-04,
        		8.159234e-04,3.016074e-04,9.984975e-02,6.072661e-02,1.040475e-01,5.827662e-02,2.610255e-03,4.707261e-03,9.777681e-04,3.397439e-03,
        		1.292552e-02,1.502855e-02,1.017521e-02,7.776232e-03,2.489455e-04,4.652802e-04,9.462066e-04,3.815349e-04,1.233373e-03,1.183613e-03,
        		7.647608e-03,1.153448e-02,2.790438e-03,1.000683e-02,1.062598e-04,3.968364e-04,1.019728e-04,1.780473e-04,3.360625e-04,2.636806e-04,
        		3.225042e-04});
    	
        RealParameter codonProb2 = new RealParameter(new Double[]{
        		0.0708144084,0.0365015955,0.0737915015,0.0350289510,0.0134057090,0.0162621033,0.0053950276,0.0110212535,0.0160822056,0.0179713013,
        		0.0154668482,0.0116245039,0.0006124190,0.0017506265,0.0039252937,0.0015115031,0.0377442932,0.0123967368,0.0393310928,0.0118965946,
        		0.0239040657,0.0260125691,0.0069675222,0.0239794665,0.0105349944,0.0155013786,0.0155620904,0.0071307393,0.0007869281,0.0021464767,
        		0.0042681935,0.0015777448,0.0767535373,0.0547719836,0.0799803160,0.0525622265,0.0063591373,0.0114678898,0.0023820515,0.0082768857,
        		0.0184226259,0.0214200491,0.0145026255,0.0110833894,0.0011512075,0.0021516120,0.0043755770,0.0017643454,0.0043196665,0.0041453910,
        		0.0104300256,0.0157310551,0.0038056787,0.0136476063,0.0007463889,0.0019162404,0.0007162761,0.0009313871,0.0019898101,0.0013793451,
        		0.0019095319});
        
        RealParameter codonProb3 = new RealParameter(new Double[]{
        		0.0735727266,0.0367896284,0.0766657815,0.0353053634,0.0128928858,0.0156400115,0.0051886457,0.0105996456,0.0161223411,0.0175263141,
        		0.0155054480,0.0113366696,0.0005105762,0.0014595043,0.0034638089,0.0012601461,0.0381622377,0.0118946590,0.0397666080,0.0114147730,
        		0.0240921782,0.0262172745,0.0070223530,0.0241681725,0.0105612860,0.0155400646,0.0156009279,0.0071485351,0.0006701706,0.0018280013,
        		0.0036349164,0.0013436529,0.0793074247,0.0556896345,0.0826415709,0.0534428551,0.0058493520,0.0105485574,0.0021910924,0.0076133626,
        		0.0178867974,0.0207970395,0.0140808116,0.0107610252,0.0009926640,0.0018552934,0.0037729753,0.0015213608,0.0038306006,0.0036760563,
        		0.0101717678,0.0153415386,0.0037114464,0.0133096780,0.0006170363,0.0016451456,0.0005921421,0.0007931960,0.0016735485,0.0011746899,
        		0.0016060297});
        
        OneStruct site1_oneStruct = new OneStruct();
        OneStruct site2_oneStruct = new OneStruct();
        OneStruct site3_oneStruct = new OneStruct();
        
        site1_oneStruct.initByName("kappa", "1", "omega", "1", "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb1);
        site2_oneStruct.initByName("kappa", "1", "omega", "1", "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb2);
        site3_oneStruct.initByName("kappa", "1", "omega", "1", "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb3);
        
        SiteModel site1_siteModel = new SiteModel();
        site1_siteModel.initByName("substModel", site1_oneStruct);
        SiteModel site2_siteModel = new SiteModel();
        site2_siteModel.initByName("substModel", site2_oneStruct);
        SiteModel site3_siteModel = new SiteModel();
        site3_siteModel.initByName("substModel", site3_oneStruct);
        
        TreeLikelihood site1_likelihood = newTreeLikelihood();
        site1_likelihood.initByName("data", site1, "tree", tree, "siteModel", site1_siteModel);
        TreeLikelihood site2_likelihood = newTreeLikelihood();
        site2_likelihood.initByName("data", site2, "tree", tree, "siteModel", site2_siteModel);
        TreeLikelihood site3_likelihood = newTreeLikelihood();
        site3_likelihood.initByName("data", site3, "tree", tree, "siteModel", site3_siteModel); 
        
        double totalLogP = site1_likelihood.calculateLogP() + site2_likelihood.calculateLogP() + site3_likelihood.calculateLogP();
        System.out.format("total logP is:" + "%f%n", totalLogP);
        //total logP is:-34.116570
    }
    
}
