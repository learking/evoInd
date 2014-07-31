package test.beast.evolution.substmodel;

import java.util.Arrays;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.YN98;
import junit.framework.TestCase;

public class YN98Test extends TestCase {
 
	public interface Instance {
	    //frequencies for pi_A, pi_T, pi_C and pi_G
        Double[] getPi();

        Double getKappa();

        Double getOmega();
        
        //double getDistance();

        //double[] getExpectedResult();
    }
	
    protected Instance test0 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.25, 0.25, 0.25, 0.25};
        }

        public Double getKappa() {
            return 2.0;
        }

        public Double getOmega() {
            return 2.0;
        }


    };
    
    protected Instance test1 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.2, 0.2, 0.3, 0.3};
        }

        public Double getKappa() {
            return 1.5;
        }

        public Double getOmega() {
            return 0.9;
        }

    };
    
    Instance[] all = {test0, test1};
    
    Double getTotalSubstRate(double[][] rateM, double[] diagMatrix) {
        double fSubst = 0.0;
        for (int i = 0; i < rateM.length; i++)
        	//please see Yang and Nielsen (2008) for explanation (symmMatrix times diagonal, while diagonal matrix represents stationary for codons)
            fSubst += -rateM[i][i] * diagMatrix[i];
    	return fSubst;
    }
    
    public void testSetupRateMatrix() throws Exception {
        for (Instance test : all) {

            RealParameter f = new RealParameter(test.getPi());

            Frequencies nucleoFrequencies = new Frequencies();
            nucleoFrequencies.initByName("frequencies", f, "estimate", false);
            
            System.out.println(Arrays.toString(nucleoFrequencies.getFreqs()));

            YN98 yn98 = new YN98();
            yn98.initByName("kappa", test.getKappa().toString(), "omega", test.getOmega().toString(), "nucleoFrequencies", nucleoFrequencies);

            //set up matrices before test
            yn98.prepareMatricesForTest();
            
            final double[][] rateM = yn98.getRateMatrix();
            System.out.print("rateMatrix:");
            System.out.println(Arrays.deepToString(rateM));
            //check entries in our rateMatrix (currently, it looks like YN98)
            
            //dimension is 61X61
            assertEquals(rateM.length, 61);
            assertEquals(rateM[0].length, 61);
     
            // make sure symmetric matrix is indeed symmetric
            final double[][] symmM = yn98.getSymmMatrix();
            System.out.print("symmMatrix:");
            System.out.println(Arrays.deepToString(symmM));
            for (int i = 0; i < symmM.length; i++) {
            	for (int j = 0; j < symmM[0].length; j++) {
        			if (i < j){
        				assertEquals(symmM[j][i],symmM[i][j], 1e-10);
        			}
            	}
            }
            
            //expected substitution should be one per unit time
            final double[] diagM = yn98.getDiagMatrix();
            System.out.print("diagMatrix:");
            System.out.println(Arrays.toString(diagM));
            assertEquals(getTotalSubstRate(rateM, diagM), 1.0, 1e-10);
            
        }
    }
    
}
