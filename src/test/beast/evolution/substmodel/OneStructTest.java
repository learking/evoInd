package test.beast.evolution.substmodel;

import java.util.Arrays;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.OneStruct;
import beast.evolution.substitutionmodel.YN98;
import test.beast.evolution.substmodel.YN98Test.Instance;
import junit.framework.TestCase;

public class OneStructTest extends TestCase {
	public interface Instance {
	    //frequencies for pi_A, pi_T, pi_C and pi_G
        Double[] getPi();

        Double getKappa();

        Double getOmega();
        
        Double[] getCodonProb();
    }
	
    protected Instance test0 = new Instance() {
        public Double[] getPi() {
            return new Double[]{0.2, 0.3, 0.1, 0.4};
        }

        public Double getKappa() {
            return 2.0;
        }

        public Double getOmega() {
            return 2.0;
        }

        //R: a = runif(61); a = a / sum(a)
		public Double[] getCodonProb() {
			Double[] codonProbTest0 = new Double[61];
			double tmpTotal = 0.0;
			for (int i = 0; i < (codonProbTest0.length - 1); i++) {
				codonProbTest0[i] = 1.0/61.0;
				tmpTotal += codonProbTest0[i];
			}
			codonProbTest0[codonProbTest0.length - 1] = 1.0 - tmpTotal;
			return codonProbTest0;
		}
    };
    
    Instance[] all = {test0};
    
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
            RealParameter codonProb = new RealParameter(test.getCodonProb());

            Frequencies nucleoFrequencies = new Frequencies();
            nucleoFrequencies.initByName("frequencies", f, "estimate", false);

            OneStruct oneStruct = new OneStruct();
            oneStruct.initByName("kappa", test.getKappa().toString(), "omega", test.getOmega().toString(), "nucleoFrequencies", nucleoFrequencies, "codonProb", codonProb);
            
            oneStruct.prepareMatricesForTest();
            
            final double[][] rateM = oneStruct.getRateMatrix();
            assertEquals(rateM.length, 61);
            assertEquals(rateM[0].length, 61);
            
            // make sure symmetric matrix is indeed symmetric
            final double[][] symmM = oneStruct.getSymmMatrix();
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
            final double[] diagM = oneStruct.getDiagMatrix();
            System.out.print("diagMatrix:");
            System.out.println(Arrays.toString(diagM));
            assertEquals(getTotalSubstRate(rateM, diagM), 1.0, 1e-10);
            
        }
    }
}
