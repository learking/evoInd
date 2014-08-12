package beast.evolution.substitutionmodel;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;

public class OneStruct extends GeneralSubstitutionModel {
    public Input<RealParameter> kappaInput = new Input<RealParameter>("kappa", "kappa parameter in YN98 model", Validate.REQUIRED);
    public Input<RealParameter> omegaInput = new Input<RealParameter>("omega", "kappa parameter in YN98 model", Validate.REQUIRED);
    //frequencies for pi_A, pi_C, pi_G and pi_T
    public Input<Frequencies> nucleoFreqInput =
            new Input<Frequencies>("nucleoFrequencies", "substitution model equilibrium state frequencies", Validate.REQUIRED);
    // one prob for each codon (dimension should be 61)
    // AAAAACAAGAATACAACCACGACTAGAAGCAGGAGTATAATCATGATTCAACACCAGCATCCACCCCCGCCTCGACGCCGGCGTCTACTCCTGCTTGAAGACGAGGATGCAGCCGCGGCTGGAGGCGGGGGTGTAGTCGTGGTTTACTATTCATCCTCGTCTTGCTGGTGTTTATTCTTGTTT
    public Input<RealParameter> codonProbInput = new Input<RealParameter>("codonProb", "probabilities for each codon in OneStruct model", Validate.REQUIRED);    
    double[] codonProb;
    
    Frequencies nucleoFrequencies;
    double[][] symmMatrix;
    double[] diagMatrix;
    
    //change frequenciesInput and ratesInput from "REQUIRED" to "OPTIONAL"
    public OneStruct() {
    	frequenciesInput.setRule(Validate.OPTIONAL);
        ratesInput.setRule(Validate.OPTIONAL);
        try {
        	frequenciesInput.setValue(null, this);
        	ratesInput.setValue(null, this);
        } catch (Exception e) {
        	e.printStackTrace();
			// TODO: handle exception
		}
    }
    
    @Override
    public void initAndValidate() throws Exception {
        if (frequenciesInput.get() != null) {
            throw new Exception("the frequencies attribute should not be used. Use the nucleoFrequencies instead");
        }
        if (ratesInput.get() != null) {
            throw new Exception("the rates attribute should not be used. Use the individual rates rateAC, rateCG, etc, instead.");
        }

        // set codonProb here (codonProb will not be changed after initialization)
        codonProb = new double[codonProbInput.get().getDimension()];
        for (int i = 0; i < codonProb.length; i++) {
        	codonProb[i] = codonProbInput.get().getValue(i);
        }
        // sanity check
        double totalCodonP = 0.0;
        for(int i=0;i<codonProb.length;i++)
        {
        	totalCodonP = totalCodonP + codonProb[i];
        }
        if (totalCodonP != 1.0) {
            throw new Exception("Codon probabilities should sum up to 1.");
        }
        
        nucleoFrequencies = nucleoFreqInput.get();
        updateMatrix = true;
        nrOfStates = 61;

        eigenSystem = createEigenSystem();
        
        symmMatrix = new double[nrOfStates][nrOfStates];
        diagMatrix = new double[nrOfStates];
        rateMatrix = new double[nrOfStates][nrOfStates];

    }
}
