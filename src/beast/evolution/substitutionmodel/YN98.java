package beast.evolution.substitutionmodel;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.Codon;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

public class YN98 extends GeneralSubstitutionModel {
    public Input<RealParameter> kappaInput = new Input<RealParameter>("kappa", "kappa parameter in YN98 model", Validate.REQUIRED);
    public Input<RealParameter> omegaInput = new Input<RealParameter>("omega", "kappa parameter in YN98 model", Validate.REQUIRED);
    //frequencies for pi_A, pi_T, pi_C and pi_G
    public Input<Frequencies> nucleoFreqInput =
            new Input<Frequencies>("nucleoFrequencies", "substitution model equilibrium state frequencies", Validate.REQUIRED);
    
    Frequencies nucleoFrequencies;
    
    //change frequenciesInput and ratesInput from "REQUIRED" to "OPTIONAL"
    public YN98() {
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

        nucleoFrequencies = nucleoFreqInput.get();
        updateMatrix = true;
        nrOfStates = 61;

        eigenSystem = createEigenSystem();
        rateMatrix = new double[nrOfStates][nrOfStates];
        //relativeRates, in our case, cannot be simplified into an array
        //do not need relativeRates or storedRelativeRates for now
        //relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        //storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)];

        //get parameters from inputs ... (can be implemented later)
        //rateAC = getParameter(rateACInput);
        //rateAG = getParameter(rateAGInput);
        //rateAT = getParameter(rateATInput);
        //rateCG = getParameter(rateCGInput);
        //rateCT = getParameter(rateCTInput);
        //rateGT = getParameter(rateGTInput);
    }
    
    @Override
    protected void setupRateMatrix() {
        double[] nucleoFreqs = nucleoFrequencies.getFreqs();
        
        /***************************************************/
        //need change from this point below
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
            for (int j = 0; j < i; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j];
            }
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j - 1];
            }
        }
        // bring in frequencies
        /*
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] *= fFreqs[j];
                rateMatrix[j][i] *= fFreqs[i];
            }
        }
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }
        // normalise rate matrix to one expected substitution per unit time
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            fSubst += -rateMatrix[i][i] * fFreqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
        */
    }
    
    @Override
    public boolean canHandleDataType(DataType dataType) throws Exception {
        if (dataType instanceof Codon) {
            return true;
        }
        throw new Exception("Can only handle codon data");
    }
    
}
