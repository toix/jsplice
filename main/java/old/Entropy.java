package old;

import java.util.Arrays;

import jsplice.data.Config;
import jsplice.data.RefGene;
import jsplice.io.Variants;
import jsplice.tools.Functions;

public class Entropy {
	
	/**
	 * Sequences of all variants
	 */
	private String[] sequences;
	/**
	 * True if this instance processes the sequences before the exon start
	 */
	private boolean intronExonJunction;
	/**
	 * True if the instance processes the alternative sequence (containing the mutation)
	 */
	private boolean useAlternative;
	/**
	 * Length of the sequence to be analyzed
	*/
	private int sequenceLength;
	/**
	 * {@link Variants} with corresponding {@link RefGene} and sequence
	 */
	private Variants varFile;
	/**
	 * Number of possible bases
	 */
	private int[] variantMap;
	public final int numberOfBases = Functions.bases.length();
	/**
	 * The entropy of each positions relative to the junction
	 */
	private double[] entropy;
	/**
	 * The probability of the bases at each position relative to the junction
	 */
	private double[][] probability;
	/**
	 * The entropy for every sequence at each position
	 */
	private double[][] baseEntropy;
	private double[][] baseHeight;
	private double[] positionHeight;
	
	/**
	 * @param varFile
	 *            {@link Variants} with corresponding {@link RefGene} and sequence
	 * @param isIntronExonJunction
	 *            True if this instance processes the sequences before the exon start
	 */
	public Entropy(Variants varFile, boolean intronExonJunction, boolean useAlternative) {
		this.varFile = varFile;
		this.intronExonJunction = intronExonJunction;
		this.useAlternative = useAlternative;
		this.sequenceLength = Config.getLengthModelSequence();
		variantMap = new int[varFile.size()];
		this.sequences = getSequences();
		run();
	}

	/**
	 * Extracts the sequences from the {@link Variants}
	 * @return Sequences of all variants
	 */
	private String[] getSequences() {
		System.out.println("Default base order ACGT");
		String[] sequences = new String[varFile.size(intronExonJunction)];
		int j =0;
		System.out.println("Intron-Exon: " + intronExonJunction);
		for (int i = 0; i < varFile.size(); i++) {
			boolean isIntronExonJunction = varFile.get(i).isAcceptorSite();
			if(isIntronExonJunction == this.intronExonJunction){
				sequences[j] = varFile.get(i).getSequence().toString();
				variantMap[j] = i;
				if (useAlternative) {
					int distanceToExon = varFile.get(i).getDistanceToExon();
					String alternative = varFile.get(i).getAlt();
					if(intronExonJunction){
						sequences[j] = sequences[j].subSequence(0, sequenceLength+distanceToExon) + alternative + sequences[j].substring(sequenceLength+distanceToExon+1);
					} else {
						sequences[j] = sequences[j].subSequence(0, 0+distanceToExon) + alternative + sequences[j].substring(0+distanceToExon+1);
					}
				}
				j++;
			}
		}
		return sequences;
	}

	/**
	 * Calculate probability and entropy
	 */
	public void run() {
		// parse bases to upper case
		for (int i = 0; i < sequences.length; i++) {
			sequences[i] = Functions.convertToUpperBases(sequences[i]);
		}
		// calculate probability
		probability = new double[sequenceLength][numberOfBases];
		String[] sequenceRows = Functions.transpose(sequences);
		for (int i = 0; i < sequenceLength; i++) {
			probability[i] = getProbability(sequenceRows[i]);
		}
		// calculate entropy
		entropy = new double[sequenceLength];
		for (int i = 0; i < sequenceLength; i++) {
			for (int j = 0; j < numberOfBases; j++) {
				if(probability[i][j] != 0){
					entropy[i] += - probability[i][j] * Math.log(probability[i][j]) / Math.log(2);
				}
			}
		}
		// calculate error
		double error = 1/Math.log(2) * (4-1)/(2*sequences.length);
		// calculate height
		positionHeight = new double[sequenceLength];
		baseHeight = new double[sequenceLength][numberOfBases];
		for (int i = 0; i < sequenceLength; i++) {
			positionHeight[i] = 2 - (entropy[i] + error);
			int roundedHeight = (int)Math.round(positionHeight[i]*1000);
			System.out.print(i + "\t" + roundedHeight + "\t");
			for (int j = 0; j < numberOfBases; j++) {
				baseHeight[i][j] = probability[i][j] * positionHeight[i];
				int roundedBH = (int)Math.round(baseHeight[i][j]*1000);
				System.out.print(roundedBH + " ");
			}
			System.out.println();
		}
		// calculate sequence entropy
		baseEntropy = new double[sequences.length][sequenceLength];
		double[] sequenceEntropy = new double[sequences.length];
		for (int i = 0; i < sequences.length; i++) {
			baseEntropy[i] = getSequenceEntropy(sequences[i]);
			sequenceEntropy[i] = Functions.sum(baseEntropy[i]);
		}
		Arrays.sort(sequenceEntropy);
		System.out.println("First: " + sequenceEntropy[0] + "\tLast: " + sequenceEntropy[sequenceEntropy.length-1]);
	}

	/**
	 * @param string
	 * @return
	 */
	public double[] getSequenceEntropy(String sequence) {
		double[] sequenceEntropy = new double[sequenceLength];
		for (int i = 0; i < sequenceLength; i++) {
			int baseNumber = Functions.mapNumber.get(sequence.charAt(i));
			sequenceEntropy[i] = baseHeight[i][baseNumber];
		}
		return sequenceEntropy;
	}

	/**
	 * Probability that a base occurs at a certain position
	 * @return
	 */
	public double[] getProbability(String sequence){
		double[] count = Functions.getInitializedDoubleArray(numberOfBases);
		double[] prob = new double[numberOfBases];
		for (int i = 0; i < sequence.length(); i++) {
			char base =sequence.charAt(i);
			Integer itgr =Functions.mapNumber.get(base);
			try {
				count[itgr]++;
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new IllegalArgumentException();
			}
		}
		for (int i = 0; i < numberOfBases; i++) {
			prob[i] = count[i] / sequence.length();
		}
		return prob;
	}
	
	@Override
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String out = "";
		for (int i = 0; i < sequenceLength; i++) {
			
		}
		return out;
	}
}