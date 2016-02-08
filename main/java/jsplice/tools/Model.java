/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;
import java.util.Arrays;

import jsplice.data.Config;
import jsplice.data.RefGene;
import jsplice.data.Sequence;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.io.Variants;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Model {

	/**
	 * Number of DNA Bases (=4)
	 */
	public static final int numberOfBases = Functions.bases.length();
	/**
	 * The probability of the bases at each position relative to the junction
	 */
	private double[][] probability;
	/**
	 * Weight matrix to calculate Individual Information by position and base
	 */
	private double[][] weightMatrix;
	// /**
	// * Weight matrix to calculate Individual Information by position and base of the changes
	// */
	// private double[][] weightMatrixChange;
	/**
	 * Pearson correlation of all locations between each other
	 */
	private double[][] correlation;
	/**
	 * true -> calculate weight matrix and Individual Information by position and base of the changes
	 */
	private boolean filtered = false;
	/**
	 * 
	 */
	Sequences sequences;

	/**
	 * 
	 */
	public Model() {
		// TODO Auto-generated constructor stub
	}

	//	/**
	//	 * TODO make all functions static Constructor to filter non-cryptic change sequences
	//	 * 
	//	 * @param variantsP
	//	 * @param intronExonJunction
	//	 * @param flag
	//	 */
	//	private Sequences(VariantFile variantsP, boolean intronExonJunctionP, String flag) {
	//		this.intronExonJunction = intronExonJunctionP;
	//		this.variants = variantsP.filterVariants(intronExonJunction);
	//		verifySequencesEquality(variants.getSequences());
	//		Sequence firstSequence = variants.get(0).getSequence();
	//		// this.sequenceLength = firstSequence.length();
	//		this.junctionPosition = firstSequence.getJunctionPosition();
	//		// sequence length
	//		// int sequenceLengthVariants = variants.get(0).getSequence().length();
	//		// if (sequenceLengthVariants < sequenceLength) {
	//		// throw new IllegalArgumentException("The parameter sequenceLength (" + sequenceLength
	//		// + ") has to be equal or smaller than the length of the seuquences of the variant(" + sequenceLengthVariants + ").");
	//		// }
	//		variants = VariantFile.extractNonCrypticVariants(variants, intronExonJunction);
	//		this.sequences = shortenSequences(variants);
	//		sequenceLength = this.sequences.get(0).length();
	//		this.junctionPosition = sequences.get(0).getJunctionPosition();
	//	}
	
		/**
	 * Constructor without filters
	 * @param variants
	 *            {@link Variants} with corresponding {@link RefGene} and sequence
	 * @param acceptorP
	 *            True -> extract Variants before that are before exon
	 */
	public Model(Variants variantsP, boolean acceptorP) {
		if (variantsP.size() < 1) {
			throw new IllegalArgumentException("The parameter contains no variants.");
		}
		this.sequences = new Sequences(variantsP, acceptorP);
		weightMatrix = calculateMatrix();
		System.out.println(matrixToString("matrix"));
		this.filtered = false;
	}

	/**
		 * Constructor with filters
		 * @param variants
		 * @param acceptor
		 * @param filtered true -> create the change model
		 */
		public Model(Variants variantsP, Model modelStd, Model modelStdOtherSide, boolean acceptorP) {
			if (variantsP.size() < 1) {
				throw new IllegalArgumentException("The parameter contains no variants.");
			}
			Variants variants = Variants.filterVariants(variantsP, isAcceptor());
			variants = Variants.extractCrypticVariants(variants, modelStd, modelStdOtherSide, isAcceptor(), false);
	//		variants = VariantFile.filterActivatingVariants(this.variants, modelStd, true);
			variants = Variants.filterNonACGT(variants);
			sequences = new Sequences(variants, acceptorP);
			weightMatrix = calculateMatrix();
			this.filtered = true;
		}

	/**
		 * Calculates the maximum intron length for the change matrix with sufficient variants per position
		 * 
		 * @param variants
		 * @return The intron length for the change matrix
		 */
		private int calculateChangeIntronLength(Variants variants) {
			// Extract all references of the Variants
			ArrayList<String> allRef = Variants.extractAllRef(variants);
			System.out.println(allRef);
			int intronLength = variants.get(0).getSequence().getStringIntron().length();
			int changeIntronLength = intronLength;
	//		System.out.println("for: " + (junctionPosition > 0) + " && " + (changeIntronLength == intronLength));
	//		System.out.println("junctionP: " + junctionPosition);
			if (isAcceptor()) {
				for (int i = sequences.getJunctionPosition(); i > 0 && changeIntronLength == intronLength; i--) {
					// System.out.println(allRef.get(i-1).length());
					if (allRef.get(i - 1).length() < 12) {
						changeIntronLength = sequences.getJunctionPosition() - i;
					}
				}
			} else {
				for (int i = sequences.getJunctionPosition(); i < sequences.length() - 1 && changeIntronLength == intronLength; i++) {
					if (allRef.get(i + 1).length() < 10) {
						changeIntronLength = i - sequences.getJunctionPosition();
					}
				}
			}
			System.out.println("change intron length: " + changeIntronLength);
			return changeIntronLength;
		}

	/**
	 * shorten and set sequences by the given sequence length
	 * 
	 * @param sequenceLengthNew
	 *            The desired length
	 * @param variants
	 */
	private ArrayList<Sequence> shortenSequences(Variants variants, int sequenceLengthNew) {
		ArrayList<Sequence> shorterSequences = new ArrayList<Sequence>();
		for (int i = 0; i < variants.size(); i++) {
			Variant variant = variants.get(i);
			Sequence sequence = variant.getSequence();
			int sequenceLengthBefore = sequence.length();
			int intronLengthBefore = sequence.getStringIntron().length();
			int intronLengthNew = (int) Functions.round(intronLengthBefore * sequenceLengthNew / sequenceLengthBefore, 0);
			int intronLengthDelta = Math.abs(intronLengthBefore - intronLengthNew);
			int exonLengthBefore = sequence.getStringExon().length();
			int exonLengthNew = sequenceLengthNew - intronLengthNew;
			int exonLengthDelta = Math.abs(exonLengthBefore - exonLengthNew);
			// String sequenceString;
			// // int junctionPositionAfter;
			// // int variantPositionAfter;
			// if (intronExonJunction) {
			// sequenceString = sequence.getIntronicPart().substring(intronLengthDelta) + sequence.getExonicPart().substring(0,
			// exonLengthNew);
			// // junctionPositionAfter = intronLenghtAfter;
			// // variantPositionAfter = variantPositionBefore - intronLengthDelta;
			// } else {
			// sequenceString = sequence.getExonicPart().substring(exonLengthDelta) + sequence.getIntronicPart().substring(0,
			// sequenceLengthNew);
			// // junctionPositionAfter = exonLengthAfter - 1;
			// // variantPositionAfter = variantPositionBefore - exonLengthDelta;
			// }
			try {
				Sequence newSeq = new Sequence(sequence.getStringExtended(), sequence.getLengthExonExtended(),
						sequenceLengthNew, exonLengthNew, variant);
				shorterSequences.add(newSeq);
	
				// System.out.println("before: " + intronLengthBefore + "\t after: " + intronLenghtAfter);
				// if(Math.random()>0.9){
				// System.out.println("new seq: " + newSeq);
				// System.out.println("old seq: " + sequence);
				// System.out.println(sequence.getIntronicPart());
				// System.out.println(sequence.getExonicPart().substring(0, exonLengthAfter));
				// System.out.println(sequenceString);
				// System.out.println("seq len before: " + sequenceLengthBefore);
				// System.out.println("intron length before: " + intronLengthBefore + "\t intronLength after: " + intronLenghtAfter);
				// System.out.println("intron delta: " + intronLengthDelta + "\t exon delta: " + exonLengthDelta);
				// System.out.println("exon len a: " + exonLengthAfter + "\t intron length a: " + intronLenghtAfter);
				// System.out.println("sequence length: " + sequenceString.length());
				// System.out.println();
				// }
			} catch (IllegalArgumentException e) {
				Log.add("Failed to shorten sequence to " + sequenceLengthNew + ": " + e.getMessage(), 2);
			}
		}
		return shorterSequences;
	}

	/**
		 * Calculate probability and weigth matrix
		 */
		public double[][] calculateMatrix() {
			// // calculate uncertainty
			// uncertainty = Functions.getDoubleArrayWithZeroes(sequenceLength);
			// for (int l = 0; l < sequenceLength; l++) {
			// for (int b = 0; b < numberOfBases; b++) {
			// if(probability[l][b] != 0){
			// uncertainty[l] += - Math.log(probability[l][b]) / Math.log(2);
			// }
			// }
			// }
			// informationContent = new double[sequenceLength];
			// calculate probability
			probability = Functions.getFrequencies(sequences.getSequences(), true);
	//		System.out.println(sequences);
			System.out.println(probabilityToString());
			// calculate error and weight matrix
			double[][] weightMatrix = new double[sequences.length()][numberOfBases];
			for (int l = 0; l < sequences.length(); l++) {
				// informationContent[l] = 2 - (uncertainty[l] + error);
				double error = (4.0 - 1) / (2 * Math.log(2) * sequences.size());
				for (int b = 0; b < numberOfBases; b++) {
					if (probability[l][b] == 0) {
						probability[l][b] = 1.0 / (sequences.size() * 2);
					}
					weightMatrix[l][b] = 2.0 - (-Math.log(probability[l][b]) / Math.log(2) + error);
				}
			}
			return weightMatrix;
		
			// calculate individual information of the natural splice sites
			// individualInformation = new double[sequences.size()][sequences.getSequenceLength()];
			// double[] sequenceEntropy = new double[sequences.size()];
			// for (int i = 0; i < sequences.size(); i++) {
			// individualInformation[i] = getIndividualInformation(sequences.get(i), getJunctionPosition(), true);
			// sequenceEntropy[i] = Functions.sum(individualInformation[i]);
			// }
			// Arrays.sort(sequenceEntropy);
			// System.out.println("First: " + sequenceEntropy[0] + "\tLast: " + sequenceEntropy[sequenceEntropy.length-1]);
		}

	

	// public void runOld() {
	// // calculate probability
	// probability = Functions.getProbabilities(sequences);
	// // calculate error
	// double error = 1/Math.log(2) * (4-1)/(2*sequences.size());
	// // calculate uncertainty
	// uncertainty = new double[sequences.getSequenceLength()];
	// for (int l = 0; l < sequences.getSequenceLength(); l++) {
	// for (int b = 0; b < numberOfBases; b++) {
	// if(probability[l][b] != 0){
	// uncertainty[l] += - probability[l][b] * Math.log(probability[l][b]) / Math.log(2);
	// }
	// }
	// }
	// // calculate information content and weight matrix
	// informationContent = new double[sequences.getSequenceLength()];
	// weightMatrix = new double[sequences.getSequenceLength()][numberOfBases];
	// for (int l = 0; l < sequences.getSequenceLength(); l++) {
	// informationContent[l] = 2 - (uncertainty[l] + error);
	// for (int b = 0; b < numberOfBases; b++) {
	// weightMatrix[l][b] = probability[l][b] * informationContent[l];
	// }
	// }
	// // calculate individual information of the natural splice sites
	// individualInformation = new double[sequences.size()][sequences.getSequenceLength()];
	// double[] sequenceEntropy = new double[sequences.size()];
	// for (int i = 0; i < sequences.size(); i++) {
	// individualInformation[i] = getIndividualInformation(sequences.get(i), getJunctionPosition(), true);
	// sequenceEntropy[i] = Functions.sum(individualInformation[i]);
	// }
	// Arrays.sort(sequenceEntropy);
	// System.out.println("First: " + sequenceEntropy[0] + "\tLast: " + sequenceEntropy[sequenceEntropy.length-1]);
	//
	// // calculate correlation
	// correlation = new double[sequences.getSequenceLength()][sequences.getSequenceLength()];
	// for (int l1 = 0; l1 < sequences.getSequenceLength(); l1++) {
	// for (int l2 = 0; l2 < sequences.getSequenceLength(); l2++) {
	// correlation[l1][l2] = getCorrelation(l1, l2);
	// }
	// }
	// }
	
	public void calculateCorrelation() {
		// calculate correlation
		correlation = new double[sequences.length()][sequences.length()];
		for (int l1 = 0; l1 < sequences.length(); l1++) {
			for (int l2 = 0; l2 < sequences.length(); l2++) {
				correlation[l1][l2] = getCorrelation(l1, l2);
			}
		}
	}

	/**
	 * @param location1
	 *            Location on the sequence
	 * @param location2
	 *            Location on the sequence
	 * @return Average Pearson correlation between the two locations
	 */
	public double getCorrelation(int location1, int location2) {
		// System.out.println("size: " + size());
		PearsonsCorrelation pearsons = new PearsonsCorrelation();
		// Generate a binary matrix for the two locations. Every sequence contains one 1 for the occurring base.
		double[][] binaryProbabilityL1 = new double[numberOfBases][sequences.size()];
		double[][] binaryProbabilityL2 = new double[numberOfBases][sequences.size()];
		double[][] correlationBases = new double[numberOfBases][numberOfBases];
		double correlationSum = 0.0;
		// initialize with 0
		for (int b = 0; b < numberOfBases; b++) {
			binaryProbabilityL2[b] = Functions.getInitializedDoubleArray(sequences.size());
			binaryProbabilityL1[b] = Functions.getInitializedDoubleArray(sequences.size());
		}
		for (int j = 0; j < sequences.size(); j++) {
			char base1 = sequences.get(j).charAt(location1);
			char base2 = sequences.get(j).charAt(location2);
			Integer baseIdx1 = Functions.mapNumber.get(base1);
			Integer baseIdx2 = Functions.mapNumber.get(base2);
			binaryProbabilityL1[baseIdx1][j] = 1;
			binaryProbabilityL2[baseIdx2][j] = 1;
		}
		for (int b1 = 0; b1 < numberOfBases; b1++) {
			for (int b2 = 0; b2 < numberOfBases; b2++) {
				// System.out.println(Functions.arrayToString(binaryProbabilityL1[b1]));
				// System.out.println(Functions.arrayToString(binaryProbabilityL2[b2]));
				if (Functions.sum(binaryProbabilityL1[b1]) == 0 || Functions.sum(binaryProbabilityL2[b2]) == 0) {
					correlationBases[b1][b2] = -1;
				} else if (Arrays.equals(binaryProbabilityL1[b1], binaryProbabilityL2[b2])) {
					correlationBases[b1][b2] = 1d;
				} else {
					correlationBases[b1][b2] = pearsons.correlation(binaryProbabilityL1[b1], binaryProbabilityL2[b2]);
				}
				// System.out.println(b1 + "," + b2 + ":\t " + correlationBases[b1][b2]);
			}
		}
	
		// System.out.println("l1: " + location1 + "\t l2: " + location2);
		for (int j = 0; j < sequences.size(); j++) {
			char base1 = sequences.get(j).charAt(location1);
			char base2 = sequences.get(j).charAt(location2);
			Integer baseIdx1 = Functions.mapNumber.get(base1);
			Integer baseIdx2 = Functions.mapNumber.get(base2);
			correlationSum += correlationBases[baseIdx1][baseIdx2];
			// System.out.println("j: " + j + "\t bx1: " + base1 + "\t bx2: " + base2 + "\t cor: " + correlationBases[baseIdx1][baseIdx2]);
			if (Double.isNaN(correlationBases[baseIdx1][baseIdx2])) {
			}
		}
		// if (location1 == 24 && location2 == 25) {
		// System.out.println("correlation Matrix");
		// for (int i = 0; i < correlationBases.length; i++) {
		// System.out.println(Functions.arrayToString(correlationBases[i]));
		// }
		// int j = 0;
		// char base1 = sequences.get(j).charAt(location1);
		// char base2 = sequences.get(j).charAt(location2);
		// Integer baseIdx1 = Functions.mapNumber.get(base1);
		// Integer baseIdx2 = Functions.mapNumber.get(base2);
		// System.out.println("b1: " + baseIdx1 + "\t b2: " + baseIdx2);
		// System.out.println("result: " + correlationSum/(double)size());
		// }
		return correlationSum / (double) sequences.size();
	}

	/**
		 * Calculates the individual information for a sequence with a certain junction position
		 * @param sequenceP
		 *            A {@link Sequence}
		 * @param junction
		 *            The position of the potential splice site on the pattern
		 * @param reference
		 *            false -> the alternate sequence is used
		 * @param limitLength true -> lengthIntronCryptic in Config defines the maximum intron length
		 * @return A summand for each position
		 */
		public Result getIndividualInformation(Sequence sequenceP, int junction, boolean reference, boolean limitLength) {
			Log.add("Calculate individual information for the " + (reference ? "reference" : "alternate") + " of "
					+ (sequenceP.isAcceptor() ? "an acceptor" : "a donor") + " sequence with " + (isAcceptor() ? "accptor" : "donor")
					+ " model and " + (limitLength ? "limited" : "unlimited") + " length at position " + junction + ".", 2);
			if (sequenceP.length() < sequences.length()) {
				throw new IllegalArgumentException("The length of the pattern Sequence (" + sequenceP.length()
						+ ") has to be equal or bigger than the length of this instance (" + sequences.length() + ").");
			}
			int min = Config.getMinAnalysisPosition(sequenceP.isAcceptor(), isAcceptor());
			int max = Config.getMaxAnalysisPosition(sequenceP.isAcceptor(), isAcceptor());
	//		System.out.println("II min: " + min + "\t max: " + max);
	//		System.out.println("seq: " + sequenceP.isAcceptor() + "\t model: " + isAcceptor());
			if (junction < min || junction > max) {
				throw new IllegalArgumentException("The position (=" + junction + ") has to be a number between " + min + " and " + max);
			}
	//		System.out.println("acc seq: " + sequenceP.isAcceptor() + "\t acc model: " + isAcceptor());
	//		System.out.println("The junction (=" + junction + ") is between " + min + " and " + max);
			int intronLengthMax = Config.lengthIntronCryptic;
			int patternStart = junction - sequences.getJunctionPosition();
			int matrixStart = 0;
			int matrixEnd = sequences.length();
			// limit intron length?
			if (limitLength && isAcceptor() && patternStart < junction - intronLengthMax) {
				patternStart = junction - intronLengthMax;
				matrixStart = sequences.getJunctionPosition() - intronLengthMax;
				matrixEnd = matrixStart + Config.getLengthModelExon() + intronLengthMax;
			} else if (limitLength && !isAcceptor() && matrixEnd > Config.getLengthModelExon() + intronLengthMax) {
				matrixEnd = matrixStart + Config.getLengthModelExon() + intronLengthMax;
			}
			System.out.println("junction: " + junction);
			System.out.println("getJunctionPosition(): " + sequences.getJunctionPosition());
			System.out.println("patternStart: " + patternStart);
			System.out.println("matrixStart: " + matrixStart);
			System.out.println("matrixEnd: " + matrixEnd);
			double[] individualInformation = Functions.getInitializedDoubleArray(sequences.length());
			int changePos = sequenceP.getPositionChange();
			for (int locM = matrixStart, locP = patternStart; locM < matrixEnd; locM++, locP++) {
	//			try {
				int baseNumber = Functions.mapNumber.get(sequenceP.charAt(locP));
				// is alternate sequence and change position?
				if (!reference && locP == changePos) {
					char alt = sequenceP.getAlt().charAt(0);
					baseNumber = Functions.mapNumber.get(alt);
				}
				individualInformation[locM] = weightMatrix[locM][baseNumber];
	//			} catch (StringIndexOutOfBoundsException e) {
	//				System.out.println(e.getMessage() + "\nl1:" + l1 + "\t l2: " + l2 + "\t patternStart: " + patternStart + "\t seqLen: " + sequenceLength);
	//				System.out.println("seq acc: " + sequenceP.isAcceptor() + "\t model acc: " + isAcceptor());
	//				System.out.println("exLen: " + sequenceP.lengthExtended() + "\t len: " + sequenceP.length());
	//			}
			}
			Log.add("Information: " + Functions.sum(individualInformation), 2);
			return new Result(individualInformation, junction, sequenceP);
		}

	/**
	 * Calculates the individual information for a sequence with the natural junction position
	 * 
	 * @param sequences
	 * @param reference true -> reverence sequence is used; false -> the alternate sequence is used
	 * @return sum of Individual Information for every sequence
	 */
	public double[] getIndividualInformation(Sequences sequences, boolean reference) {
		double[] indInfo = new double[sequences.size()];
		for (int s = 0; s < sequences.size(); s++) {
			Sequence sequence = sequences.get(s);
			int junction = sequence.getPositionJunction();
			indInfo[s] = getIndividualInformation(sequence, junction, reference, false).getTotalInformation();
			if (filtered && Math.abs(sequence.getPositionChangeRelative()) > 3) {
				indInfo[s] = sequence.getMaxPatternQty(AlgorithmAdministrator.quantityRelative, reference);
			}
		}
		return indInfo;
	}

	public String probabilityToString() {
		System.out.println();
		System.out.println("probability");
		String probStr = 0 + " \t" + Functions.arrayToString(probability[0], 3);
		for (int i = 1; i < probability.length; i++) {
			probStr = probStr + "\n" + i + " \t" + Functions.arrayToString(probability[i], 3);
		}
		return probStr;
	}

	public String correlationToString() {
		System.out.println();
		System.out.println("correlation");
		String corStr = 0 + " \t" + Functions.arrayToString(correlation[0], 4);
		for (int i = 1; i < correlation.length; i++) {
			corStr = corStr + "\n" + i + " \t" + Functions.arrayToString(correlation[i], 4);
		}
		return corStr;
	}

	/**
	 * @return weight matrix as String
	 */
	public String matrixToString(String title) {
		String matrixString = title + "\n" + 0 + " \t" + Functions.arrayToString(weightMatrix[0], 3);
		for (int i = 1; i < weightMatrix.length; i++) {
			matrixString = matrixString + "\n" + i + " \t" + Functions.arrayToString(weightMatrix[i], 3);
		}
		return matrixString;
	}

	/**
	 * @param junction for relative enumeration
	 * @return weight matrix as String
	 */
	public String matrixToString(String title, int junction) {
		String matrixString = title + "\n" + (0-junction) + " \t" + Functions.arrayToString(weightMatrix[0], 5);
		for (int i = 1; i < weightMatrix.length; i++) {
			if (i == junction) {
				junction--;
			}
			matrixString = matrixString + "\n" + (i-junction) + " \t" + Functions.arrayToString(weightMatrix[i], 5);
		}
		return matrixString;
	}

	public boolean isAcceptor() {
		return sequences.isAcceptor();
	}

	/**
	 * @return
	 */
	public int getJunctionPosition() {
		return sequences.getJunctionPosition();
	}
	
	public Sequences getSequences() {
		return sequences;
	}
}