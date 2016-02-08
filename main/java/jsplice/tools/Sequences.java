package jsplice.tools;

import java.util.ArrayList;

import jsplice.data.RefGene;
import jsplice.data.Sequence;
import jsplice.io.Variants;
import sun.security.util.Length;

/**
 * 
 * @author Tobias Gresser (gresserT@gmail.com)
 * 
 *         A set of {@link Sequence}s extracted from a {@link Variants} <br/>
 *         TODO separate Model and Sequence
 *         length
 */
public class Sequences {

	/**
	 * {@link Variants} with corresponding {@link RefGene} and sequence
	 */
	private Variants variants;
	/**
	 * Sequences of all variants
	 */
	private ArrayList<Sequence> sequences;
	/**
	 * Is true if variant is related to an exon start junction
	 */
	private boolean acceptor;
	/**
	 * Number of bases in the sequences
	 */
	private int sequenceLength;
	/**
	 * The position of the natural splice site
	 */
	private int junctionPosition;
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
	 * Constructor for standard model
	 * @param variants
	 *            {@link Variants} with corresponding {@link RefGene} and sequence
	 * @param acceptorP
	 *            True -> extract Variants before that are before exon
	 */
	public Sequences(Variants variantsP, boolean acceptorP) {
		if (variantsP.size() < 1) {
			throw new IllegalArgumentException("The parameter contains no variants.");
		}
		this.acceptor = acceptorP;
		this.variants = Filter.filterVariantType(variantsP, acceptor);
		verifySequencesEquality(variants.getSequences());
		Sequence firstSequence = variants.get(0).getSequence();
		this.sequenceLength = firstSequence.length();
		this.junctionPosition = firstSequence.getPositionJunction();
		this.sequences = variants.getSequences();
	}
	
	/**
	 * TODO static check junction position, sequence length and for equal junction type (intron-exon or exon-intron)
	 * 
	 * @param sequences
	 */
	private static void verifySequencesEquality(ArrayList<Sequence> sequences) {
		Sequence firstSequence = sequences.get(0);
		boolean intronExonJunction = firstSequence.isAcceptor();
		int sequenceLength = firstSequence.length();
		int junctionPosition = firstSequence.getPositionJunction();
		if (sequences.size() < 2) {
			throw new IllegalArgumentException("The number of sequences is to small.");
		}
		for (int i = 1; i < sequences.size(); i++) {
			Sequence sequence = sequences.get(i);
			if (intronExonJunction != sequence.isAcceptor()) {
				throw new IllegalArgumentException("The sequence must be \"" + intronExonJunction
						+ "\" for belonging to an intron exon junction: \n" + sequence.getVariant());
			}
			if (sequenceLength != sequence.length()) {
				throw new IllegalArgumentException("The sequence must have a sequence lenght of " + sequenceLength + ": \n"
						+ sequence.getVariant());
			}
			if (junctionPosition != sequence.getPositionJunction()) {
				throw new IllegalArgumentException("The Sequence must have its junction position at " + junctionPosition + ": \n"
						+ sequence.getVariant());
			}
		}

	}

	

	
	
	

	// public void runOld() {
	// // calculate probability
	// probability = Functions.getProbabilities(sequences);
	// // calculate error
	// double error = 1/Math.log(2) * (4-1)/(2*sequences.size());
	// // calculate uncertainty
	// uncertainty = new double[sequenceLength];
	// for (int l = 0; l < sequenceLength; l++) {
	// for (int b = 0; b < numberOfBases; b++) {
	// if(probability[l][b] != 0){
	// uncertainty[l] += - probability[l][b] * Math.log(probability[l][b]) / Math.log(2);
	// }
	// }
	// }
	// // calculate information content and weight matrix
	// informationContent = new double[sequenceLength];
	// weightMatrix = new double[sequenceLength][numberOfBases];
	// for (int l = 0; l < sequenceLength; l++) {
	// informationContent[l] = 2 - (uncertainty[l] + error);
	// for (int b = 0; b < numberOfBases; b++) {
	// weightMatrix[l][b] = probability[l][b] * informationContent[l];
	// }
	// }
	// // calculate individual information of the natural splice sites
	// individualInformation = new double[sequences.size()][sequenceLength];
	// double[] sequenceEntropy = new double[sequences.size()];
	// for (int i = 0; i < sequences.size(); i++) {
	// individualInformation[i] = getIndividualInformation(sequences.get(i), getJunctionPosition(), true);
	// sequenceEntropy[i] = Functions.sum(individualInformation[i]);
	// }
	// Arrays.sort(sequenceEntropy);
	// System.out.println("First: " + sequenceEntropy[0] + "\tLast: " + sequenceEntropy[sequenceEntropy.length-1]);
	//
	// // calculate correlation
	// correlation = new double[sequenceLength][sequenceLength];
	// for (int l1 = 0; l1 < sequenceLength; l1++) {
	// for (int l2 = 0; l2 < sequenceLength; l2++) {
	// correlation[l1][l2] = getCorrelation(l1, l2);
	// }
	// }
	// }

	@Override
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String out = sequences.get(0).toString();
		for (int i = 1; i < size(); i++) {
			out += "\n" + sequences.get(i).toString();
		}
		return out;
	}

	/**
	 * @return {@link jsplice.data.Sequence#junctionPosition}
	 */
	public int getJunctionPosition() {
		return junctionPosition;
	}

	public int length() {
		return sequenceLength;
	}

	/**
	 * @return Number of sequences
	 */
	public int size() {
		return sequences.size();
	}

	public Sequence get(int i) {
		return sequences.get(i);
	}

	public ArrayList<Sequence> getSequences() {
		return sequences;
	}

	public boolean isAcceptor() {
		return acceptor;
	}
}