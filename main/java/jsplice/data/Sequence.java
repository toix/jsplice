/**
 * 
 */
package jsplice.data;

import java.util.HashMap;

import jsplice.tools.Cluster;
import jsplice.tools.Functions;
import jsplice.tools.Result;

/**
 * TODO extend the sequence internally to the left and the right without changes to the methods
 * @author Tobias Gresser (gresserT@gmail.com)
 * A sequence String corresponding to a splice variant
 */
public class Sequence {

	/**
	 * The main field of the Object
	 */
	private String sequenceString;
	/**
	 * The main field of the Object
	 */
	private String intronicExtensionSequence;
	/**
	 * The main field of the Object
	 */
	private String exonicExtensionSequence;
	/**
	 * Position of the variant in the reference sequence
	 */
	private int changePosition;
	/**
	 * Position of the exon splice site junction relative to the sequence start
	 */
	private int junctionPosition;
	/**
	 * The variant the sequence belongs to
	 */
	private Variant variant;
	
///**
//		 * Create a {@link Sequence} with the fields of the variant and a specific string length
//		 * @param sequenceStr The Bases of the Sequence as String
//		 * @param desiredLength of the sequence. Has to be at most 1/3 of the length of sequenceStr
//		 * @param exonLength Length of the exonic part of the sequence
//		 * @param variantP The variant that the sequence belongs to
//		 * @return
//		 */
//		public Sequence(String sequenceStr, int desiredLength, int exonLength, Variant variantP) throws IllegalArgumentException{
//			// check sequence validity
//			if(!Functions.isValidDNA(sequenceStr)){
//				throw new IllegalArgumentException("The sequence is not a valid:\n" + sequenceStr);
//			}
//			// check length validity
//			if (desiredLength*3 > sequenceStr.length()){
//				throw new IllegalArgumentException("The length has to be at most 1/3 of the lengt of the sequence: " + desiredLength + " / " + sequenceStr.length());
//			}
//			int extensionLengths = (sequenceStr.length() - length()) / 2;
//			// set junctionPosition
//			if(variantP.isIntronExonJunction()){
//				this.junctionPosition = sequenceStr.length() - exonLength;
//			} else {
//				this.junctionPosition = exonLength - 1;
//			}
//			this.changePosition = junctionPosition + variantP.getDistanceToExon();
//			this.sequenceString = Functions.convertToUpperBases(sequenceStr);
//			this.variant = variantP;
//			
//			if (junctionPosition < 0 || junctionPosition >= length()) {
//				throw new IllegalArgumentException("Junction position at " + junctionPosition + " outside of sequence String range.");
//			}
//			if (changePosition < 0 || changePosition >= length()) {
//	//			System.out.println("Variant position out of range: " + changePosition + "\t exonLen: " + exonLength);
//				throw new IllegalArgumentException("Variant position at" + changePosition + " outside of sequence String range.");
//			}
//		}

//	/**
//	 * @param sequence 
//	 * @param junctionPosition 
//	 * 
//	 */
//	public Sequence(String sequence, int junctionPosition, int variantPosition, Variant variant, int i) {
//		sequenceString = Functions.convertToUpperBases(sequence);
//		this.junctionPosition = junctionPosition;
//		this.changePosition = variantPosition;
//		this.variant = variant;
//	}

	/**
	 * Create a {@link Sequence} with the fields of the variant and a specific string length
	 * @param sequenceP The Bases of the Sequence as String
	 * @param extendedExonLength exon length of the sequence parameter
	 * @param length desired sequence length; must not be longer than the sequence parameter
	 * @param exonLength Length of the exonic part of the sequence
	 * @param variantP The variant that the sequence belongs to
	 * @return
	 */
	public Sequence(String sequenceP, int extendedExonLength, int length, int exonLength, Variant variantP) throws RuntimeException{
		// check sequence validity
		if(!Functions.isValidDNA(sequenceP)){
			throw new IllegalArgumentException("The sequence is not a valid:\n" + sequenceP);
		}
		// check length validity
		if (length > sequenceP.length()){
			throw new IllegalArgumentException("The desired length has to be shorter than the lengt of the sequence: " + length + " / " + sequenceP.length());
		}
		this.variant = variantP;
		// length of the extending sequences
		int exonicExtensionLength = extendedExonLength - exonLength;
		int intronicExtensionLength = (sequenceP.length() - extendedExonLength) - (length - exonLength);
		// sequence Strings
		sequenceP = Functions.convertToUpperBases(sequenceP);
		if (isAcceptor()) {
			intronicExtensionSequence = sequenceP.substring(0, intronicExtensionLength);
			this.sequenceString = sequenceP.substring(intronicExtensionLength, intronicExtensionLength + length);
			exonicExtensionSequence = sequenceP.substring(intronicExtensionLength + length, sequenceP.length());
		} else {
			exonicExtensionSequence = sequenceP.substring(0, exonicExtensionLength);
			this.sequenceString = sequenceP.substring(exonicExtensionLength, exonicExtensionLength + length);
			intronicExtensionSequence = sequenceP.substring(exonicExtensionLength + length, sequenceP.length());
		}
		// set junctionPosition
		if(isAcceptor()){
			this.junctionPosition = length - exonLength;
		} else {
			this.junctionPosition = exonLength - 1;
		}
		// change position
		this.changePosition = junctionPosition + variantP.getDistanceToExon();
		// check junction position validity
		if (junctionPosition < 0 || junctionPosition > length()) {
//			System.out.println("\t junPos: " + junctionPosition + "\t dist: " + variantP.getDistanceToExon() + "\t exLen: " + exonLength);
			throw new IllegalArgumentException("Junction position at " + junctionPosition + " outside of sequence String range.");
		}
//		// check change position validity
//		if (changePosition < 0 || changePosition >= length()) {
//			System.out.println("Variant position out of range: " + changePosition + "\t exonLen: " + exonLength + "\t seqLen: " + length());
//			throw new VariantToFarAwayException("Variant position at " + changePosition + " outside of sequence String range.");
//		}
	}
	
	@Override
	public String toString() {
		if (isAcceptor()) {
			return getStringIntron() + "_" + getStringExon() + " junPos: " + junctionPosition + " varPos: " + getPositionChangeRelative()
					+ " alt: " + getAlt() + " intronExon: " + isAcceptor();
		} else {
			return getStringExon() + "_" + getStringIntron() + " junPos: " + junctionPosition + " varPos: +" + getPositionChangeRelative()
			+ " alt: " + getAlt() + " intronExon: " + isAcceptor();
		}
	}
	
	/**
	 * TODO many scores are 0 instead of -1
	 * @param quantityRelative
	 * @param reference
	 * @return
	 */
	public double getMaxPatternQty(HashMap<String,Cluster> quantityRelative, boolean reference) {
		double deltaPattern = 0;//11.5;
		int lengthPatternMax = Config.lengthIntronPatternMax;
		int distanceJunctionMin = Config.getDistanceClusterMin(isAcceptor());
		int distanceJunctionMax = Config.getDistanceClusterMax();
		int posChangeRel = Math.abs(getPositionChangeRelative());
		int posChangeAbs = getPositionChange();
		double quantityMax = - 1;
		// find max pattern for ref and alt
		if (posChangeRel >= distanceJunctionMin && posChangeRel <= distanceJunctionMax) {
			for (int length = 1; length <= lengthPatternMax; length++) {
				for (int shift = 0; shift < length; shift++) {
					int from = posChangeAbs + shift - length + 1;
					int to = posChangeAbs + shift;
					int junction = getPositionJunction();
	//				System.out.println("Math.abs(" + from + " - " + junction + ") > 3	&& Math.abs(" + to + " - " + junction + ") > 3)");
					if (Math.abs(from - junction) > distanceJunctionMin	&& Math.abs(to - junction) > distanceJunctionMin) {
						String pattern = substring(from, to + 1, reference);
						double quantity = 0;
						if (quantityRelative.containsKey(pattern)) {
	//						System.out.println("new max: " + pattern);
							quantity = quantityRelative.get(pattern).getInformation(pattern);
						} else {
	//						 System.out.println("No quantity found for pattern " + pattern);
						}
						if (quantityMax < quantity) {
							quantityMax = quantity;
						}
					}
				}
			}
			quantityMax -= deltaPattern;
		} else {
			quantityMax = 0;
		}
		return quantityMax;
	}
	
	/**
	 * Rate the sequence without knowledge about the variant
	 * TODO write function
	 * @param quantityRelative
	 * @param reference
	 * @return
	 */
	public double getMaxPatternQtyNoVariant(HashMap<String,Cluster> quantityRelative, boolean reference) {
		double deltaPattern = 11.5;
		int lengthPatternMax = Config.lengthIntronPatternMax;
		int distanceJunctionMin = Config.getDistanceClusterMin(isAcceptor());
		int distanceJunctionMax = Config.getDistanceClusterMax();
		int posChangeRel = Config.getDistanceClusterMax();//Math.abs(getPositionChangeRelative());
		int posChangeAbs = getPositionChange();
		double quantityMax = - 1;
		// find max pattern for ref and alt
		if (posChangeRel >= distanceJunctionMin && posChangeRel <= distanceJunctionMax) {
			for (int length = 1; length <= lengthPatternMax; length++) {
				for (int shift = 0; shift < length; shift++) {
					int from = posChangeAbs + shift - length + 1;
					int to = posChangeAbs + shift;
					int junction = getPositionJunction();
	//				System.out.println("Math.abs(" + from + " - " + junction + ") > 3	&& Math.abs(" + to + " - " + junction + ") > 3)");
					if (Math.abs(from - junction) > distanceJunctionMin	&& Math.abs(to - junction) > distanceJunctionMin) {
						String pattern = substring(from, to + 1, reference);
						double quantity = 0;
						if (quantityRelative.containsKey(pattern)) {
	//						System.out.println("new max: " + pattern);
							quantity = quantityRelative.get(pattern).getInformation(pattern);
						} else {
	//						 System.out.println("No quantity found for pattern " + pattern);
						}
						if (quantityMax < quantity) {
							quantityMax = quantity;
						}
					}
				}
			}
			quantityMax -= deltaPattern;
		} else {
			quantityMax = 0;
		}
		return quantityMax;
	}

	/**
	 * @return length of the first part of the sequence (intron or exon)
	 */
	public int getLengthFirst(){
		if (isAcceptor()) {
			return getPositionJunction();
		}
		else {
			return getPositionJunction()+1;
		}
	}
	
	/**
	 * @return length of the first part of the extended sequence (intron or exon)
	 */
	public int getLengthFirstExtended(){
		return getLengthPreExtension() + getLengthFirst();
	}
	
	/**
	 * @return length of the sequence extending string before the sequence string
	 */
	public int getLengthPreExtension(){
		if (isAcceptor()) {
			return intronicExtensionSequence.length();
		}
		else {
			return exonicExtensionSequence.length();
		}
	}

	public int getPositionChange() {
		return changePosition;
	}
	
	public int getPositionChangeRelative() {
		return changePosition - junctionPosition;
//		return variant.getDistanceToExon();
	}
	
	public int getPositionChangeExtended() {
		return changePosition + getLengthPreExtension();
	}

	/**
	 * @return {@link #junctionPosition}
	 */
	public int getPositionJunction() {
		return junctionPosition;
	}
	public int getPositionJunctionExtended() {
		return junctionPosition + getLengthPreExtension();
	}
	
	public String getReference(){
		return sequenceString.substring(changePosition, changePosition+1);
	}
	
	public String getStringIntron(){
		if(variant.isAcceptorSite()){
			return sequenceString.substring(0, junctionPosition);
		} else {
			return sequenceString.substring(junctionPosition+1);
		}
	}
	public String getStringExon(){
		if(variant.isAcceptorSite()){
			return sequenceString.substring(junctionPosition);
		} else {
			return sequenceString.substring(0, junctionPosition+1);
		}
	}
	
	/**
	 * @return The base wise sequence as a String
	 */
	public String getStringReference() {
		return sequenceString;
	}
	/**
	 * @return The base wise sequence as a String with extended ends
	 */
	public String getStringExtended() {
		if (isAcceptor()) {
			return intronicExtensionSequence + sequenceString + exonicExtensionSequence;
		}else {
			return exonicExtensionSequence + sequenceString + intronicExtensionSequence;
		}
	}

	public String getStringAlternate(){
		String alternate = sequenceString.substring(0, changePosition) + getAlt() + sequenceString.substring(changePosition+1);
		return alternate;
	}

	/**
	 * @param l location
	 * @return The base at a certain Location of the sequence
	 */
	public char charAt(int l) {
		return stringAt(l).charAt(0);
	}

	/**
	 * @param l location
	 * @return The base at a certain Location of the sequence
	 */
	public String stringAt(int l) {
		return substring(l, l+1);
	}

	/**
	 * @param begin sequence index; inclusive
	 * @param end sequence index; exclusive
	 * @return parameterized part of the whole sequence
	 * @throws IndexOutOfBoundsException if the beginIndex is negative, or endIndex is larger than the length of this String object, or beginIndex is larger than endIndex
	 */
	public String substring(int begin, int end) throws IndexOutOfBoundsException {
		return getStringExtended().substring(begin + getLengthPreExtension(), end + getLengthPreExtension());
	}

	/**
	 * @param begin
	 * @param end
	 * @param reference
	 * @return
	 */
	public String substring(int begin, int end, boolean reference) {
		if (!reference) {
//			System.out.println(getStringExtended());
//			System.out.println(getPositionChangeExtended());
//			System.out.println(begin + getLengthPreExtension());
//			System.out.println(end + getLengthPreExtension());
			String alternateExtended = getStringExtended().substring(0, getPositionChangeExtended()) + getAlt() + getStringExtended().substring(getPositionChangeExtended() + 1);
			return alternateExtended.substring(begin + getLengthPreExtension(), end + getLengthPreExtension());
		} else {
			return substring(begin, end);
		}
	}

	public int length() {
		return sequenceString.length();
	}
	
	/**
	 * @return
	 */
	public int lengthExtended() {
		return intronicExtensionSequence.length() + sequenceString.length() + exonicExtensionSequence.length();
	}

	/**
	 * @return
	 */
	public int getLengthExon() {
		return getStringExon().length();
	}

	/**
	 * @return
	 */
	public int getLengthExonExtended() {
		return getStringExon().length() + exonicExtensionSequence.length();
	}

	/**
	 * @return
	 */
	public int getLengthIntron() {
		return length() - getLengthExon();
	}

	/**
	 * @return
	 */
	public int getLengthIntronExtended() {
		return lengthExtended() - getLengthExonExtended();
	}

	public boolean isAcceptor(){
		return variant.isAcceptorSite();
	}
	
	/**
	 * @return The reference char of the variant as String
	 */
	public String getRef(){
		return variant.getRef();
	}

	/**
	 * @return The alternate char of the variant as String
	 */
	public String getAlt() {
		return variant.getAlt();
	}

	public Variant getVariant() {
		return variant;
	}
}
