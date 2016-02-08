package jsplice.io;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import jsplice.data.Sequence;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.tools.Functions;
import jsplice.tools.Sequences;

/**
 * Reads a variant file and stores the variants
 * 
 * @author Tobias Gresser (gresserT@gmail.com) TODO validate some variants TODO store Sequences and getSequences function
 */
public class Variants {

	/**
	 * All variants, the core data of the Object
	 */
	ArrayList<Variant> variants = new ArrayList<Variant>();
	/**
	 * The IDs of the variants that are not stored in the variants field because some part of the variant isn't valid
	 */
	ArrayList<Integer> missingVariantIDs = new ArrayList<Integer>();
	/**
	 * @author Tobias Gresser (gresserT@gmail.com) The class handles different variant files
	 */
	// /**
	// * Row names
	// */
	// String[] header;
	// /**
	// * HashMap to get the number by the row name
	// */
	// HashMap<String, Integer> annotationRowNumber = new HashMap<String, Integer>();
	// /**
	// * HashMap to get the name by the row number
	// */
	// HashMap<Integer, String> annotationRowTitle = new HashMap<Integer, String>();
	/**
	 * HashMap that stores the Variants by start position, chromosome and alternate
	 */
	HashMap<Integer, HashMap<String, HashMap<String, Variant>>> variantHash = new HashMap<Integer, HashMap<String, HashMap<String, Variant>>>();
	/**
	 * The number of variants that are related to an exon start junction
	 */
	int numberOfVariantsBeforeExon = 0;
	/**
	 * Filter pattern that is applied on the DNA change pattern to find splice variants
	 */
	final static Pattern relevantLines = Pattern.compile(".*[0-9][+-][0-9].*");

	/**
	 * @param variantFileName
	 *            Name and path of the file
	 */
	public Variants() {
	}

	/**
	 * DON'T USE
	 * 
	 * @param variantFileName
	 *            Name and path of the file
	 */
	Variants(String variantFileName) {
	}

	/**
	 * Make a soft copy of the VariantFile (no Variants are duplicated)
	 * 
	 * @param variantFileName
	 *            Name and path of the file
	 * @throws UnexpectedException 
	 */
	public Variants(Variants varFile) {
		for (int i = 0; i < varFile.size(); i++) {
			variants.add(varFile.get(i));
			addToHash(variants.get(i));
		}
		countPreExonVariants();
	}

	void initial() {
		for (int i = 0; i < variants.size(); i++) {
			addToHash(variants.get(i));
		}
		countPreExonVariants();
	}

	/**
	 * Count the number of variants that are related to an exon start junction
	 */
	void countPreExonVariants() {
		for (int i = 0; i < variants.size(); i++) {
			if (variants.get(i).isAcceptorSite()) {
				numberOfVariantsBeforeExon++;
			}
		}
	}

	/**
	 * Creates a HashMap that stores the Variants by start position
	 * 
	 * @param variant
	 * @return false - if the hash already contained a variant with the same parameters
	 * @throws UnexpectedException 
	 */
	private boolean addToHash(Variant variant) {
		// if(variant.getStart() == 31816338 && variant.getChromosome().equals("11") && variant.getAlt().equals("G"))
		// System.out.println("adding matching variant:\n" + variant);
		// if (variant.getGene().equals("pax6tv2")) {
		// System.out.println("adding variant to hash....\n" + variant);
		// }
		int start = variant.getStart();
		String chrom = variant.getChromosome();
		String alt = variant.getAlt();
		// Start
		if (variantHash.containsKey(start)) {
			HashMap<String, HashMap<String, Variant>> hashByChr = variantHash.get(start);
			// Chromosome
			if (hashByChr.containsKey(chrom)) {
				HashMap<String, Variant> hashByAlt = hashByChr.get(chrom);
				// Alternative
				if (hashByAlt.containsKey(alt)) {
					Variant oldVariant = hashByAlt.get(alt);
					// Replace the old variant if the new variant contains more information
					if (oldVariant.getRefSeq() != variant.getRefSeq() && variant.getRefSeqAccession() != null
							&& oldVariant.getRefSeqAccession() == null) {
						boolean replaced = hashByAlt.put(alt, variant) == oldVariant;
						int idx = variants.indexOf(oldVariant);
						variants.set(idx, variant);
					}
					Log.add("There are two variants that describe the same sequence change:\nold:\t" + oldVariant + "\nnew:\t" + variant, 1);
					return false;
				} else {
					hashByAlt.put(alt, variant);
					return true;
				}
			} else {
				HashMap<String, Variant> hashByAlt = new HashMap<String, Variant>();
				hashByAlt.put(alt, variant);
				hashByChr.put(chrom, hashByAlt);
				return true;
			}
		} else {
			HashMap<String, Variant> hashByAlt = new HashMap<String, Variant>();
			HashMap<String, HashMap<String, Variant>> hashByChr = new HashMap<String, HashMap<String, Variant>>();
			hashByAlt.put(alt, variant);
			hashByChr.put(chrom, hashByAlt);
			variantHash.put(start, hashByChr);
			return true;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String variantString = get(0).toString();
		for (int i = 1; i < size(); i++) {
			variantString += "\n" + get(i).toString();
		}
		return variantString;
	}

	public int size() {
		return variants.size();
	}

	/**
	 * Returns the number of variants before or after the exon
	 * 
	 * @param intronExonJunction
	 *            True -> Number of variants before an exon
	 * @return Number of variants
	 */
	public int size(boolean intronExonJunction) {
		if (intronExonJunction) {
			return numberOfVariantsBeforeExon;
		} else {
			return variants.size() - numberOfVariantsBeforeExon;
		}
	}

	public ArrayList<Variant> getVariants() {
		return variants;
	}

	public String[] getRefSeqAccessions() {
		String[] refSeqAcessions = new String[variants.size()];
		for (int i = 0; i < size(); i++) {
			refSeqAcessions[i] = variants.get(i).getRefSeqAccession();
		}
		return refSeqAcessions;
	}

	// void setHashHeader(String headerLine) {
	// header = headerLine.split("\t");
	// // create hash map for annotation titles and row numbers
	// for (int i = 0; i < header.length; i++) {
	// annotationRowNumber.put(header[i], i);
	// annotationRowTitle.put(i, header[i]);
	// }
	// }

	public ArrayList<Integer> getMissingVariantIDs() {
		return missingVariantIDs;
	}

	/**
	 * Returns the variant specified by the parameter
	 * 
	 * @param index
	 *            of the desired variant
	 * @return The specified Variant
	 */
	public Variant get(int i) {
		return variants.get(i);
	}

	/**
	 * Remove a variant from this VariantFile
	 * 
	 * @param i
	 *            index of the variant
	 * @return The removed variant
	 */
	public Variant remove(int i) {
		Variant variant = variants.get(i);
		if (variant.isAcceptorSite()) {
			numberOfVariantsBeforeExon--;
		}
		// if(variant.getStart() == 31816338 && variant.getChromosome().equals("11") && variant.getAlt().equals("G"))
		// System.out.println("removing matching variant:\n" + variant);
		// if (variant.getGene().equals("pax6tv2")) {
		// System.out.println("removing....\n" + variant);
		// }
		boolean removed = variantHash.get(variant.getStart()).get(variant.getChromosome()).remove(variant.getAlt()) == variant;
		if (!removed) {
			Log.add("The variant was not found in the variant hash and could not be removed:\n" + variant, 2);
		}
		return variants.remove(i);
	}

	/**
	 * @param variant
	 */
	public boolean add(Variant variant) {
		// if (variant.getGene().contains("PAX6") || variant.getGene().contains("pax6")) {
		// System.out.println("adding....\n" + variant);
		// }
		// variant already stored?
		if (addToHash(variant)) {
			variants.add(variant);
			if (variant.isAcceptorSite()) {
				numberOfVariantsBeforeExon++;
			}
			return true;
		} else {
			Variant oldVariant = variantHash.get(variant.getStart()).get(variant.getChromosome()).get(variant.getAlt());
			if (!oldVariant.getGene().equals(variant.getGene())) {
				Log.add("Variants with start at " + variant.getStart() + " are equal." + "\n old variant: " + oldVariant.getChromosome()
						+ " " + oldVariant.getGene() + "\n new variant: " + variant.getChromosome() + " " + variant.getGene(), 2);
			}
			return false;
		}
	}

	/**
	 * @return All Sequences
	 */
	public ArrayList<Sequence> getSequences() {
		ArrayList<Sequence> sequences = new ArrayList<Sequence>();
		for (int i = 0; i < variants.size(); i++) {
			sequences.add(variants.get(i).getSequence());
		}
		return sequences;
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
	
		// /**
		// * Extract all references of the Variants.
		// * @param variant true-> take variants field; false -> take sequences field
		// * @return All references of the Variants
		// */
		// private ArrayList<String> extractAllRef(boolean variant) {
		// ArrayList<String> allRef = new ArrayList<String>();
		// ArrayList<String> allAlt = new ArrayList<String>();
		// // initialize list with empty strings
		// int sequenceLength;
		// int size;
		// if (variant) {
		// sequenceLength = variants.get(0).getSequence().length();
		// size = variants.size();
		// } else{
		// sequenceLength = get(0).length();
		// size = sequences.size();
		// }
		// for (int l = 0; l < sequenceLength; l++) {
		// allRef.add("");
		// allAlt.add("");
		// }
		// // add all reference/alternate characters
		// for (int i = 0; i < size; i++) {
		// Sequence sequence;
		// int position;
		// if (variant) {
		// sequence = variants.get(i).getSequence();
		// position = sequence.getChangePosition();
		// } else {
		// sequence = sequences.get(i);
		// position = sequence.getChangePosition();
		// }
		// String refBase;
		// String altBase;
		// refBase = sequence.getRef();
		// altBase = sequence.getAlt();
		// String refReplace = allRef.get(position) + refBase;
		// String altReplace = allAlt.get(position) + altBase;
		// allRef.set(position, refReplace);
		// allAlt.set(position, altReplace);
		// }
		// Log.add("Sequence length for change matrix is now: " + sequenceLength);
		// return allRef;
		// }
		/**
		 * Extract all references of the Variants.
		 * 
		 * @param variant
		 *            true-> take variants field; false -> take sequences field
		 * @return All references of the Variants
		 */
		public static ArrayList<String> extractAllRef(Variants variants) {
			ArrayList<String> allRef = new ArrayList<String>();
			ArrayList<String> allAlt = new ArrayList<String>();
			// initialize list with empty strings
			int sequenceLength;
			int size;
			sequenceLength = variants.get(0).getSequence().length();
			size = variants.size();
			for (int l = 0; l < sequenceLength; l++) {
				allRef.add("");
				allAlt.add("");
			}
			// add all reference/alternate characters
			for (int i = 0; i < size; i++) {
				Sequence sequence;
				sequence = variants.get(i).getSequence();
				int position = sequence.getPositionChange();
				if (position >= 0 && position < sequenceLength) {
					String refBase;
					String altBase;
					refBase = sequence.getRef();
					altBase = sequence.getAlt();
					String refReplace = allRef.get(position) + refBase;
					String altReplace = allAlt.get(position) + altBase;
					allRef.set(position, refReplace);
					allAlt.set(position, altReplace);
				}
			}
			Log.add("Sequence length for change matrix is now: " + sequenceLength);
			return allRef;
		}

	/**
	 * 
	 * @param
	 * @return
	 */
	public static Variants concat(Variants variants1, Variants variants2) {
		String log = "Concatenation: " + variants1.size() + " variants + " + variants2.size() + " variants";
		Variants variants = new Variants(variants2);
		for (int i = 0; i < variants1.size(); i++) {
			boolean added = variants.add(variants1.get(i));
		}
		Log.add(log + " = " + variants.size() + " variants", 3);
		return variants;
	}

	

//	/**
//	 * Defines whether the change is probably a cryptic splice site.
//	 * 
//	 * @param totalInformationAlt
//	 *            Sum of individual Information of the alternative sequence at all locations
//	 * @param totalInformationRef
//	 *            Sum of individual Information of the reference sequence at all locations
//	 * @param naturalPos
//	 *            position where the sequence is naturally spliced
//	 * @param penalty true -> model and sequence belong to different sites (acceptor, donor)
//	 * @return whether the position mutation
//	 */
//	public static boolean isCryptic(double[] totalInformationRef, double[] totalInformationAlt, int naturalPos, boolean penalty) {
//		if (!penalty) {
//			return isCryptic(totalInformationRef, totalInformationAlt, naturalPos);
//		} else {
//			System.out.println("nat pos: " + naturalPos);
//			// decrease of individual information at the natural splice site
//			double totalNaturalRef = 0;
//			double totalNaturalAlt = 0;
//			double naturalDecrease = totalNaturalRef - totalNaturalAlt;
//			for (int i = 0; i < totalInformationAlt.length; i++) {
//				// increase of individual information at the potential site
//				double totalCrypticAlt = totalInformationAlt[i] - 20.0 / Math.abs(i - naturalPos);
//				double totalCrypticRef = totalInformationRef[i];
//				double crypticIncrease = totalCrypticAlt - totalCrypticRef;
//				if (	totalCrypticAlt > totalNaturalAlt * 3/4
//						&& totalCrypticAlt > totalNaturalRef * 1/2
//						&& naturalDecrease + crypticIncrease >= 1) {
//					System.out.println((i-naturalPos) + " IS cry: " + totalCrypticAlt + " > " + totalNaturalAlt + " * 3/4"
//							+ "\t && " + totalCrypticAlt + " > " + totalNaturalRef + " * 1/2"
//							+ "\t && " + naturalDecrease + " + " + crypticIncrease + " >= 1");
//					return true;
//				}
//				if (crypticIncrease > 0 || totalCrypticAlt > 0) {
//					System.out.println((i-naturalPos) + " not cry: " + totalCrypticAlt + " > " + totalNaturalAlt + " * 3/4"
//							+ "\t && " + totalCrypticAlt + " > " + totalNaturalRef + " * 1/2"
//							+ "\t && " + naturalDecrease + " + " + crypticIncrease + " >= 1");
//				}
//			}
//			return false;
//		}
//	}

	/**
	 * @param variantsPathogene
	 * @return
	 */
	public static Variants changeRandomBase(Variants variantsP) {
		Variants variants = new Variants();
		for (int i = 0; i < variantsP.size(); i++) {
			Variant variant = variantsP.get(i);
			Variant variantCopy = new Variant(variant);
			variantCopy.changeRandomBase();
//			System.out.println(variant.getAlt() + " > " + variantCopy.getAlt());
			variants.add(variantCopy);
		}
		return variants;
	}

	/*
	 * public void writeFile(String fileName) { File file = new File(fileName); FileWriter fWriter = null; try { // delete old file and
	 * create file if (file.exists()) file.delete(); if (file.createNewFile()) { file.setReadable(true); fWriter = new FileWriter(file);
	 * String title = new String(header[0]); // write title as first line for (int i = 0; i < header.length; i++) title =
	 * title.concat("\t").concat(header[i]); fWriter.write(title); fWriter.append("%n"); // write variants lines to file for
	 * (Iterator<Variant> iterator = variants.iterator(); iterator .hasNext();) { fWriter.write(iterator.next().getString());
	 * fWriter.append("%n"); } } else System.out.println("Unable to create file."); } catch (IOException e) { e.printStackTrace(); } finally
	 * { if (fWriter != null) try { fWriter.close(); } catch (IOException e) { e.printStackTrace(); } } //
	 * System.out.println(variants.get(14).getAnnotationByName("Type") .getString()); }
	 */
}
