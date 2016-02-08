package jsplice.io;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import jsplice.data.Config;
import jsplice.data.Sequence;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.tools.Functions;
import jsplice.tools.Model;
import jsplice.tools.Result;
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

	/**
	 * Filter variants by reference sequence, alternate sequence and variantPosition range
	 * 
	 * @param varFile
	 * @param refFilter
	 *            String with RegEx or null
	 * @param altFilter
	 *            String with RegEx or null
	 * @param from
	 *            Position relative to exon junction
	 * @param to
	 *            Position relative to exon junction
	 * @return A copy of the variant file containing the filtered variants
	 */
	public static Variants filterVariants(Variants variants, String refFilter, String altFilter, int from, int to) {
		Variants newVarFile = new Variants(variants);
		if (refFilter == null || refFilter.isEmpty())
			refFilter = ".";
		if (altFilter == null || altFilter.isEmpty())
			altFilter = ".";
		for (int i = 0; i < newVarFile.size(); i++) {
			Variant variant = newVarFile.get(i);
			if (!Pattern.matches(refFilter, variant.getRef()) || !Pattern.matches(altFilter, variant.getAlt())
					|| from > variant.getDistanceToExon() || to < variant.getDistanceToExon()) {
				newVarFile.remove(i);
				i--;
			}
		}
		return newVarFile;
	}

	/**
	 * Filter variants by reference sequence, alternate sequence and variantPosition
	 * 
	 * @param varFile
	 * @param refFilter
	 *            String with RegEx or null
	 * @param altFilter
	 *            String with RegEx or null
	 * @param position
	 *            Position relative to exon junction
	 * @return A copy of the variant file containing the filtered variants
	 */
	public static Variants filterVariants(Variants variants, String refFilter, String altFilter, int position) {
		return filterVariants(variants, refFilter, altFilter, position, position);
	}

	/**
	 * Filter variants by reference sequence, alternate sequence and weather the variant is before the exon
	 * 
	 * @param varFile
	 * @param refFilter
	 *            String with RegEx or null
	 * @param altFilter
	 *            String with RegEx or null
	 * @param acceptor
	 *            Has to be true to filter variants before an exon
	 * @return A copy of the variant file containing the filtered variants
	 */
	public static Variants filterVariants(Variants variants, String refFilter, String altFilter, boolean acceptor) {
		if (acceptor) {
			return filterVariants(variants, refFilter, altFilter, -Config.getLengthTrainingSequence(), -1);
		} else {
			return filterVariants(variants, refFilter, altFilter, 1, Config.getLengthTrainingSequence());
		}
	}

	/**
	 * Filter variants by weather the variant is before or after the exon
	 * 
	 * @param acceptor
	 *            Has to be true to filter variants before an exon
	 * @return A copy of the variant file containing the filtered variants
	 */
	public static Variants filterVariants(Variants variants, boolean acceptor) {
		return filterVariants(variants, ".", ".", acceptor);
	}

	/**
	 * Removes all Variants that are related to a GT or AG junction
	 * 
	 * @return A modified copy of the VariantFile
	 */
	public static Variants filterNonACGT(Variants variants) {
		Variants filteredVariants = new Variants(variants);
		int count = 0;
		for (int i = 0; i < filteredVariants.size(); i++) {
			String sequence = filteredVariants.get(i).getSequence().getStringIntron();
			boolean isIntronExonJunction = filteredVariants.get(i).isAcceptorSite();
			if (isIntronExonJunction) {
				int last = sequence.length() - 1;
				if (!Pattern.matches("[aA][gG]", sequence.substring(last - 1, last + 1))) {
					filteredVariants.remove(i);
					i--;
					count++;
				}
			} else {
				if (!Pattern.matches("[gG][tT]", sequence.substring(0, 2))) {
					filteredVariants.remove(i);
					i--;
					count++;
				}
			}
		}
		Log.add("Number of non-GT-AG variants found:\t" + count, 3);
		return filteredVariants;
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

	// public String matrixChangeToString(){
	// // copy everything into a String
	// System.out.println();
	// System.out.println("changeMatrix");
	// String returnString = 0 + " \t" + Functions.arrayToString(weightMatrix[0]);
	// for (int l = 1; l < sequenceLength; l++) {
	// returnString = returnString + "\n" + l + " \t" + Functions.arrayToString(weightMatrixChange[l]);
	// }
	// return returnString;
	// }
	/**
	 * Remove variants that are out of analysis range <br/>
	 * @param variantsInRange
	 * @param min
	 *            relative position of the sequence
	 * @param max
	 *            relative position of the sequence
	 * @return
	 */
	public static Variants extractVariantsInRange(Variants variants, int min, int max) {
		Variants variantsInRange = new Variants(variants);
		int count = 0;
		for (int i = 0; i < variantsInRange.size(); i++) {
			int variantPosition = variantsInRange.get(i).getSequence().getPositionChangeRelative();
			// System.out.print("varPos: " + variantPosition + "...");
			if (variantPosition < min || variantPosition > max) {
				variantsInRange.remove(i);
				i--;
				count++;
				// System.out.println("removed");
			}
			// System.out.println();
		}
		Log.add("Number of variants that are out of the analysis range: \t" + count, 3);
		return variantsInRange;
	}

	// public String matrixChangeToString(){
	// // copy everything into a String
	// System.out.println();
	// System.out.println("changeMatrix");
	// String returnString = 0 + " \t" + Functions.arrayToString(weightMatrix[0]);
	// for (int l = 1; l < sequenceLength; l++) {
	// returnString = returnString + "\n" + l + " \t" + Functions.arrayToString(weightMatrixChange[l]);
	// }
	// return returnString;
	// }
	/**
	 * Remove all variants that are not within a range and in the intron <br/>
	 * 
	 * @param variantsInRange
	 * @param range
	 *            relative range from the junction position
	 * @return
	 */
	public static Variants extractVariantsInRelativeRange(Variants variants, int range) {
		Variants variantsInRange = new Variants(variants);
		int count = 0;
		for (int i = 0; i < variantsInRange.size(); i++) {
			Variant variant = variantsInRange.get(i);
			Sequence sequence = variant.getSequence();
			int variantPosition = sequence.getPositionChange();
			int junctionPosition = sequence.getPositionJunction();
			int min, max;
			if (sequence.isAcceptor()) {
				min = junctionPosition - range;
				max = junctionPosition - 1;
			} else {
				min = junctionPosition + 1;
				max = junctionPosition + range;
			}
			// System.out.print("varPos: " + variantPosition + "...");
			if (variantPosition < min || variantPosition > max) {
				variantsInRange.remove(i);
				i--;
				count++;
				// System.out.println("removed");
			}
			// System.out.println();
		}
		Log.add("Number of variants that are out of the relative analysis range: \t" + count, 3);
		return variantsInRange;
	}

	/**
	 * Remove variants that create a new cryptic splice site
	 * @param variants Variants that have to be filtered
	 * @param modelAcc Sequences containing a Matrix that is used to find cryptic acceptor sites
	 * @param modelDon Sequences containing a Matrix that is used to find cryptic donor sites
	 * @param returnCryptic true -> return all cryptic; false -> return all non cryptic
	 * @return A filtered copy of the variants
	 */
	public static Variants extractCrypticVariants(Variants variants, Model modelAcc, Model modelDon, boolean returnCryptic) {
		// System.out.println("ie? " + intronExonJunction + "\t seq len bef: " + sequenceLengthPattern);
		Variants variantsAcc = filterVariants(variants, true);
		if (variantsAcc.size() > 0) {
			variantsAcc = extractCrypticVariants(variantsAcc, modelAcc, modelDon, true, returnCryptic);
		}
		Variants variantsDon = filterVariants(variants, false);
		if (variantsDon.size() > 0) {
			variantsDon = extractCrypticVariants(variantsDon, modelAcc, modelDon, false, returnCryptic);
		}
		
		return Variants.concat(variantsAcc, variantsDon);
	}

	/**
	 * Remove variants that create a new cryptic splice site
	 * TODO individual information for sequence part next to exon (5-10 bases), not the whole model
	 * @param variantsP Variants that have to be filtered
	 * @param modelStd Sequences containing a Matrix that is used to find cryptic acceptor sites
	 * @param modelStdOtherSide Sequences containing a Matrix that is used to find cryptic donor sites
	 * @param acceptorVariants true -> all variants will be treated as acceptor variants
	 * @param returnCryptic true -> return all cryptic; false -> return all non cryptic
	 * @return A filtered copy of the variants
	 */
	public static Variants extractCrypticVariants(Variants variantsP, Model modelStd, Model modelStdOtherSide, boolean acceptorVariants, boolean returnCryptic) {
		// acceptor model is donor model?
		if (!modelStd.isAcceptor()) {
			Model modelTemp = modelStd;
			modelStd = modelStdOtherSide;
			modelStdOtherSide = modelTemp;
		}
		Variants variants = new Variants();
		// acceptor
		int minAcc = Config.getMinAnalysisPosition(acceptorVariants, true);
		int maxAcc = Config.getMaxAnalysisPosition(acceptorVariants, true);
//		System.out.println("Acc min: " + minAcc + "\t max: " + maxAcc);
//		System.out.println("seq: " + acceptorVariants + "\t model: " + true);
		Result[][] totalInformationRefAcc = new Result[variantsP.size()][maxAcc - minAcc +1];
		Result[][] totalInformationAltAcc = new Result[variantsP.size()][maxAcc - minAcc +1];
		// donor
		int minDon = Config.getMinAnalysisPosition(acceptorVariants, false);
		int maxDon = Config.getMaxAnalysisPosition(acceptorVariants, false);
//		System.out.println("Don min: " + minDon + "\t max: " + maxDon);
		Result[][] totalInformationRefDon = new Result[variantsP.size()][maxDon - minDon +1];
		Result[][] totalInformationAltDon = new Result[variantsP.size()][maxDon - minDon +1];
		
		int count = 0;
	
		System.out.println("calc cryptic Inf...");
		// calculate Individual Information for all sequences at all positions
		for (int v = 0; v < variantsP.size(); v++) {
			Variant variant = variantsP.get(v);
			Log.add(variant.toString(), 2);
			Sequence sequence = variant.getSequence();
			if (sequence.isAcceptor() != acceptorVariants) {
				throw new IllegalArgumentException("Variant is " + (sequence.isAcceptor() ? "acceptor" : "donor")
						+ " variant while parameter declares it was a(n) " + (acceptorVariants ? "acceptor" : "donor") + " variant\n"
						+ variant);
			}
			// with acceptor model
			for (int k = minAcc, j = 0; k <= maxAcc; k++, j++) {
				try {
					totalInformationRefAcc[v][j] = modelStd.getIndividualInformation(sequence, k, true, true);
					totalInformationAltAcc[v][j] = modelStd.getIndividualInformation(sequence, k, false, true);
//					System.out.print(k+": ");
//					System.out.println(totalInformationRefAcc[v][j] + " and " + totalInformationAltAcc[v][j]);
				} catch (IllegalArgumentException e) {
					throw new IllegalArgumentException(e.getMessage() + "\n" + variant);
				}
			}
//			System.out.println(Functions.arrayToString(totalInformationRefAcc[i], 2));
			// with donor model
			for (int k = minDon, j = 0; k <= maxDon; k++, j++) {
				try {
					totalInformationRefDon[v][j] = modelStdOtherSide.getIndividualInformation(sequence, k, true, true);
					totalInformationAltDon[v][j] = modelStdOtherSide.getIndividualInformation(sequence, k, false, true);
				} catch (IllegalArgumentException e) {
					throw new IllegalArgumentException(e.getMessage() + "\n" + variant);
				}
			}
			// System.out.println("nat info ref: " + individualInformationRef[i][naturalModelJunction]);
	
			int junctionPosition = sequence.getPositionJunction() - minAcc;
			boolean isCrypticAcceptorSite;
			boolean isCrypticDonorSite;
			// penalty
			if (acceptorVariants) {
				Log.add("acceptor", 2);
				isCrypticAcceptorSite = isCryptic(totalInformationRefAcc[v], totalInformationAltAcc[v], junctionPosition, false) != null;
				Log.add("donor", 2);
				isCrypticDonorSite = isCryptic(totalInformationRefDon[v], totalInformationAltDon[v], junctionPosition, true) != null;
			} else {
				Log.add("acceptor", 2);
				isCrypticAcceptorSite = isCryptic(totalInformationRefAcc[v], totalInformationAltAcc[v], junctionPosition, true) != null;
				Log.add("donor", 2);
				isCrypticDonorSite = isCryptic(totalInformationRefDon[v], totalInformationAltDon[v], junctionPosition, false) != null;
			}
			boolean isCrypticSite = isCrypticAcceptorSite || isCrypticDonorSite;
			variantsP.get(v).setCryptic(isCrypticSite);
			if (isCrypticSite == returnCryptic) {
				variants.add(variant);
			} else {
				count++;
			}
		}
		System.out.println("calculated");
		Log.add("Number of " + (acceptorVariants ? "ACCEPTOR" : "DONOR") + " variants containing " + (returnCryptic ? "NO" : "A")
				+ " cryptic splice site: " + count + "\t Number of Variants left: " + variants.size(), 3);
		return variants;
	}

	/**
	 * Defines whether the change is probably a cryptic splice site.
	 * @param totalInformationAlt
	 *            Sum of individual Information of the alternative sequence at all locations
	 * @param totalInformationRef
	 *            Sum of individual Information of the reference sequence at all locations
	 * @param naturalPos
	 *            position where the sequence is naturally spliced
	 * @param penalty true -> model and sequence belong to different site types (acceptor, donor)
	 * @return whether the sequence is probably cryptic
	 */
	public static Integer isCryptic(Result[] totalInformationRef, Result[] totalInformationAlt, int naturalPos, boolean penalty) {
//		int naturalPos = totalInformationRef[0].sequence.getPositionJunction();
		Log.add("nat pos: " + naturalPos, 1);
		int posCryptic = naturalPos;
		double informationCrypticMax = -100;
		// decrease of individual information at the natural site
		double informationNaturalRef = totalInformationRef[naturalPos].getTotalInformation();
		double informationNaturalAlt = totalInformationAlt[naturalPos].getTotalInformation();
		if (penalty) {
			informationNaturalAlt = 0;
			informationNaturalRef = 0;
		} else if (informationNaturalAlt < 0) {
			// total cryptic alt has to be at least 0
			Log.add("totalNaturalAlt log: " + informationNaturalAlt, 2) ;
			informationNaturalAlt = -Math.log10(-informationNaturalAlt+1);
			Log.add(" -> " + informationNaturalAlt, 2);
		}
		double naturalDecrease = informationNaturalRef - informationNaturalAlt;
		for (int l = 0; l < totalInformationAlt.length; l++) {
			// increase of individual information at the potential site
			double informationCrypticAlt = totalInformationAlt[l].getTotalInformation();
			if (penalty) {
				informationCrypticAlt -= 20.0 / Math.abs(l - naturalPos);
			}
			double informationCrypticRef = totalInformationRef[l].getTotalInformation();
			double crypticIncrease = informationCrypticAlt - informationCrypticRef;
			if (	informationCrypticAlt > informationNaturalAlt * 3./4
					&& informationCrypticAlt > informationNaturalRef * 1./2
					&& naturalDecrease + crypticIncrease >= 1) {
				Log.add("Position " + (l - naturalPos) + " IS  cryptic " + (penalty ? "with penalty" : "without penalty") + ": natural "
						+ informationNaturalRef + " > " + informationNaturalAlt + "\t cryptic " + informationCrypticRef + " > "
						+ informationCrypticAlt, 2);
				if (informationCrypticMax < informationCrypticAlt) {
					informationCrypticMax = informationCrypticAlt;
					posCryptic =l-naturalPos;
				}
			}
//			if (l != naturalPos && informationCrypticAlt > -1 && crypticIncrease > 0 && naturalDecrease >= 0) {
			if (true) {
				Log.add("Position " + (l - naturalPos) + " NOT cryptic " + (penalty ? "with penalty" : "without penalty") + ": natural "
						+ informationNaturalRef + " > " + informationNaturalAlt + "\t cryptic " + informationCrypticRef + " > "
						+ informationCrypticAlt, 2);
			}
		}
		return informationCrypticMax == -100 ? null : posCryptic;
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
	 * Extract variants that are probably deactivating the variant rather than activating it
	 * @param variantsP
	 * @param model
	 * @param activating true -> return activating variants; false -> return deactivating variants
	 * @return
	 */
	public static Variants filterActivatingVariants(Variants variantsP, Model model, boolean activating)  {
		variantsP = Variants.filterVariants(variantsP, model.isAcceptor());
		Variants filteredVariants = new Variants();
		int counter = 0;
		for (int i = 0; i < variantsP.size(); i++) {
			Sequence sequence = variantsP.get(i).getSequence();
			int junction = sequence.getPositionJunction();
			double ref = model.getIndividualInformation(sequence, junction, true, false).getTotalInformation();
			double alt = model.getIndividualInformation(sequence, junction, false, false).getTotalInformation();
			if (activating == ref > alt) {
				filteredVariants.add(variantsP.get(i));
			} else {
				counter ++;
			}
		}
		if (activating) {
			Log.add("Number of deactivating variants removed: " + counter + "\t Number of variants left: " + filteredVariants.size(), 3);
		} else {
			Log.add("Number of activating variants removed: " + counter + "\t Number of variants left: " + filteredVariants.size(), 3);
		}
		return filteredVariants;
	}

	/**
	 * Extract variants that are probably activating the variant rather than deactivating it
	 * @TODO write funciton
	 * @param sequences
	 * @return
	 */
	public static Variants removeStopCodon(Variants variants, boolean intronExonJunction) {
		Variants strangeVariants = new Variants();
		Model model = new Model(variants, intronExonJunction);
		for (int i = 0; i < variants.size(); i++) {
			Sequence sequence = variants.get(i).getSequence();
			int junction = sequence.getPositionJunction();
			double ref = model.getIndividualInformation(sequence, junction, true, false).getTotalInformation();
			double alt = model.getIndividualInformation(sequence, junction, false, false).getTotalInformation();
			if (ref < alt) {
				strangeVariants.add(variants.get(i));
			}
		}
		return strangeVariants;
	}

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

	/**
	 * @param variants2
	 * @return
	 */
	public static Variants deleteDuplicateJunctions(Variants variants)  {
		Variants variantsUniqueJunction = new Variants();
		HashMap<String, Variant> variantsByJunction = new HashMap<String, Variant>();
		for (int i = 0; i < variants.size(); i++) {
			String key = variants.get(i).getChromosome() + variants.get(i).getJunctionPosition();
			if (!variantsByJunction.containsKey(key)) {
				variantsByJunction.put(key, variants.get(i));
				variantsUniqueJunction.add(variants.get(i));
			}
		}
		return variantsUniqueJunction;
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
