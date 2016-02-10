/**
 * 
 */
package jsplice.tools;

import java.util.HashMap;
import java.util.regex.Pattern;

import jsplice.data.Config;
import jsplice.data.Sequence;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.io.Variants;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Filter {

	/**
	 * 
	 */
	public Filter() {
		super();
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
	public static Variants filterVariantType(Variants variants, boolean acceptor) {
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
	 * @param minRange
	 *            relative range from the junction position
	 * @return
	 */
	public static Variants extractVariantsInRelativeRange(Variants variants, int minRange, int maxRange) {
		Variants variantsInRange = new Variants();
		int count = 0;
		for (int i = 0; i < variants.size(); i++) {
			Variant variant = variants.get(i);
			Sequence sequence = variant.getSequence();
			int variantPosition = sequence.getPositionChange();
			int junctionPosition = sequence.getPositionJunction();
			if (Math.abs(variantPosition - junctionPosition) >= minRange && Math.abs(variantPosition - junctionPosition) <= maxRange) {
				variantsInRange.add(variants.get(i));
			} else {
				count++;
			}
		}
		Log.add("Number of variants that are out of the relative analysis range:\t " + count, 3);
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
		Variants variantsAcc = filterVariantType(variants, true);
		if (variantsAcc.size() > 0) {
			variantsAcc = extractCrypticVariants(variantsAcc, modelAcc, modelDon, true, returnCryptic);
		}
		Variants variantsDon = filterVariantType(variants, false);
		if (variantsDon.size() > 0) {
			variantsDon = extractCrypticVariants(variantsDon, modelAcc, modelDon, false, returnCryptic);
		}
		Log.writeToFile();
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
		
			Log.add("calc cryptic Inf...", 2);
			// calculate Individual Information for all sequences at all positions
			for (int v = 0; v < variantsP.size(); v++) {
				Variant variant = variantsP.get(v);
				Log.add(variant.toString(), 1);
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
		
				boolean isCrypticAcceptorSite;
				boolean isCrypticDonorSite;
				// penalty
				if (acceptorVariants) {
					int junctionPosition = sequence.getPositionJunction() - minAcc;
					Log.add("acceptor", 1);
					isCrypticAcceptorSite = isCryptic(totalInformationRefAcc[v], totalInformationAltAcc[v], junctionPosition, false) != null;
					Log.add("donor", 1);
					isCrypticDonorSite = isCryptic(totalInformationRefDon[v], totalInformationAltDon[v], junctionPosition, true) != null;
				} else {
					int junctionPosition = sequence.getPositionJunction() - minDon;
					Log.add("acceptor", 1);
					isCrypticAcceptorSite = isCryptic(totalInformationRefAcc[v], totalInformationAltAcc[v], junctionPosition, true) != null;
					Log.add("donor", 1);
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
			Log.add("calculated", 2);
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
		 * @param naturalJunction
		 *            position where the sequence is naturally spliced
		 * @param penalty true -> model and sequence belong to different site types (acceptor, donor)
		 * @return whether the sequence is probably cryptic
		 */
		public static Integer isCryptic(Result[] totalInformationRef, Result[] totalInformationAlt, int naturalJunction, boolean penalty) {
			Integer posCryptic = null;
			double informationCrypticMax = -100;
			// decrease of individual information at the natural site
			double informationNaturalRef = totalInformationRef[naturalJunction].getTotalInformation();
			double informationNaturalAlt = totalInformationAlt[naturalJunction].getTotalInformation();
			if (penalty) {
				informationNaturalAlt = 0;
				informationNaturalRef = 0;
			} else if (informationNaturalAlt < 0) {
				// total cryptic alt has to be at least 0
				Log.add("totalNaturalAlt log: " + informationNaturalAlt, 1) ;
				informationNaturalAlt = -Math.log10(-informationNaturalAlt+1);
				Log.add(" -> " + informationNaturalAlt, 1);
			}
			double naturalDecrease = informationNaturalRef - informationNaturalAlt;
			for (int l = 0; l < totalInformationAlt.length; l++) {
				// increase of individual information at the potential site
				double informationCrypticAlt = totalInformationAlt[l].getTotalInformation();
				if (penalty) {
					informationCrypticAlt -= 20. / Math.abs(l - naturalJunction);
				}
				double informationCrypticRef = totalInformationRef[l].getTotalInformation();
				double crypticIncrease = informationCrypticAlt - informationCrypticRef;
				if (	informationCrypticAlt > informationNaturalAlt * 3./4
						&& informationCrypticAlt > informationNaturalRef * 1./3
						&& naturalDecrease + crypticIncrease >= 1) {
					Log.add("Position " + totalInformationRef[l].position + " IS  cryptic " + (penalty ? "with penalty" : "without penalty") + ": natural "
							+ informationNaturalRef + " > " + informationNaturalAlt + "\t cryptic " + informationCrypticRef + " > "
							+ informationCrypticAlt, 1);
					if (informationCrypticMax < informationCrypticAlt) {
						informationCrypticMax = informationCrypticAlt;
						posCryptic = totalInformationRef[l].position;
					}
				}
				else if (l != naturalJunction && informationCrypticAlt > -1 && crypticIncrease > -3 && naturalDecrease >= 0) {
					Log.add("Position " + totalInformationRef[l].position + " NOT cryptic " + (penalty ? "with penalty" : "without penalty") + ": natural "
							+ informationNaturalRef + " > " + informationNaturalAlt + "\t cryptic " + informationCrypticRef + " > "
							+ informationCrypticAlt, 1);
				}
			}
			return posCryptic;
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
			variantsP = Filter.filterVariantType(variantsP, model.isAcceptor());
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

}
