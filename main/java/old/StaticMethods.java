/**
 * 
 */
package old;

import java.rmi.UnexpectedException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import jsplice.data.Sequence;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.io.Variants;
import jsplice.tools.AlgorithmAdministrator;
import jsplice.tools.Filter;
import jsplice.tools.Functions;
import jsplice.tools.Sequences;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class StaticMethods {

	/**
	 * 
	 */
	public StaticMethods() {
		// TODO Auto-generated constructor stub
	}

	/**
	 * Extract variants that are probably deactivating the variant rather than activating it
	 * @param variants
	 * @param intronExonJunction2
	 * @return
	 * @throws UnexpectedException 
	 */
//	public static Variants extractDeactivatingVariants(Variants variants, Variants variantsBenign, boolean intronExonJunction) throws UnexpectedException {
//		Variants deactivatingVariants = new Variants();
//		Sequences model = new Sequences(variants, intronExonJunction);
//		double[][] deactivatingIndicator = calculateDeactivatingDifference(variants, variantsBenign);
//		int counter = 0;
//		for (int i = 0; i < variants.size(); i++) {
//			Sequence sequence = variants.get(i).getSequence();
//			int changePosition = sequence.getPositionChange();
//			// is the sequence change in relevant range?
//			if (changePosition > 0 && changePosition < deactivatingIndicator.length) {
//				int junction = sequence.getPositionJunction();
//				double Rref = Functions.sum(model.getIndividualInformation(sequence, junction, true, false));
//				double Ralt = Functions.sum(model.getIndividualInformation(sequence, junction, false, false));
//				
//				int changeIndex = Functions.mapNumber.get(sequence.getRef().charAt(0)) + Functions.mapNumber.get(sequence.getAlt().charAt(0)) * 4;
//				double indicator = deactivatingIndicator[changePosition][changeIndex];
//				
//				// is the change probably activating a deactivated splice site
//				if (Rref > Ralt && indicator >= -5) {
//					deactivatingVariants.add(variants.get(i));
//				} else {
//					counter ++;
//				}
//			}
//		}
//		Log.add("Number of activating variants: " + counter + "\t Number of variants left: " + deactivatingVariants.size(), 3);
//		return deactivatingVariants;
//	}

	/**
	 * @param variants2
	 * @param variantsBenign
	 * @return indicator > 0 -> probably an deactivating change
	 */
	/**
	 * @param variants
	 * @param variantsBenign
	 * @return
	 */
	private static double[][] calculateDeactivatingDifference(Variants variants, Variants variantsBenign) {
		int sequenceLength = variants.get(0).getSequence().length();
		// 16 change types {A->A, A->C, A->G, ..., T->T}
		double[][] deactivatingDifference = new double[sequenceLength][16];
		int[][] benignCount = new int[sequenceLength][16];
		int[][] pathogeneCount = new int[sequenceLength][16];
		int[] benignCountPosition = new int[sequenceLength];
		int[] pathogeneCountPosition = new int[sequenceLength];
		// initialize all with 0
		for (int l = 0; l < pathogeneCount.length; l++) {
			benignCountPosition[l] = 0;
			pathogeneCountPosition[l] = 0;
			for (int c = 0; c < pathogeneCount[l].length; c++) {
				benignCount[l][c] = 0;
				pathogeneCount[l][c] = 0;
			}
		}
		// count benign changes
		for (int i = 0; i < variantsBenign.size(); i++) {
			Variant variant = variantsBenign.get(i);
			int changePosition = variant.getSequence().getPositionChange();
			int changeIndex = Functions.mapNumber.get(variant.getRef().charAt(0)) + Functions.mapNumber.get(variant.getAlt().charAt(0)) * 4;
			if (changePosition > 0 && changePosition < sequenceLength) {
				benignCount[changePosition][changeIndex]++;
				benignCountPosition[changePosition]++;
			}
		}
		// count pathogene changes
		for (int i = 0; i < variants.size(); i++) {
			Variant variant = variants.get(i);
			int changePosition = variant.getSequence().getPositionChange();
			int changeIndex = Functions.mapNumber.get(variant.getRef().charAt(0)) + Functions.mapNumber.get(variant.getAlt().charAt(0)) * 4;
			if (changePosition > 0 && changePosition < sequenceLength) {
				pathogeneCount[changePosition][changeIndex]++;
				pathogeneCountPosition[changePosition]++;
			}
		}
		// calculate difference
		for (int l = 0; l < pathogeneCount.length; l++) {
			for (int c = 0; c < pathogeneCount[l].length; c++) {
				deactivatingDifference[l][c] = (200.0*pathogeneCount[l][c]) / variants.size() - (100.0*benignCount[l][c]) / variantsBenign.size();
			}
			System.out.println(Functions.arrayToString(deactivatingDifference[l], 5));
		}
		return deactivatingDifference;
	}
	
	/**
		 * TODO donor site
		 * @param variants
		 * @param acceptor
	 * @throws UnexpectedException 
		 */
		private static HashMap<String, Double> findPattern(Variants variants, boolean acceptor) throws UnexpectedException {
			int rangeLength = 5;
			int rangeShift = 3;
			variants = Filter.filterVariantType(variants, acceptor);
			variants = Filter.extractVariantsInRange(variants, -20, -5);
			variants = Filter.deleteDuplicateJunctions(variants);
			HashMap<String, Integer> quantityAbs = new HashMap<String, Integer>();
			HashMap<String, Integer> quantityCondition = new HashMap<String, Integer>();
			int numOfPattern[] = Functions.getInitializedIntArray(rangeLength/2 + 1);
			for (int i = 0; i < variants.size(); i++) {
				Sequence sequence = variants.get(i).getSequence();
				int posChange = sequence.getPositionChangeRelative();
				// count pattern on the variant position and in other sequences at the same position
				if (posChange >= -25 && posChange <= -5) {
					for (int shift = -rangeShift; shift <= rangeShift; shift++) {
						for (int add = Math.abs(shift); add <= rangeLength/2; add++) {
							if (posChange + shift + add < sequence.getPositionJunction() - 2) {
								String key = sequence.substring(posChange + shift - add, 1+posChange + shift + add);
	//							if (key.equals("GGACT")) {
	//								System.out.println("GAGAC:\n" + sequence.getVariant());
	//							}
								if (!quantityAbs.containsKey(key)) {
									quantityAbs.put(key, 1);
								} else {
									quantityAbs.put(key, quantityAbs.get(key) + 1);
								}
								// condition
								for (int j = 0; j < variants.size(); j++) {
									String subSequence = variants.get(j).getSequence().substring(posChange + shift - add, 1+posChange + shift + add, false);
									if (!quantityCondition.containsKey(subSequence)) {
										quantityCondition.put(subSequence, 1 / variants.size());
									} else {
										quantityCondition.put(subSequence, quantityCondition.get(subSequence) + 1 / variants.size());
									}
								}
								numOfPattern[add]++;
							}
						}
					}
				}
			}
			
			HashMap<String, Double> quantityRelative = AlgorithmAdministrator.relativeQuantity(quantityAbs, quantityCondition);
	////		System.out.println(relativeQuantity);
	//		ArrayList<Cluster> cluster = new ArrayList<Cluster>();
	//		Iterator<Entry<String, Double>> patternIt = quantityRelative.entrySet().iterator();
	//		// crate a cluster for relevant pattern and remove short ones
	//		while (patternIt.hasNext()) {
	//			Entry<String, Double> entry = patternIt.next();
	//			double quantity = (double) entry.getValue();
	//			String pattern = entry.getKey();
	//			if (quantity > 1.5) {
	//				Cluster clusterEntry = new Cluster(pattern);
	//				cluster.add(clusterEntry);
	//			}
	//			if (pattern.length() < 5) {
	//				quantityAbs.remove(pattern);
	//				quantityCondition.remove(pattern);
	//				patternIt.remove();
	//			}
	//		}
	//		// add similar pattern to cluster
	//		patternIt = quantityRelative.entrySet().iterator();
	//		while (patternIt.hasNext()) {
	//			Entry<String, Double> entry = patternIt.next();
	//			String pattern = entry.getKey();
	//		    boolean found = false;
	//		    if (entry.getValue() > 1) {
	//		    	for (int i = 0; i < cluster.size() && !found; i++) {
	//		    		if (equalsShifted(cluster.get(i).patternMain, pattern)) {
	//		    			cluster.get(i).add(pattern, quantityAbs.get(pattern), quantityCondition.get(pattern));
	////						patternIt.remove();
	////						found = true;
	//		    		}
	//		    	}
	//			}
	//		}
			return quantityRelative;
		}

	/**
		 * TODO donor site
		 * @param variants
		 * @param acceptor
	 * @throws UnexpectedException 
		 */
		private static HashMap<String, Double> findPatternGood(Variants variants, boolean acceptor) throws UnexpectedException {
			int lengthPatternMax = 5;
			variants = Filter.filterVariantType(variants, acceptor);
			variants = Filter.extractVariantsInRange(variants, -20, -5);
			variants = Filter.deleteDuplicateJunctions(variants);
			HashMap<String, Integer> quantityAbs = new HashMap<String, Integer>();
			HashMap<String, Integer> quantityCondition = new HashMap<String, Integer>();
			int numOfPattern[] = Functions.getInitializedIntArray(lengthPatternMax + 1);
			for (int i = 0; i < variants.size(); i++) {
				Sequence sequence = variants.get(i).getSequence();
				int posChangeRel = sequence.getPositionChangeRelative();
				int posChangeAbs = sequence.getPositionChange();
				// count pattern on the variant position and in other sequences at the same position
				if (posChangeRel >= -35 && posChangeRel <= -5) {
					for (int length = 1; length <= lengthPatternMax; length++) {
						for (int shift = 0; shift < length; shift++) {
							int from = posChangeAbs + shift - length + 1;
							int to = posChangeAbs + shift;
							int junction = sequence.getPositionJunction();
							if (Math.abs(from - junction) > 3	&& Math.abs(to - junction) > 3) {
								String pattern = sequence.substring(from, to + 1);
	//							if (key.equals("GGACT")) {
	//								System.out.println("GAGAC:\n" + sequence.getVariant());
	//							}
								if (!quantityAbs.containsKey(pattern)) {
									quantityAbs.put(pattern, 1);
								} else {
									quantityAbs.put(pattern, quantityAbs.get(pattern) + 1);
								}
								// condition
								String subSequence = sequence.substring(from, to + 1, false);
								if (!quantityCondition.containsKey(subSequence)) {
									quantityCondition.put(subSequence, 1);
								} else {
									quantityCondition.put(subSequence, quantityCondition.get(subSequence) + 1);
								}
								numOfPattern[length]++;
							}
						}
					}
				}
			}
			HashMap<String, Double> quantityRelative = AlgorithmAdministrator.relativeQuantity(quantityAbs, quantityCondition);
	////		System.out.println(relativeQuantity);
	//		ArrayList<Cluster> cluster = new ArrayList<Cluster>();
	//		Iterator<Entry<String, Double>> patternIt = quantityRelative.entrySet().iterator();
	//		// crate a cluster for relevant pattern and remove short ones
	//		while (patternIt.hasNext()) {
	//			Entry<String, Double> entry = patternIt.next();
	//			double quantity = (double) entry.getValue();
	//			String pattern = entry.getKey();
	//			if (quantity > 1.5) {
	//				Cluster clusterEntry = new Cluster(pattern);
	//				cluster.add(clusterEntry);
	//			}
	//			if (pattern.length() < 5) {
	//				quantityAbs.remove(pattern);
	//				quantityCondition.remove(pattern);
	//				patternIt.remove();
	//			}
	//		}
	//		// add similar pattern to cluster
	//		patternIt = quantityRelative.entrySet().iterator();
	//		while (patternIt.hasNext()) {
	//			Entry<String, Double> entry = patternIt.next();
	//			String pattern = entry.getKey();
	//		    boolean found = false;
	//		    if (entry.getValue() > 1) {
	//		    	for (int i = 0; i < cluster.size() && !found; i++) {
	//		    		if (equalsShifted(cluster.get(i).patternMain, pattern)) {
	//		    			cluster.get(i).add(pattern, quantityAbs.get(pattern), quantityCondition.get(pattern));
	////						patternIt.remove();
	////						found = true;
	//		    		}
	//		    	}
	//			}
	//		}
			return quantityRelative;
		}

}
