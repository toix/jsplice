/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import jsplice.data.Config;
import jsplice.data.PatternMotif;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Cluster implements Comparable<Cluster> {

	private ArrayList<PatternMotif> pattern = new ArrayList<PatternMotif>();
	private ArrayList<PatternMotif> patternSub = new ArrayList<PatternMotif>();
	private double[][] weightMatrix;
	private int lengthCluster;
	private int lengthOverlapMax;
	private double InformationCore;
	private PatternMotif patternCore;
	
	/**
	 * 
	 * @param patternP
	 * @param quantityAbsP
	 * @param quantityConditionP
	 */
	public Cluster(String patternP, int quantityAbsP, int quantityConditionP) {
		add(patternP, quantityAbsP, quantityConditionP);
		this.patternCore = new PatternMotif(pattern.get(0));
	}
	
	/**
	 * 
	 * @param patternP
	 */
	public Cluster(PatternMotif patternP) {
		add(patternP);
		this.patternCore = new PatternMotif(pattern.get(0));
	}
	
	public Cluster(Cluster clusterP) {
		pattern = PatternMotif.clone(clusterP.pattern);
		patternSub = PatternMotif.clone(clusterP.patternSub);
		weightMatrix = clusterP.weightMatrix;
		lengthCluster = clusterP.lengthCluster;
		lengthOverlapMax = clusterP.lengthOverlapMax;
		InformationCore = clusterP.InformationCore;
		this.patternCore = findPatternCore(clusterP.getPatternCore());
	}

	/**
	 * @param pattern2
	 * @return
	 */
	private PatternMotif findPatternCore(PatternMotif patternP) {
		for (int i = 0; i < pattern.size(); i++) {
			if (pattern.get(i).equals(patternP)) {
				return new PatternMotif(pattern.get(i));
			}
		}
		return null;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "\n" + getPatternCore() + ":\n" + pattern + ",\n sub: " + patternSub;
	}
	
	public PatternMotif getPatternCore() {
		return patternCore;
	}

	public double getQuantityRef() {
		double sum = 0;
		for (int i = 0; i < pattern.size(); i++) {
			sum += getPattern(i).quantityRef;
		}
		return sum;
	}
	
	public double getQuantityAlt() {
		double sum = 0;
		for (int i = 0; i < pattern.size(); i++) {
			sum += getPattern(i).quantityAlt;
		}
		return sum;
	}
	
	public boolean addSub(String patternP, int quantityAbsP, int quantityConditionP) {
		return this.patternSub.add(new PatternMotif(patternP, quantityAbsP, quantityConditionP));
	}
	
	public boolean add(String patternP, int quantityAbsP, int quantityConditionP){
		return this.pattern.add(new PatternMotif(patternP, quantityAbsP, quantityConditionP));
	}
	
	public boolean add(Cluster clusterP){
		boolean success = true;
		for (int c = 0; c < clusterP.size(); c++) {
			success &= add(clusterP.getPattern(c).pattern, clusterP.getPattern(c).quantityRef, clusterP.getPattern(c).quantityAlt);
		}
		for (int c = 0; c < clusterP.sizeSub(); c++) {
			PatternMotif sub = clusterP.getPatternSub(c);
			double limit = Config.quantityRelLimit;
			if (sub.contains(getPatternCore()) && sub.getQuantityRelative() > limit) {
				success &= add(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef, clusterP.getPatternSub(c).quantityAlt);
			} else {
				success &= addSub(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef, clusterP.getPatternSub(c).quantityAlt);
			}
		}
		return success;
	}
	
	public int size() {
		return pattern.size();
	}
	
	public int sizeSub() {
		return patternSub.size();
	}
	
	public PatternMotif getPattern(int i) {
		return pattern.get(i);
	}
	
	public PatternMotif getPatternSub(int i) {
		return patternSub.get(i);
	}
	
	public boolean contains(String patternP) {
		return pattern.contains(patternP);
	}

	public double[][] calculateInformationMatrix(){
		// count quantities by position
		int lengthPattern = getPatternCore().length();
		lengthOverlapMax = Config.lengthIntronPatternMax - lengthPattern;
		lengthCluster = lengthOverlapMax + lengthPattern + lengthOverlapMax;
		int numberOfBases = Functions.bases.length();
		double[][] count = new double[lengthCluster][numberOfBases];
		int[] sum = new int[lengthCluster];
		for (int i = 0; i < pattern.size(); i++) {
			String patternEntry = getPattern(i).pattern;
			int align = patternEntry.indexOf(getPatternCore().pattern);
//			System.out.println(patternEntry + "\t " + getPattern());
			double quantityUnique = getPattern(i).quantityBen;
			for (int lp = 0, lm = lengthOverlapMax - align; lp < patternEntry.length(); lp++, lm++) {
				int baseIdx = Functions.mapNumber.get(patternEntry.charAt(lp));
//				System.out.println("lm: " + lm + "\t idx: " + baseIdx);
				count[lm][baseIdx] += quantityUnique;
				sum[lm] += quantityUnique;
			}
		}
		// invert count with count sum
		int[] sumInverted = new int[lengthCluster];
		int sumMax = sum[lengthOverlapMax + 1];
		for (int l = 0; l < lengthCluster; l++) {
			sumInverted[l] = sumMax - sum[l];
		}
		// weighted probability
		double[][] countWeighted = new double[lengthCluster][numberOfBases];
		double[] sumWeighted = new double[lengthCluster];
		for (int i = 0; i < pattern.size(); i++) {
			String patternEntry = getPattern(i).pattern;
			int align = patternEntry.indexOf(getPatternCore().pattern);
			double quantityRelativeLog = Math.log(getPattern(i).getQuantityRelative()) / Math.log(2);
			for (int lp = 0, lm = lengthOverlapMax - align; lp < patternEntry.length(); lp++, lm++) {
				int baseIdx = Functions.mapNumber.get(patternEntry.charAt(lp));
				countWeighted[lm][baseIdx] += quantityRelativeLog;
				sumWeighted[lm] += quantityRelativeLog;
			}
		}
		// divide probability sums through sum and calculate Information Matrix
		weightMatrix = new double[lengthCluster][numberOfBases];
		for (int l = 0; l < lengthCluster; l++) {
//			System.out.println(Functions.arrayToString(probability[l], 2));
			for (int b = 0; b < numberOfBases; b++) {
				double error = (4.0 - 1) / (2 * Math.log(2) * (sum[l]-2));	// TODO check calculation
				double probability = count[l][b] / sum[l];
				weightMatrix[l][b] = 2.0 - (-Math.log(probability) / Math.log(2) + error);
			}
//			Log.add(Functions.arrayToString(weightMatrix[l], 1), 2);
		}
		InformationCore = getInformation(getPatternCore().pattern);
		// Debug
//		if (getPattern().equals("CAAT")) {
//			System.out.println(toString());
//			System.out.println("\nprob");
//			for (int l = 0; l < lengthCluster; l++) {
//				System.out.println(Functions.arrayToString(countWeighted[l], 2));
//			}
//			System.out.println("\nmatrix");
//			for (int l = 0; l < lengthCluster; l++) {
//				System.out.println(Functions.arrayToString(weightMatrix[l], 2));
//			}
//		}
		return weightMatrix;
	}

	/**
	 * @param patternP
	 * @return
	 */
	public double getInformation(String patternP) {
		String patternCore = getPatternCore().pattern;
		if (weightMatrix == null) {
			calculateInformationMatrix();
		}
		if (patternP.length() > lengthCluster) {
			throw new IllegalArgumentException("The length of the pattern Sequence (" + patternP.length()
					+ ") has to be equal or bigger than the length of this instance (" + lengthCluster + ").");
		}
		int matrixStart = 0;
		if (getPatternCore().contains(patternP)) {
			matrixStart = lengthOverlapMax + patternCore.indexOf(patternP);
		} else if (patternP.contains(patternCore)) {
			matrixStart = lengthOverlapMax - patternP.indexOf(patternCore);
		} else {
			System.out.println("Bad pattern: " + patternP + " for cluster " + getPatternCore());
		}
		double[] individualInformation = Functions.getInitializedDoubleArray(patternP.length());
		for (int l1 = 0, l2 = matrixStart; l1 < patternP.length(); l1++, l2++) {
			int baseNumber = Functions.mapNumber.get(patternP.charAt(l1));
			individualInformation[l1] = weightMatrix[l2][baseNumber];
		}
//		return Functions.sum(individualInformation);
//		return Functions.sum(individualInformation) * Math.log(getPatternCore().getQuantityRelative()) / Math.log(2);
		return Functions.sum(individualInformation) * Math.log(getQuantityBenign()) / Math.log(2);
	}
	
	/**
	 * @return
	 */
	public int getQuantityBenign() {
		int quantity = 0;
		for (PatternMotif patternCurrent : pattern) {
			quantity += patternCurrent.quantityBen;
		}
		return quantity;
	}

	/**
	 * find the longest pattern in the sorted cluster
	 * contains(String) will sort
	 * @param patternP
	 */
	public boolean addQuantityBen(String patternP){
		boolean added = false;
		for (int p = 0; p < pattern.size() && !added; p++) {
			String patternCurrent = pattern.get(p).pattern;
			if (patternP.contains(patternCurrent)) {
				pattern.get(p).quantityBen++;
				getPatternCore().quantityBen++;
				added = true;
			}
		}
		return added;
	}

	/**
	 * sort pattern by (1)length and (2)quantityRel
	 */
	public void sortPattern() {
		Collections.sort(pattern);
		Collections.sort(patternSub);
	}

	/**
	 * @param patternP
	 */
	public boolean add(PatternMotif patternP) {
		return pattern.add(patternP);
	}

	/**
	 * @param patternP
	 */
	public void addSub(PatternMotif patternP) {
		patternSub.add(patternP);
		
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public int compareTo(Cluster clusterP) {
		int quantityBenDelta = Integer.compare(clusterP.patternCore.quantityBen, patternCore.quantityBen);
		int quantityRelDelta = Double.compare(clusterP.patternCore.getQuantityRelative(), patternCore.getQuantityRelative());
		return quantityBenDelta != 0 ? quantityBenDelta : quantityRelDelta;
	}

	/**
	 * @param clusterP
	 * @return
	 */
	public static ArrayList<Cluster> clone(ArrayList<Cluster> list) {
		ArrayList<Cluster> clone = new ArrayList<Cluster>(list.size());
	    for(Cluster item: list) {
	    	clone.add(new Cluster(item));
	    }
	    return clone;
	}

	/**
	 * @param list
	 * @return
	 */
	public static HashMap<Cluster, Cluster> cloneToMap(ArrayList<Cluster> list) {
		HashMap<Cluster, Cluster> clone = new HashMap<Cluster, Cluster>(list.size());
	    for(Cluster item: list) {
	    	clone.put(item, new Cluster(item));
	    }
	    return clone;
	}
}