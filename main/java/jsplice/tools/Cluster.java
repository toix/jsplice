/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;
import java.util.Collections;

import jsplice.data.Config;
import jsplice.data.Pattern;
import jsplice.exception.Log;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Cluster {

	private ArrayList<Pattern> patternSub = new ArrayList<Pattern>();
	private ArrayList<Pattern> pattern = new ArrayList<Pattern>();
	
	private double[][] weightMatrix;
	private int lengthCluster;
	private int lengthOverlapMax;
	private double InformationCore;
	private Pattern patternCore;
	
	/**
	 * 
	 * @param patternP
	 * @param quantityAbsP
	 * @param quantityConditionP
	 */
	public Cluster(String patternP, int quantityAbsP, int quantityConditionP) {
		add(patternP, quantityAbsP, quantityConditionP);
		this.patternCore = getPattern(0);
	}
	
	/**
	 * 
	 * @param patternP
	 */
	public Cluster(Pattern patternP) {
		add(patternP);
		this.patternCore = patternP;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "\n" + getPatternCore() + ", " + getPatternCore().getQuantityRelative() + ":\n" + pattern + ",\n sub: " + patternSub;
	}
	
	public Pattern getPatternCore() {
		return patternCore;
	}
	
	public double getQuantityRelCore() {
		return getPatternCore().getQuantityRelative();
	}

	public double getQuantityAbs() {
		double sum = 0;
		for (int i = 0; i < pattern.size(); i++) {
			sum += getPattern(i).quantityAbs;
		}
		return sum;
	}
	
	public double getQuantityCondition() {
		double sum = 0;
		for (int i = 0; i < pattern.size(); i++) {
			sum += getPattern(i).quantityCon;
		}
		return sum;
	}
	
	public boolean addSub(String patternP, int quantityAbsP, int quantityConditionP) {
		return this.patternSub.add(new Pattern(patternP, quantityAbsP, quantityConditionP));
	}
	
	public boolean add(String patternP, int quantityAbsP, int quantityConditionP){
		return this.pattern.add(new Pattern(patternP, quantityAbsP, quantityConditionP));
	}
	
	public boolean add(Cluster clusterP){
		boolean success = true;
		for (int c = 0; c < clusterP.size(); c++) {
			success &= add(clusterP.getPattern(c).pattern, clusterP.getPattern(c).quantityAbs, clusterP.getPattern(c).quantityCon);
		}
		for (int c = 0; c < clusterP.sizeSub(); c++) {
			Pattern sub = clusterP.getPatternSub(c);
			if (sub.contains(getPatternCore())) {
				success &= add(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityAbs, clusterP.getPatternSub(c).quantityCon);
			} else {
				success &= addSub(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityAbs, clusterP.getPatternSub(c).quantityCon);
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
	
	public Pattern getPattern(int i) {
		return pattern.get(i);
	}
	
	public Pattern getPatternSub(int i) {
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
			double quantityUnique = getPattern(i).quantityUnique;
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
		return Functions.sum(individualInformation);// * Math.log(getPatternCore().getQuantityRelative()) / Math.log(2);
	}
	
	/**
	 * find the longest pattern in the sorted cluster
	 * contains(String) will sort
	 * @param patternP
	 */
	public void addQuantity(String patternP){
		boolean added = false;
		for (int p = 0; !added; p++) {
			String patternCurrent = pattern.get(p).pattern;
			if (patternP.contains(patternCurrent)) {
				pattern.get(p).quantityUnique++;
				added = true;
			}
		}
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
	public boolean add(Pattern patternP) {
		return pattern.add(patternP);
	}

	/**
	 * @param patternP
	 */
	public void addSub(Pattern patternP) {
		patternSub.add(patternP);
		
	}
	
//	public void addQuantityCondition(String patternP) {
//		if (pattern.contains(patternP)) {
//			int idx = pattern.indexOf(patternP);
//			quantityCondition.set(idx, quantityCondition.get(idx) + 1);
//		} else if (patternSub.contains(patternP)) {
//			int idx = patternSub.indexOf(patternP);
//			quantityConditionSub.set(idx, quantityConditionSub.get(idx) + 1);
//		}
//	}
}