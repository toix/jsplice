/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import jsplice.data.Config;
import jsplice.data.PatternMotif;
import jsplice.exception.Log;

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
  private double[][] probability;

  /**
   * 
   * @param patternP
   * @param quantityAbsP
   * @param quantityConditionP
   */
  public Cluster(String patternP, int quantityAbsP, int quantityConditionP) {

    add(patternP, quantityAbsP, quantityConditionP, 0);
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
   * @param patternP
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

  /**
   * sort pattern by (1)length and (2)quantityRel
   */
  public void sortPattern() {
    Collections.sort(pattern);
    Collections.sort(patternSub);
  }

  public double[][] calculateInformationMatrix() {
    // count quantities by position
    int lengthPattern = getPatternCore().length();
    lengthOverlapMax = Config.lengthIntronPatternMax; // TODO proper length
    lengthCluster = lengthOverlapMax + lengthPattern + lengthOverlapMax;
    int numberOfBases = Functions.bases.length();
    double[][] count = new double[lengthCluster][numberOfBases];
    int[] sum = new int[lengthCluster];
    for (int p = 0; p < pattern.size(); p++) {
      PatternMotif patternEntry = getPattern(p);
      int align = patternEntry.shift + lengthOverlapMax;
//      Log.add("align: " + align, 2);
//      System.out.println("align: " + align + "\t " + patternEntry);
      // System.out.println(patternEntry + "\t " + getPattern());
      double quantityBen = patternEntry.quantityBen;
      for (int b = 0, lm = align; b < patternEntry.length(); b++, lm++) {
        int baseIdx = Functions.mapNumber.get(patternEntry.pattern.charAt(b));
        // System.out.println("lm: " + lm + "\t idx: " + baseIdx);
        count[lm][baseIdx] += quantityBen;
        sum[lm] += quantityBen;
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
    probability = new double[lengthCluster][numberOfBases];
    weightMatrix = new double[lengthCluster][numberOfBases];
//    Log.add("count", 2);
    for (int l = 0; l < lengthCluster; l++) {
//      Log.add(Functions.arrayToString(count[l], 2));
      for (int b = 0; b < numberOfBases; b++) {
        double error = (4.0 - 1) / (2 * Math.log(2) * (sum[l] - 2)); // TODO check calculation
        probability[l][b] = count[l][b] / sum[l];
        weightMatrix[l][b] = 2.0 - (-Math.log(probability[l][b]) / Math.log(2) + error);
      }
//      Log.add("prob", 2);
//      Log.add(Functions.arrayToString(probability[l], 1), 2);
    }
    InformationCore = getInformation(getPatternCore().pattern);
//     Debug
    if (getPatternCore().equals("ATC")) {
      System.out.println(toString());
      System.out.println("\nprob");
//      for (int l = 0; l < lengthCluster; l++) {
//        System.out.println(Functions.arrayToString(countWeighted[l], 2));
//      }
      System.out.println("\nmatrix");
      for (int l = 0; l < lengthCluster; l++) {
        System.out.println(Functions.arrayToString(weightMatrix[l], 2));
      }
    }
    return weightMatrix;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public int compareTo(Cluster clusterP) {
    int quantityBenDelta =
        Integer.compare(clusterP.patternCore.quantityBen, patternCore.quantityBen);
    int quantityRelDelta =
        Double.compare(clusterP.patternCore.getQuantityRelative(),
            patternCore.getQuantityRelative());
    return quantityBenDelta != 0 ? quantityBenDelta : quantityRelDelta;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return "\n" + getPatternCore() + ":\n" + pattern + ",\n sub: " + patternSub;
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
   * @param patternP
   * @param quantityAbsP
   * @param quantityConditionP
   * @param shift
   * @return
   */
  public boolean add(String patternP, int quantityAbsP, int quantityConditionP, int shift) {
    return add(new PatternMotif(patternP, quantityAbsP, quantityConditionP, shift));
  }

  /**
   * @param patternP
   */
  public boolean add(PatternMotif patternP) {
    if (this.pattern.contains(patternP)){
      return true;
    } else {
      return pattern.add(patternP);
    }
  }

  public boolean add(Cluster clusterP) {
    if (this == clusterP) {
      throw new IllegalArgumentException("Not possible to add cluster to itself.");
    }
    Log.add("align " + this.getPatternCore() + " to " + clusterP.getPatternCore(), 2);
    int shift = Cluster.align(this, clusterP);
    boolean success = true;
    // add pattern
    for (int c = 0; c < clusterP.size(); c++) {
      PatternMotif patterCurrent = clusterP.getPattern(c);
      success &=
          add(patterCurrent.pattern, patterCurrent.quantityRef,
              patterCurrent.quantityAlt, patterCurrent.shift + shift);
    }
    // add sub pattern
    for (int c = 0; c < clusterP.sizeSub(); c++) {
      PatternMotif patternSub = clusterP.getPatternSub(c);
      double limit = Config.quantityRelLimit;
      if (patternSub.contains(getPatternCore()) && patternSub.getQuantityRelative() > limit) {
        success &=
            add(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef,
                clusterP.getPatternSub(c).quantityAlt, patternSub.shift + shift);
      } else {
        success &=
            addSub(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef,
                clusterP.getPatternSub(c).quantityAlt);
      }
    }
    if(getPatternCore().quantityBen > 0){
      calculateInformationMatrix();
    }
    return success;
  }

  public boolean addSub(String patternP, int quantityAbsP, int quantityConditionP) {
    return addSub(new PatternMotif(patternP, quantityAbsP, quantityConditionP, 0));
  }

  /**
   * @param patternP
   */
  public boolean addSub(PatternMotif patternP) {
    if (this.pattern.contains(patternP)) {
      return true;
    } else {
      return patternSub.add(patternP);
    }

  }

  /**
   * find the longest pattern in the sorted cluster contains(String) will sort
   * 
   * @param patternP
   */
  public boolean addQuantityBen(String patternP) {
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
   * @param sequenceAlt
   * @return
   */
  public boolean addQuantityPat(String patternP) {
    boolean added = false;
    for (int p = 0; p < pattern.size() && !added; p++) {
      String patternCurrent = pattern.get(p).pattern;
      if (patternP.contains(patternCurrent)) {
        pattern.get(p).quantityPat++;
        getPatternCore().quantityPat++;
        added = true;
      }
    }
    return added;
  }

  public int size() {
    return pattern.size();
  }

  public int sizeSub() {
    return patternSub.size();
  }

  /**
   * @return
   */
  public int length() {
    return lengthCluster;
  }

  public boolean contains(String patternP) {
    return pattern.contains(patternP);
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

  public PatternMotif getPattern(int i) {
    return pattern.get(i);
  }

  public PatternMotif getPatternSub(int i) {
    return patternSub.get(i);
  }
  

//  /**
//   * @return
//   */
//  public int getQuantityBenign() {
//    int quantity = 0;
//    for (PatternMotif patternCurrent : pattern) {
//      quantity += patternCurrent.quantityBen;
//    }
//    return quantity;
//  }


  /**
   * TODO precalculate info for all pattern in HashMap
   * @param patternP
   * @return
   */
  public double getInformation(String patternP) {
    if (weightMatrix == null) {
      calculateInformationMatrix();
    }
    if (patternP.length() > lengthCluster) {
      throw new IllegalArgumentException("The length of the pattern Sequence (" + patternP.length()
          + ") has to be equal or bigger than the length of this instance (" + lengthCluster + ").");
    }
    // align pattern
    int idx = pattern.indexOf(new PatternMotif(patternP, 1, 0, 0));
    if (idx == -1) {
      throw new IllegalArgumentException(this.getPatternCore() + "\n does not coontain " + patternP);
    }
    PatternMotif patternMatching = getPattern(idx);
    int matrixStart = lengthOverlapMax + patternMatching.shift;
    double[] individualInformation = Functions.getInitializedDoubleArray(patternP.length());
    for (int lp = 0, lm = matrixStart; lp < patternP.length(); lp++, lm++) {
      int baseNumber = Functions.mapNumber.get(patternP.charAt(lp));
      individualInformation[lp] = weightMatrix[lm][baseNumber];
    }
    return Functions.sum(individualInformation);
//    return Functions.sum(individualInformation) * Math.log(getPatternCore().getQuantityRelative()) / Math.log(2);
//    return Functions.sum(individualInformation) * Math.log(getPatternCore().quantityBen) / Math.log(2);
  }

  /**
   * @param cluster1
   * @param cluster2
   * @return shifts to align cluste2 to cluster1
   */
  public static int align(Cluster cluster1, Cluster cluster2) {
    PatternMotif patternC1 = cluster1.getPatternCore();
    PatternMotif patternC2 = cluster2.getPatternCore();
    int posBest = PatternMotif.align(patternC1.pattern, patternC2.pattern);
    if (patternC1.contains(patternC2)) {
      posBest = patternC1.pattern.indexOf(patternC2.pattern);
    } else if (patternC2.contains(patternC1)) {
      posBest = - patternC2.pattern.indexOf(patternC1.pattern);
    } else if (patternC1.quantityBen > 0 && patternC2.quantityBen > 0) {
      if (cluster1.weightMatrix == null) {
        cluster1.calculateInformationMatrix();
        Log.add(cluster1.matrixToString("matrix1"), 2);
      }
      if (cluster2.weightMatrix == null) {
        cluster2.calculateInformationMatrix();
        Log.add(cluster2.matrixToString("matrix2"), 2);
      }
      double similarityBest = 0;
      //     String str1 = pattern1.pattern;
      //     String str2 = pattern2.pattern;
      for (int pos = cluster1.length() - 1; pos > -cluster2.length(); pos--) {
        int from1 = pos > 0 ? pos : 0;
        int to1 = pos + cluster2.length() < cluster1.length() ? pos + cluster2.length() : cluster1.length();
        int from2 = -pos > 0 ? -pos : 0;
        int to2 = -pos + cluster1.length() < cluster2.length() ? -pos + cluster1.length() : cluster2.length();
        //      System.out.println(from1 + " > " + to1 + " of " + str1);
        //      System.out.println(from2 + " > " + to2 + " of " + str2);
        double[][] subArray1 = Arrays.copyOfRange(cluster1.probability, from1, to1);
        double[][] subArray2 = Arrays.copyOfRange(cluster2.probability, from2, to2);
        double similarity = Functions.probabilityMatrixSimilarity(subArray1, subArray2);
        if (similarity > similarityBest) {
          similarityBest = similarity;
          posBest = pos;
        } else if (similarity == similarityBest) {
          //TODO save pos and add to log if matches = matchesBest best after loop
          Log.add("alignment has same quality: " + similarity, 4);
        }
      }
    }
    Log.add("posBest: " + posBest, 2);
    return posBest;
  }

  /**
   * @param clusterP
   * @return
   */
  public static ArrayList<Cluster> clone(ArrayList<Cluster> list) {
    ArrayList<Cluster> clone = new ArrayList<Cluster>(list.size());
    for (Cluster item : list) {
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
    for (Cluster item : list) {
      clone.put(item, new Cluster(item));
    }
    return clone;
  }
}
