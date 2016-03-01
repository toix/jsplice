/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

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
  private int variantsSize = 1;

  /**
   * 
   * @param patternP
   * @param quantityAbsP
   * @param quantityConditionP
   */
  public Cluster(String patternP, int quantityAbsP, int quantityConditionP, int variantsSizeP) {
    add(patternP, quantityAbsP, quantityConditionP, 0);
    this.variantsSize = variantsSizeP;
    this.patternCore = new PatternMotif(pattern.get(0));
  }

  /**
   * 
   * @param patternP
   */
  public Cluster(PatternMotif patternP, int variantsSizeP) {
    add(patternP);
    this.variantsSize = variantsSizeP;
    this.patternCore = new PatternMotif(pattern.get(0));
  }

  public Cluster(Cluster clusterP) {
    pattern = PatternMotif.clone(clusterP.pattern);
    patternSub = PatternMotif.clone(clusterP.patternSub);
    weightMatrix = clusterP.weightMatrix;
    lengthCluster = clusterP.lengthCluster;
    lengthOverlapMax = clusterP.lengthOverlapMax;
    InformationCore = clusterP.InformationCore;
    this.variantsSize = clusterP.variantsSize;
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

  /**
   * TODO NaN
   * @return
   */
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
      double quantityBen = patternEntry.getQuantityBen();
      for (int l = 0, lm = align; l < patternEntry.length(); l++, lm++) {
        int baseIdx = Functions.mapNumber.get(patternEntry.pattern.charAt(l));
        // System.out.println("lm: " + lm + "\t idx: " + baseIdx);
        count[lm][baseIdx] += quantityBen;
        sum[lm] += quantityBen;
      }
    }
    // divide probability sums through sum and calculate Information Matrix
    probability = new double[lengthCluster][numberOfBases];
    weightMatrix = new double[lengthCluster][numberOfBases];
//    Log.add("count", 2);
    for (int l = 0; l < lengthCluster; l++) {
//      Log.add(Functions.arrayToString(count[l], 2));
      for (int b = 0; b < numberOfBases; b++) {
        double error = (4.0 - 1) / (2 * Math.log(2) * (sum[l] - 1)); // TODO check calculation
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
    int quantityRelDelta =
        Double.compare(clusterP.patternCore.getQuantityRefRelative(),
            patternCore.getQuantityRefRelative());
    int quantityRefDelta =
        Integer.compare(clusterP.patternCore.quantityRef, patternCore.quantityRef);
    return quantityRelDelta != 0 ? quantityRelDelta : quantityRefDelta;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return "\n" + getPatternCore().pattern + " " + getQuantityBenRel(variantsSize) + " " + getPatternCore().quantityRef + " " + getPatternCore().getQuantityRefRelative() + ":\n" + pattern;
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

  /**
   * @param clusterP
   * @return
   */
  public boolean add(Cluster clusterP) {
    if (this == clusterP) {
      throw new IllegalArgumentException("Not possible to add cluster to itself.");
    }
//    Log.add("align " + this.getPatternCore() + " to " + clusterP.getPatternCore(), 2);
    int shift = Cluster.align(this, clusterP);
    if (shift == Integer.MIN_VALUE) {
      return false;
    }
    if (Math.abs(shift) > 5 || pattern.size() > 15 || clusterP.size() > 15) {
      shift = shift;
    }
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
      if (patternSub.contains(getPatternCore()) && patternSub.getQuantityRefRelative() > limit) {
        success &=
            add(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef,
                clusterP.getPatternSub(c).quantityAlt, patternSub.shift + shift);
      } else {
        success &=
            addSub(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef,
                clusterP.getPatternSub(c).quantityAlt);
      }
    }
    if(getPatternCore().getQuantityBen() > 0){
      calculateInformationMatrix();
    }
    sortPattern();
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
   * @param sequenceP
   * @param posRel 
   */
  public boolean addQuantityBen(String sequenceP, int posRel, int variantsSize) {
    boolean added = false;
    sortPattern();
    for (int p = 0; p < pattern.size() && !added; p++) {
      String patternStr = pattern.get(p).pattern;
      if (sequenceP.contains(patternStr)) {
        pattern.get(p).addQuantityBen();
        getPatternCore().addQuantityBen();
        pattern.get(p).addPos(posRel);
        getPatternCore().addPos(posRel);
        added = true;
      }
    }
    return added;
  }

//  /**
//   * @param sequenceAlt
//   * @return
//   */
//  public boolean addQuantityPat(String patternP) {
//    boolean added = false;
//    for (int p = 0; p < pattern.size() && !added; p++) {
//      String patternCurrent = pattern.get(p).pattern;
//      if (patternP.contains(patternCurrent)) {
//        pattern.get(p).quantityPat++;
//        getPatternCore().quantityPat++;
//        added = true;
//      }
//    }
//    return added;
//  }

  public double getQuantityBenRel(int variantsSize) {
    double quantityBenRel = 0;
    for (int p = 0; p < pattern.size(); p++) {
      quantityBenRel += pattern.get(p).getQuantityBenRelative(variantsSize);
    }
    return quantityBenRel;
  }

  /**
   * 
   */
  public void resetQuantityBen() {
    getPatternCore().resetQuantityBen();
    getPatternCore().resetPositions();
    for (int p = 0; p < pattern.size(); p++) {
      pattern.get(p).resetQuantityBen();
      pattern.get(p).resetPositions();
    }
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


//  /**
//   * @return
//   */
//  public int getPosition() {
//    ArrayList<Integer> positions = new ArrayList<Integer>();
//    for (int p = 0; p < pattern.size(); p++) {
//      positions.addAll(pattern.get(p).getPositions());
//    }
//    double pos = Functions.mean(positions);
//    return (int) Functions.round(pos, 0);
//  }

  /**
   * TODO precalculate info for all pattern (contained in the HashMap)
   * @param patternP
   * @return
   */
  public double getInformation(String patternP) {
    boolean multiRel = Config.multiClusterRel;
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
//    return Functions.sum(individualInformation) * Math.log(getPatternCore().quantityBen + 1) / Math.log(2);
    if (multiRel) {
      return Functions.sum(individualInformation) * Math.log(getPatternCore().getQuantityRefRelative() + 1) / Math.log(2);
    } else {
      return Functions.sum(individualInformation);
    }
  }

  /**
   * @param sequence
   * @return
   */
  public boolean isContainedBy(String sequence) {
    if (sequence.contains(getPatternCore().pattern)) {
      return true;
    }
    for (PatternMotif patternMotif : pattern) {
      if (sequence.contains(patternMotif.pattern)) {
        return true;
      }
    }
    return false;
  }

  /**
     * @param cluster1
     * @param cluster2
     * @return shifts to align cluste2 to cluster1
     */
    public static int align(Cluster cluster1, Cluster cluster2) {
      PatternMotif patternC1 = cluster1.getPatternCore();
      PatternMotif patternC2 = cluster2.getPatternCore();
      int posBest = Integer.MIN_VALUE;
      // pattern main matches?
      if (patternC1.contains(patternC2)) {
        posBest = patternC1.pattern.indexOf(patternC2.pattern);
      } else if (patternC2.contains(patternC1)) {
        posBest = - patternC2.pattern.indexOf(patternC1.pattern);
      } else if (patternC1.getQuantityBen() > 1 && patternC2.getQuantityBen() > 1) {
        // align with redundancy matrix
        cluster1.calculateInformationMatrix();
        cluster2.calculateInformationMatrix();
        double similarityBest = 0;
        //     String str1 = pattern1.pattern;
        //     String str2 = pattern2.pattern;
        for (int pos = -cluster2.length() + 1; pos < cluster1.length(); pos++) {
          int from1 = pos > 0 ? pos : 0;
          int to1 = pos + cluster2.length() < cluster1.length() ? pos + cluster2.length() : cluster1.length();
          int from2 = -pos > 0 ? -pos : 0;
          int to2 = -pos + cluster1.length() < cluster2.length() ? -pos + cluster1.length() : cluster2.length();
          //      System.out.println(from1 + " > " + to1 + " of " + str1);
          //      System.out.println(from2 + " > " + to2 + " of " + str2);
          double[][] subArray1 = Arrays.copyOfRange(cluster1.probability, from1, to1);
          double[][] subArray2 = Arrays.copyOfRange(cluster2.probability, from2, to2);
          double similarity = Functions.probabilityMatrixSimilarity(subArray1, subArray2);
          if (similarity > similarityBest && similarity > 0) {
            similarityBest = similarity;
            posBest = pos;
          } else if (similarity == similarityBest  && similarity > 0) {
            double pos1 = cluster1.getPatternCore().getPosition();
            double pos2 = cluster2.getPatternCore().getPosition();
            double posDif = pos + pos2 - pos1;
            if (posDif < 0) {
              similarityBest = similarity;
              posBest = pos;
            } else if (posDif > 0) {
            } else {
              Log.add(cluster1.getPatternCore() + " and " + cluster2.getPatternCore() + " have same similarity: " + similarity, 4);
            }
          }
        }
      } else if(cluster1.size() > 1) {
        // align by position
        HashMap<Integer, Integer> posCount = new HashMap<Integer, Integer>();
        for (int p = 0; p < cluster1.size(); p++) {
          if (cluster1.getPattern(p).contains(patternC1)) {
            int pos = cluster1.getPattern(p).pattern.indexOf(patternC1.pattern);
            if (posCount.containsKey(pos)) {
              posCount.put(pos, posCount.get(pos) + 1);
            } else {
              posCount.put(pos, 1);
            }
          } else if (cluster1.getPattern(p).contains(patternC1)) {
            int pos = - patternC1.pattern.indexOf(cluster1.getPattern(p).pattern);
            if (posCount.containsKey(pos)) {
              posCount.put(pos, posCount.get(pos) + 1);
            } else {
              posCount.put(pos, 1);
            }
          }
        }
        int maxValue = Collections.max(posCount.values());
        if (maxValue > 1) {
          for (Integer pos : posCount.keySet()) {
            if (posCount.get(pos).equals(maxValue)) {
              return pos;
            }
          }                                                                        ;
        } else {
          if (patternC1.contains(patternC2)) {
            posBest = patternC1.pattern.indexOf(patternC2.pattern);
          } else if (patternC2.contains(patternC1)) {
            posBest = - patternC2.pattern.indexOf(patternC1.pattern);
          }
        }
      }
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
   * @param clusterP
   * @return
   */
  public static ArrayList<Cluster> cloneSoft(ArrayList<Cluster> list) {
    ArrayList<Cluster> clone = new ArrayList<Cluster>(list.size());
    for (Cluster item : list) {
      clone.add(item);
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
