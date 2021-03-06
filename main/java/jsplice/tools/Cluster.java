/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
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
  private double InformationCore;
  private PatternMotif patternCore;
  private double[][] probability;
  private int variantSize = 1;

//  public Cluster(String patternP, int quantityAbsP, int quantityConditionP, int variantsSizeP) {
//    add(patternP, quantityAbsP, quantityConditionP);
//    this.variantSize = variantsSizeP;
//    this.patternCore = new PatternMotif(pattern.get(0));
//  }

  public Cluster(PatternMotif patternP, int variantsSizeP) {
    if (patternP == null) {
      throw new IllegalArgumentException("PatternMotif must not be null.");
    }
    this.patternCore = new PatternMotif(patternP);
    patternP.setCluster(this);
    this.variantSize = variantsSizeP;
    add(patternP, 0);
  }

  public Cluster(Cluster clusterP) {
    pattern = PatternMotif.clone(clusterP.pattern);
    patternSub = PatternMotif.clone(clusterP.patternSub);
    weightMatrix = clusterP.weightMatrix;
    lengthCluster = clusterP.lengthCluster;
    InformationCore = clusterP.InformationCore;
    this.variantSize = clusterP.variantSize;
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
  public void sortPattern(Comparator<PatternMotif> comp1, Comparator<PatternMotif> comp2) {
    Collections.sort(pattern, PatternComp.order(comp1, comp2));
    Collections.sort(patternSub, PatternComp.order(comp1, comp2));
  }

  /**
   * @return
   */
  public double[][] calculateInformationMatrix() {
    if (patternCore.getQuantityBen() < 1) {
      throw new IllegalStateException("There is no pattern quantityBen to calculate a matrix");
    }
    // count quantities by position
    int lengthPattern = getPatternCore().length();
    int lengthOverlapMax = Config.lengthIntronPatternMax; // TODO proper length
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
        if (sum[l] == 0 || count[l][b] == 0) {
          probability[l][b] = - 1000000;
          weightMatrix[l][b] = - 1000000;
        } else {
          probability[l][b] = count[l][b] / sum[l];
          weightMatrix[l][b] = 2.0 - (-Math.log(probability[l][b]) / Math.log(2) + error);
        }
      }
      //      Log.add("prob", 2);
      //      Log.add(Functions.arrayToString(probability[l], 1), 2);
    }
    InformationCore = getInformation(getPatternCore().pattern);
    return weightMatrix;
  }

  /**
   * TODO precalculate info for all pattern (contained in the HashMap)
   * @param patternP
   * @return
   */
  public double getInformation(String patternP) {
//    boolean multiRel = Config.multiClusterRel;
//    if (weightMatrix == null) {
//      calculateInformationMatrix();
//    }
//    if (patternP.length() > lengthCluster) {
//      throw new IllegalArgumentException("The length of the pattern Sequence (" + patternP.length()
//          + ") has to be equal or bigger than the length of this instance (" + lengthCluster + ").");
//    }
//    int lengthOverlapMax = Config.lengthIntronPatternMax;
//    // align pattern
//    int idx = pattern.indexOf(new PatternMotif(patternP, 1, 0, 0, null));
//    int align;
//    if (idx > -1) {
//      align = getPattern(idx).shift;
//    } else {
//      //        throw new IllegalArgumentException(this.getPatternCore() + "\n does not coontain " + patternP);
//      align = align(this, new PatternMotif(patternP, 1, 0, 0, null));
//    }
//    int matrixStart = lengthOverlapMax + align;
//    double[] individualInformation = Functions.getInitializedDoubleArray(patternP.length());
//    for (int lp = 0, lm = matrixStart; lp < patternP.length(); lp++, lm++) {
//      int baseNumber = Functions.mapNumber.get(patternP.charAt(lp));
//      individualInformation[lp] = weightMatrix[lm][baseNumber];
//    }
//    //    return Functions.sum(individualInformation) * Math.log(getPatternCore().quantityBen + 1) / Math.log(2);
//    if (multiRel) {
//      return Functions.sum(individualInformation) * Math.log10(getPatternCore().getQuantityBenRelative() + 1);
//    } else {
//      return Functions.sum(individualInformation);
//    }
    int len = getPatternCore().length();
    double rel = getPatternCore().getQuantityRefRelative();
    // error
//    double result = getPatternCore().length() * (2.0 - (4.0 - 1) / (2 * Math.log(2) * (rel - 1)));
    // similar log
//    double result = getPatternCore().length() * (-1+(Math.log(getPatternCore().getQuantityRefRelative()))/(Math.log(2)));
    // log at 1
    double result = len * ((Math.log(rel))/(Math.log(2)));
    return result;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public int compareTo(Cluster clusterP) {
    int quantityRelDelta =
        Double.compare(clusterP.getPatternCore().getQuantityRefRelative(),
            getPatternCore().getQuantityRefRelative());
    int quantityRefDelta =
        Integer.compare(clusterP.getPatternCore().quantityRef, getPatternCore().quantityRef);
    return quantityRelDelta != 0 ? quantityRelDelta : quantityRefDelta;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return "\n" + getPatternCore().pattern + " " + getQuantityBenRelative() + " " + getPatternCore().quantityRef + " " + getPatternCore().getQuantityRefRelative() + ":\n" + pattern;
  }

  public int size() {
    return pattern.size();
  }

  /**
   * @param patternP
   * @param quantityAbsP
   * @param quantityConditionP
   * @param shift
   * @return
   */
  public boolean add(String patternP, int quantityAbsP, int quantityConditionP) {
    return add(new PatternMotif(patternP, quantityAbsP, quantityConditionP, 0, this));
  }

  /**
   * @param patternP
   */
  public boolean add(PatternMotif patternP) {
    if (this.pattern.contains(patternP)) {
      return true;
    }
    if (pattern.size() > 0) {
      patternP.shift = align(this, patternP);
    } else if (patternP.containsOnce(getPatternCore())) {
      patternP.shift = getPatternCore().pattern.indexOf(patternP.pattern);
    } else {
      throw new IllegalArgumentException("Can't add PatternMotif " + patternP + " to Cluster "
          + this);
    }
    boolean success = pattern.add(patternP);
    if (success) {
      patternP.setCluster(this);
    }
    return success;
  }
  
  /**
   * @param patternP
   */
  public boolean add(PatternMotif patternP, int shift) {
    if (this.pattern.contains(patternP)) {
      return true;
    }
    patternP.shift = shift;
    boolean success = pattern.add(patternP);
    if (success) {
      patternP.setCluster(this);
    }
    return success;
  }

  /**
   * @param clusterP
   * @return
   */
  public boolean add(Cluster clusterP) {
    if (this == clusterP) {
      throw new IllegalArgumentException("Not possible to add cluster to itself.");
    }
    //    int align = Cluster.align(this, clusterP);
    //    if (align == Integer.MIN_VALUE) {
    //      return false;
    //    }
    if (pattern.size() > 20 || clusterP.size() > 15) {
      pattern.size();
    }
    boolean success = true;
    // add pattern
    for (int c = 0; c < clusterP.size(); c++) {
      PatternMotif patterCurrent = clusterP.getPattern(c);
      success &=
          add(patterCurrent);
    }
    // add sub pattern
    for (int c = 0; c < clusterP.sizeSub(); c++) {
      PatternMotif patternSub = clusterP.getPatternSub(c);
      double limit = Config.quantityRelLimit;
      if (patternSub.contains(getPatternCore()) && patternSub.getQuantityRefRelative() > limit) {
        success &=
            add(patternSub);
      } else {
        success &=
            addSub(clusterP.getPatternSub(c).pattern, clusterP.getPatternSub(c).quantityRef,
                clusterP.getPatternSub(c).quantityAlt);
      }
    }
    if(getPatternCore().getQuantityBen() > 0){
      calculateInformationMatrix();
    }
    return success;
  }

  public boolean addSub(String patternP, int quantityAbsP, int quantityConditionP) {
    return addSub(new PatternMotif(patternP, quantityAbsP, quantityConditionP, 0, this));
  }

  /**
   * @param patternP
   */
  public boolean addSub(PatternMotif patternP) {
    if (this.pattern.contains(patternP)) {
      return true;
    } else {
      boolean success = patternSub.add(patternP);
      if (success) {
        patternP.setCluster(this);
      }
      return success;
    }

  }

  /**
   * TODO remove variantSize param
   * find the longest pattern in the sorted cluster contains(String) will sort
   * @param sequenceP
   * @param posRel 
   */
  public boolean addQuantityBen(String sequenceP, int posRel) {
    for (int p = 0; p < pattern.size(); p++) {
      String patternStr = pattern.get(p).pattern;
      if (sequenceP.contains(patternStr)) {
        pattern.get(p).addQuantityBen(posRel);
        getPatternCore().addQuantityBen(posRel);
        return true;
      }
    }
    return false;
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

  public double getQuantityBenRelative() {
    if (getPatternCore().getQuantityBen() == 0) {
      return 0;
    }
    ArrayList<Integer> len = new ArrayList<>();
    for (int p = 0; p < pattern.size(); p++) {
      len.add(pattern.get(p).length() * pattern.get(p).getQuantityBen());
    }
    int lengthPatternMax = Config.lengthIntronPatternMax;
    double lengthPattern = Functions.mean(len) / getPatternCore().getQuantityBen();
    double possibleOccurances = lengthPatternMax - (lengthPattern - 1);
    double possibleCombinations = Math.pow(4, lengthPattern);
    double quantityBenRel = ((getPatternCore().getQuantityBen() - 1.) / variantSize) * (possibleCombinations / possibleOccurances);
    return quantityBenRel;
  }

  /**
   * 
   */
  public void resetQuantityBen() {
    getPatternCore().resetQuantityBenPositions();
    for (int p = 0; p < pattern.size(); p++) {
      pattern.get(p).resetQuantityBenPositions();
    }
  }

  /**
   * @param sequence
   * @return
   */
  public boolean isContainedBy(String sequenceStr) {
    for (PatternMotif patternMotif : pattern) {
      if (sequenceStr.contains(patternMotif.pattern)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Find the first matching pattern and return its relevance
   * @param sequenceStr
   * @return The relevance of the pattern depending on its length
   */
  public double getRatingMatch(String sequenceStr) {
    for (PatternMotif patternMotif : pattern) {
      if (sequenceStr.contains(patternMotif.pattern)) {
        return patternMotif.getRandomOccuranceProbability();
      }
    }
    return 0;
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

  /**
   * @return
   */
  public int getQuantityBen() {
    return getPatternCore().getQuantityBen();
  }

  /**
   * @return
   */
  public int getPosition() {
    ArrayList<Integer> positions = new ArrayList<Integer>();
    for (int p = 0; p < pattern.size(); p++) {
      positions.addAll(pattern.get(p).getPositions());
    }
    double pos = Functions.mean(positions);
    return (int) Functions.round(pos, 0);
  }

  public void realign() {
    ArrayList<PatternMotif> patternTemp = pattern;
    Collections.sort(patternTemp,
        PatternComp.order(PatternComp.desc(AttrPat.LENGTH), PatternComp.desc(AttrPat.QTY_BEN_REL)));
    pattern = new ArrayList<PatternMotif>();
    for (int p = 0; p < patternTemp.size(); p++) {
      try {
        add(patternTemp.get(p));
      } catch (IllegalArgumentException e) {
        patternTemp.add(patternTemp.get(p));
      }
    }
    calculateInformationMatrix();
  }


  /**
   * @return
   */
  public double variantSize() {
    return variantSize;
  }

  /**
   * @param cluster1
   * @param cluster2
   * @return shifts to align cluster2 to cluster1
   */
  public static int align(Cluster cluster1, PatternMotif pattern2) {
    int lengthOverlapMax = Config.lengthIntronPatternMax;
    PatternMotif patternC1 = cluster1.getPatternCore();
    Integer posBest = null;
    // pattern matches main pattern?
    if (patternC1.containsOnce(pattern2)) {
      posBest = patternC1.pattern.indexOf(pattern2.pattern);
    } else if (pattern2.containsOnce(patternC1)) {
      posBest = - pattern2.pattern.indexOf(patternC1.pattern);
    } else if (patternC1.getQuantityBen() > 0) {
      // align with redundancy matrix
      cluster1.calculateInformationMatrix();

      double similarityMax = Double.NEGATIVE_INFINITY;
      for (int pos = -pattern2.length() + 1; pos < patternC1.length(); pos++) {
        int from1 = (pos > 0 ? pos : 0) + lengthOverlapMax;
        int to1 = (pos + pattern2.length() < patternC1.length() ? pos + pattern2.length() : patternC1.length()) + lengthOverlapMax;
        int from2 = -pos > 0 ? -pos : 0;
        int to2 = -pos + patternC1.length() < pattern2.length() ? -pos + patternC1.length() : pattern2.length();
        double[][] subArray = Arrays.copyOfRange(cluster1.weightMatrix, from1, to1);
        String subString = pattern2.pattern.substring(from2, to2);
        double similarity = patternMatrixSimilarity(subString, subArray);
        if (similarity > similarityMax) {
          similarityMax = similarity;
          posBest = pos;
        } else if (similarity == similarityMax && pattern2.getPositions() != null) {
          try {
            double pos1 = cluster1.getPatternCore().getPosition();
            double pos2 = pattern2.getPosition();
            double difBest = pos2 - pos1;
            double oldDif = Math.abs(posBest - difBest);
            double newDif = Math.abs(pos - difBest);
            if (newDif < oldDif) {
              similarityMax = similarity;
              posBest = pos;
            }
          } catch (Exception e) {
//            if (similarity > 0 ) {
//              Log.add(cluster1.getPatternCore().pattern + " and " + pattern2.pattern + " have same similarity: " + similarity, 4);
//            }
          }
        }
      }
    } else if(cluster1.size() > 1) {
      // align by position
      HashMap<Integer, Integer> posCount = new HashMap<Integer, Integer>();
      for (int p = 0; p < cluster1.size(); p++) {
        if (cluster1.getPattern(p).contains(patternC1)) {
          PatternMotif pat = cluster1.getPattern(p);
          int pos = pat.pattern.indexOf(patternC1.pattern) + pat.shift;
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
            posBest = pos;
          }
        }                                                                        ;
      }
    }
    if (posBest == null) {
      if (patternC1.contains(pattern2)) {
        posBest = patternC1.pattern.indexOf(pattern2.pattern);
      } else if (pattern2.contains(patternC1)) {
        posBest = - pattern2.pattern.indexOf(patternC1.pattern);
      } else {
        posBest = PatternMotif.align(patternC1.pattern, pattern2.pattern);
      }
    }
    if (posBest == null) {
      throw new IllegalArgumentException("No good position found");
    }
    if (Math.abs(posBest) > 8){
      int x;
    }
    return posBest;
  }


  /**
   * @param matrix1
   * @param matrix2
   * @return
   */
  public static double patternMatrixSimilarity(String pattern1, double[][] matrix2) {
    if (pattern1.length() != matrix2.length) {
      throw new IllegalArgumentException("Both strings must have same length.");
    }
    double similarity = 0;
    for (int l = 0; l < pattern1.length(); l++) {
      int baseNumber = Functions.mapNumber.get(pattern1.charAt(l));
      double info = matrix2[l][baseNumber] + 100;
      if (info > 0) {
        similarity += info;
      }
    }
    return similarity;
  }

    /**
   * @param matrix1
   * @param matrix2
   * @return
   */
  public static double matrixMatrixSimilarity(double[][] matrix1, double[][] matrix2) {
    if (matrix1.length != matrix2.length) {
      throw new IllegalArgumentException("Both strings must have same length.");
    }
    double similarity = 0;
    for (int l = 0; l < matrix1.length; l++) {
      double similaritySum = 0;
      for (int b = 0; b < matrix1[l].length; b++) {
        double sumMatrix = matrix1[l][b] + matrix2[l][b];
        if (sumMatrix > 0) {
          similaritySum += sumMatrix;
        }
      }
      // System.out.println(Functions.arrayToString(matrix1[p], 1));
      //      System.out.println(Functions.arrayToString(matrix2[p], 1));
      //      System.out.println("add: " + distanceSum / matrix1[p].length);
      similarity += similaritySum / matrix1[l].length;
    }
    return similarity;
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