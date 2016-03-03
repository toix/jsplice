/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import jsplice.exception.Log;
import jsplice.tools.Cluster;
import jsplice.tools.Functions;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class PatternMotif implements Comparable<PatternMotif>{

  public String pattern;
  public int quantityRef;
  public int quantityAlt;
  public int quantityPat;
  public int shift;
  private ArrayList<Integer> positions;
  private Cluster cluster;

  /**
   * 
   */
  public PatternMotif(String patternP, int quantityAbsP, int quantityConP, int shiftP, Cluster clusterP) {
    this.pattern = patternP;
    this.quantityRef = quantityAbsP;
    this.quantityAlt = quantityConP;
    this.shift = shiftP;
    this.cluster = clusterP;
    positions = new ArrayList<Integer>();
  }

  @SuppressWarnings("unchecked")
  public PatternMotif(PatternMotif patternP) {
    this.pattern = patternP.pattern;
    this.quantityRef = patternP.quantityRef;
    this.quantityAlt = patternP.quantityAlt;
    this.quantityPat = patternP.quantityPat;
    this.shift = patternP.shift;
    this.cluster = patternP.cluster;
    if (positions != null) {
      this.positions = (ArrayList<Integer>) patternP.positions.clone();
    } else {
      positions = new ArrayList<Integer>();
    }
  }
  /* (non-Javadoc)
   * @see java.lang.Object#toString()
   */
  @Override
  public String toString() {
    return pattern + " " + getQuantityBen() + " " + getQuantityBenRelative() + " " + quantityRef + " " + getQuantityRefRelative();
  }
  
  /* (non-Javadoc)
   * @see java.lang.Object#equals()
   */
  @Override
  public boolean equals(Object patternP) {
    if (patternP == null) {
      return false;
    } else if (patternP == this) {
      return true;
    }
    else if ((patternP instanceof PatternMotif)) {
      PatternMotif patternMotif = (PatternMotif) patternP;
      return this.pattern.equals(patternMotif.pattern);
    }
    else if ((patternP instanceof String)) {
      String str = (String) patternP;
      return this.pattern.equals(str);
    }
    return super.equals(patternP);
  }

  /* (non-Javadoc)
   * @see java.lang.Object#toString()
   */
  @Override
  public int compareTo(PatternMotif patternP) {
    int lengthDelta = Integer.compare(patternP.length(), length());
    int quantytyRelDelta = Double.compare(patternP.getQuantityRefRelative(), getQuantityRefRelative());
    return lengthDelta != 0 ? lengthDelta : quantytyRelDelta;
  }

  public int length(){
    return pattern.length();
  }

  public void addQuantityBen(int posRel) {
    this.positions.add(posRel);
  }

  public double getQuantityRefRelative() {
    return (double) quantityRef / (quantityAlt + 1);
  }

  public double getQuantityBenRelative() {
    if (getQuantityBen() < 2) {
      return 0;
    }
    double relativeQuantity = ((getQuantityBen() - 1.) / cluster.variantSize()) * getRandomOccuranceProbability();
    if (relativeQuantity < 0) {
      System.out.println();
    }
    return relativeQuantity;
  }
  
  /**
   * @return The relevance of the pattern depending on its length
   */
  public double getRandomOccuranceProbability(){
    int lengthPatternMax = Config.lengthIntronPatternMax;
    int lengthPattern = pattern.length();
    int possibleOccurances = lengthPatternMax - (lengthPattern - 1);
    double possibleCombinations = Math.pow(4, lengthPattern);
    return possibleCombinations / possibleOccurances;
  }

  /**
   * @param patternP
   * @return
   */
  public boolean containsOnce(PatternMotif patternP) {
    int pos = pattern.indexOf(patternP.pattern);
    if (pos == -1) {
      return false;
    }
    int posLast = pattern.lastIndexOf(patternP.pattern);
    return pos == posLast;
  }
  
  /**
   * @param patternP
   * @return
   */
  public boolean contains(PatternMotif patternP) {
    return pattern.contains(patternP.pattern);
  }

  /**
   * @param patternP
   * @return
   */
  public boolean contains(String patternP) {
    return pattern.contains(patternP);
  }

  public int getQuantityBen() {
    return positions.size();
  }
  
  public double getPosition() {
    double pos = Functions.mean(positions);
    return pos;
  }

  public void resetQuantityBenPositions() {
    this.positions = new ArrayList<Integer>();
  }

  /* (non-Javadoc)
   * @see java.lang.Object#clone()
   */
  @Override
  public PatternMotif clone() {
    return new PatternMotif(this);
  }

  public static ArrayList<PatternMotif> clone(List<PatternMotif> list) {
    ArrayList<PatternMotif> clone = new ArrayList<PatternMotif>(list.size());
    for(PatternMotif item: list) {
      clone.add(new PatternMotif(item));
    }
    return clone;
  }

  /**
   * @param list
   * @return
   */
  public static HashMap<PatternMotif, PatternMotif> cloneToMap(ArrayList<PatternMotif> list) {
    HashMap<PatternMotif, PatternMotif> clone = new HashMap<PatternMotif, PatternMotif>(list.size());
    for(PatternMotif item: list) {
      clone.put(item, new PatternMotif(item));
    }
    return clone;
  }

  /**
   * TODO if no character matches throw an exception 
   * @param pat1
   * @param pat2
   * @return
   */
  private static int align(PatternMotif pat1, PatternMotif pat2){
    int matchesMax = 0;
    int posMax = Integer.MIN_VALUE;
//    String str1 = pattern1.pattern;
//    String str2 = pattern2.pattern;
    for (int pos = pat1.length()-1; pos > -pat2.length(); pos--) {
      int from1 = pos > 0 ? pos : 0;
      int to1 = pos + pat2.length() < pat1.length() ? pos + pat2.length() : pat1.length();
      int from2 = -pos > 0 ? -pos : 0;
      int to2 = -pos + pat1.length() < pat2.length() ? -pos + pat1.length() : pat2.length();
//      System.out.println(from1 + " > " + to1 + " of " + str1);
//      System.out.println(from2 + " > " + to2 + " of " + str2);
      int matches = Functions.countMatchingChar(pat1.pattern.substring(from1, to1), pat2.pattern.substring(from2, to2));
      //            System.out.println(matches);
      if (matches > 0 && matches > matchesMax) {
        matchesMax = matches;
        posMax = pos;
      } else if (matches > 0 && matches == matchesMax) {
        double pos1 = pat1.getPosition();
        double pos2 = pat2.getPosition();
        double posDif = pos + pos2 - pos1;
        if (posDif < 0) {
          matchesMax = matches;
          posMax = pos;
        } else if (posDif > 0) {
        } else {
          Log.add(posMax + " and " + pos + " are both matching " + matches + " characters.", 2);
        }
      }
    }
    if (Math.abs(posMax) > 9) {
      posMax = posMax;
    }
    return posMax;
  }
  
  /**
   * TODO if no character matches throw an exception 
   * @param str1
   * @param str2
   * @return
   */
  public static int align(String str1, String str2){
    int matchesMax = 0;
    int posMax = Integer.MIN_VALUE;
//    String str1 = pattern1.pattern;
//    String str2 = pattern2.pattern;
    for (int pos = str1.length()-1; pos > -str2.length(); pos--) {
      int from1 = pos > 0 ? pos : 0;
      int to1 = pos + str2.length() < str1.length() ? pos + str2.length() : str1.length();
      int from2 = -pos > 0 ? -pos : 0;
      int to2 = -pos + str1.length() < str2.length() ? -pos + str1.length() : str2.length();
//      System.out.println(from1 + " > " + to1 + " of " + str1);
//      System.out.println(from2 + " > " + to2 + " of " + str2);
      int matches = Functions.countMatchingChar(str1.substring(from1, to1), str2.substring(from2, to2));
      //			System.out.println(matches);
      if (matches > 0 && matches > matchesMax) {
        matchesMax = matches;
        posMax = pos;
      } else if (matches > 0 && matches == matchesMax) {
//        Log.add(posMax + " and " + pos + " are both matching " + matches + " characters.", 2);
      }
    }
    if (Math.abs(posMax) > 9) {
      int x;
    }
    return posMax;
  }

  /**
   * @return
   */
  public ArrayList<Integer> getPositions() {
    return positions;
  }

  /**
   * @return
   */
  public double getQuantityRef() {
    return quantityRef;
  }

  /**
   * @param clusterP
   * @return
   */
  public boolean setCluster(Cluster clusterP) {
    this.cluster = clusterP;
    return cluster != null;
  }
}
