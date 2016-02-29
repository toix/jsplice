/**
 * 
 */
package jsplice.tools;

import java.time.Clock;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Timer;

import jsplice.data.Config;
import jsplice.data.PatternMotif;
import jsplice.data.RefGene;
import jsplice.data.Sequence;
import jsplice.data.Sequences;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.io.Variants;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import com.google.common.primitives.Doubles;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Model implements Runnable {

  /**
   * Number of DNA Bases (=4)
   */
  public static final int numberOfBases = Functions.bases.length();
  private HashMap<String, Cluster> clusterHash;
  private HashMap<String, Cluster> clusterHashPatho;
  /**
   * The probability of the bases at each position relative to the junction
   */
  private double[][] probability;
  /**
   * Weight matrix to calculate Individual Information by position and base
   */
  private double[][] weightMatrix;
  // /**
  // * Weight matrix to calculate Individual Information by position and base of the changes
  // */
  // private double[][] weightMatrixChange;
  /**
   * Pearson correlation of all locations between each other
   */
  private double[][] correlation;
  /**
   * 
   */
  Sequences sequences;
  /**
   * Constructor without filters
   * 
   * @param variants
   *            {@link Variants} with corresponding {@link RefGene} and sequence
   * @param acceptorP
   *            True -> extract Variants before that are before exon
   */
  public Model(Variants variantsP, boolean acceptorP) {
    if (variantsP.size() < 1) {
      throw new IllegalArgumentException("The parameter contains no variants.");
    }
    this.sequences = new Sequences(variantsP, acceptorP);
    Log.add("Number of pathogene training sequences for " + (acceptorP? "acceptor" : "donor") + " site: " + sequences.getSequences().size(), 3);
  }

  /**
   * Calculates the maximum intron length for the change matrix with sufficient variants per position
   * 
   * @param variants
   * @return The intron length for the change matrix
   */
  private int calculateChangeIntronLength(Variants variants) {
    // Extract all references of the Variants
    ArrayList<String> allRef = Variants.extractAllRef(variants);
    System.out.println(allRef);
    int intronLength = variants.get(0).getSequence().getStringIntron().length();
    int changeIntronLength = intronLength;
    // System.out.println("for: " + (junctionPosition > 0) + " && " + (changeIntronLength == intronLength));
    // System.out.println("junctionP: " + junctionPosition);
    if (isAcceptor()) {
      for (int i = sequences.getJunctionPosition(); i > 0 && changeIntronLength == intronLength; i--) {
        // System.out.println(allRef.get(i-1).length());
        if (allRef.get(i - 1).length() < 12) {
          changeIntronLength = sequences.getJunctionPosition() - i;
        }
      }
    } else {
      for (int i = sequences.getJunctionPosition(); i < sequences.length() - 1 && changeIntronLength == intronLength; i++) {
        if (allRef.get(i + 1).length() < 10) {
          changeIntronLength = i - sequences.getJunctionPosition();
        }
      }
    }
    System.out.println("change intron length: " + changeIntronLength);
    return changeIntronLength;
  }

  /**
   * shorten and set sequences by the given sequence length
   * 
   * @param sequenceLengthNew
   *            The desired length
   * @param variants
   */
  private ArrayList<Sequence> shortenSequences(Variants variants, int sequenceLengthNew) {
    ArrayList<Sequence> shorterSequences = new ArrayList<Sequence>();
    for (int i = 0; i < variants.size(); i++) {
      Variant variant = variants.get(i);
      Sequence sequence = variant.getSequence();
      int sequenceLengthBefore = sequence.length();
      int intronLengthBefore = sequence.getStringIntron().length();
      int intronLengthNew = (int) Functions.round(intronLengthBefore * sequenceLengthNew / sequenceLengthBefore, 0);
      int intronLengthDelta = Math.abs(intronLengthBefore - intronLengthNew);
      int exonLengthBefore = sequence.getStringExon().length();
      int exonLengthNew = sequenceLengthNew - intronLengthNew;
      int exonLengthDelta = Math.abs(exonLengthBefore - exonLengthNew);
      // String sequenceString;
      // // int junctionPositionAfter;
      // // int variantPositionAfter;
      // if (intronExonJunction) {
      // sequenceString = sequence.getIntronicPart().substring(intronLengthDelta) + sequence.getExonicPart().substring(0,
      // exonLengthNew);
      // // junctionPositionAfter = intronLenghtAfter;
      // // variantPositionAfter = variantPositionBefore - intronLengthDelta;
      // } else {
      // sequenceString = sequence.getExonicPart().substring(exonLengthDelta) + sequence.getIntronicPart().substring(0,
      // sequenceLengthNew);
      // // junctionPositionAfter = exonLengthAfter - 1;
      // // variantPositionAfter = variantPositionBefore - exonLengthDelta;
      // }
      try {
        Sequence newSeq = new Sequence(sequence.getStringExtended(), sequence.getLengthExonExtended(), sequenceLengthNew,
            exonLengthNew, variant);
        shorterSequences.add(newSeq);

        // System.out.println("before: " + intronLengthBefore + "\t after: " + intronLenghtAfter);
        // if(Math.random()>0.9){
        // System.out.println("new seq: " + newSeq);
        // System.out.println("old seq: " + sequence);
        // System.out.println(sequence.getIntronicPart());
        // System.out.println(sequence.getExonicPart().substring(0, exonLengthAfter));
        // System.out.println(sequenceString);
        // System.out.println("seq len before: " + sequenceLengthBefore);
        // System.out.println("intron length before: " + intronLengthBefore + "\t intronLength after: " + intronLenghtAfter);
        // System.out.println("intron delta: " + intronLengthDelta + "\t exon delta: " + exonLengthDelta);
        // System.out.println("exon len a: " + exonLengthAfter + "\t intron length a: " + intronLenghtAfter);
        // System.out.println("sequence length: " + sequenceString.length());
        // System.out.println();
        // }
      } catch (IllegalArgumentException e) {
        Log.add("Failed to shorten sequence to " + sequenceLengthNew + ": " + e.getMessage(), 2);
      }
    }
    return shorterSequences;
  }

  /**
   * Calculate probability and weigth matrix
   */
  public double[][] calculateMatrix() {
    // // calculate uncertainty
    // uncertainty = Functions.getDoubleArrayWithZeroes(sequenceLength);
    // for (int l = 0; l < sequenceLength; l++) {
    // for (int b = 0; b < numberOfBases; b++) {
    // if(probability[l][b] != 0){
    // uncertainty[l] += - Math.log(probability[l][b]) / Math.log(2);
    // }
    // }
    // }
    // informationContent = new double[sequenceLength];
    // calculate probability
    probability = Functions.getFrequencies(sequences.getSequences(), true);
    // calculate error and weight matrix
    double[][] weightMatrix = new double[sequences.length()][numberOfBases];
    for (int l = 0; l < sequences.length(); l++) {
      // informationContent[l] = 2 - (uncertainty[l] + error);
      double error = (4.0 - 1) / (2 * Math.log(2) * sequences.size());
      for (int b = 0; b < numberOfBases; b++) {
        if (probability[l][b] == 0) {
          probability[l][b] = 1.0 / (sequences.size() * 2);
        }
        weightMatrix[l][b] = 2.0 - (-Math.log(probability[l][b]) / Math.log(2) + error);
      }
    }
    return weightMatrix;

    // calculate individual information of the natural splice sites
    // individualInformation = new double[sequences.size()][sequences.getSequenceLength()];
    // double[] sequenceEntropy = new double[sequences.size()];
    // for (int i = 0; i < sequences.size(); i++) {
    // individualInformation[i] = getIndividualInformation(sequences.get(i), getJunctionPosition(), true);
    // sequenceEntropy[i] = Functions.sum(individualInformation[i]);
    // }
    // Arrays.sort(sequenceEntropy);
    // System.out.println("First: " + sequenceEntropy[0] + "\tLast: " + sequenceEntropy[sequenceEntropy.length-1]);
  }

  // public void runOld() {
  // // calculate probability
  // probability = Functions.getProbabilities(sequences);
  // // calculate error
  // double error = 1/Math.log(2) * (4-1)/(2*sequences.size());
  // // calculate uncertainty
  // uncertainty = new double[sequences.getSequenceLength()];
  // for (int l = 0; l < sequences.getSequenceLength(); l++) {
  // for (int b = 0; b < numberOfBases; b++) {
  // if(probability[l][b] != 0){
  // uncertainty[l] += - probability[l][b] * Math.log(probability[l][b]) / Math.log(2);
  // }
  // }
  // }
  // // calculate information content and weight matrix
  // informationContent = new double[sequences.getSequenceLength()];
  // weightMatrix = new double[sequences.getSequenceLength()][numberOfBases];
  // for (int l = 0; l < sequences.getSequenceLength(); l++) {
  // informationContent[l] = 2 - (uncertainty[l] + error);
  // for (int b = 0; b < numberOfBases; b++) {
  // weightMatrix[l][b] = probability[l][b] * informationContent[l];
  // }
  // }
  // // calculate individual information of the natural splice sites
  // individualInformation = new double[sequences.size()][sequences.getSequenceLength()];
  // double[] sequenceEntropy = new double[sequences.size()];
  // for (int i = 0; i < sequences.size(); i++) {
  // individualInformation[i] = getIndividualInformation(sequences.get(i), getJunctionPosition(), true);
  // sequenceEntropy[i] = Functions.sum(individualInformation[i]);
  // }
  // Arrays.sort(sequenceEntropy);
  // System.out.println("First: " + sequenceEntropy[0] + "\tLast: " + sequenceEntropy[sequenceEntropy.length-1]);
  //
  // // calculate correlation
  // correlation = new double[sequences.getSequenceLength()][sequences.getSequenceLength()];
  // for (int l1 = 0; l1 < sequences.getSequenceLength(); l1++) {
  // for (int l2 = 0; l2 < sequences.getSequenceLength(); l2++) {
  // correlation[l1][l2] = getCorrelation(l1, l2);
  // }
  // }
  // }

  public void calculateCorrelation() {
    // calculate correlation
    correlation = new double[sequences.length()][sequences.length()];
    for (int l1 = 0; l1 < sequences.length(); l1++) {
      for (int l2 = 0; l2 < sequences.length(); l2++) {
        correlation[l1][l2] = getCorrelation(l1, l2);
      }
    }
  }

  /**
   * @param location1
   *            Location on the sequence
   * @param location2
   *            Location on the sequence
   * @return Average Pearson correlation between the two locations
   */
  public double getCorrelation(int location1, int location2) {
    // System.out.println("size: " + size());
    PearsonsCorrelation pearsons = new PearsonsCorrelation();
    // Generate a binary matrix for the two locations. Every sequence contains one 1 for the occurring base.
    double[][] binaryProbabilityL1 = new double[numberOfBases][sequences.size()];
    double[][] binaryProbabilityL2 = new double[numberOfBases][sequences.size()];
    double[][] correlationBases = new double[numberOfBases][numberOfBases];
    double correlationSum = 0.0;
    // initialize with 0
    for (int b = 0; b < numberOfBases; b++) {
      binaryProbabilityL2[b] = Functions.getInitializedDoubleArray(sequences.size());
      binaryProbabilityL1[b] = Functions.getInitializedDoubleArray(sequences.size());
    }
    for (int j = 0; j < sequences.size(); j++) {
      char base1 = sequences.get(j).charAt(location1);
      char base2 = sequences.get(j).charAt(location2);
      Integer baseIdx1 = Functions.mapNumber.get(base1);
      Integer baseIdx2 = Functions.mapNumber.get(base2);
      binaryProbabilityL1[baseIdx1][j] = 1;
      binaryProbabilityL2[baseIdx2][j] = 1;
    }
    for (int b1 = 0; b1 < numberOfBases; b1++) {
      for (int b2 = 0; b2 < numberOfBases; b2++) {
        // System.out.println(Functions.arrayToString(binaryProbabilityL1[b1]));
        // System.out.println(Functions.arrayToString(binaryProbabilityL2[b2]));
        if (Functions.sum(binaryProbabilityL1[b1]) == 0 || Functions.sum(binaryProbabilityL2[b2]) == 0) {
          correlationBases[b1][b2] = -1;
        } else if (Arrays.equals(binaryProbabilityL1[b1], binaryProbabilityL2[b2])) {
          correlationBases[b1][b2] = 1d;
        } else {
          correlationBases[b1][b2] = pearsons.correlation(binaryProbabilityL1[b1], binaryProbabilityL2[b2]);
        }
        // System.out.println(b1 + "," + b2 + ":\t " + correlationBases[b1][b2]);
      }
    }

    // System.out.println("l1: " + location1 + "\t l2: " + location2);
    for (int j = 0; j < sequences.size(); j++) {
      char base1 = sequences.get(j).charAt(location1);
      char base2 = sequences.get(j).charAt(location2);
      Integer baseIdx1 = Functions.mapNumber.get(base1);
      Integer baseIdx2 = Functions.mapNumber.get(base2);
      correlationSum += correlationBases[baseIdx1][baseIdx2];
      // System.out.println("j: " + j + "\t bx1: " + base1 + "\t bx2: " + base2 + "\t cor: " + correlationBases[baseIdx1][baseIdx2]);
      if (Double.isNaN(correlationBases[baseIdx1][baseIdx2])) {
      }
    }
    // if (location1 == 24 && location2 == 25) {
    // System.out.println("correlation Matrix");
    // for (int i = 0; i < correlationBases.length; i++) {
    // System.out.println(Functions.arrayToString(correlationBases[i]));
    // }
    // int j = 0;
    // char base1 = sequences.get(j).charAt(location1);
    // char base2 = sequences.get(j).charAt(location2);
    // Integer baseIdx1 = Functions.mapNumber.get(base1);
    // Integer baseIdx2 = Functions.mapNumber.get(base2);
    // System.out.println("b1: " + baseIdx1 + "\t b2: " + baseIdx2);
    // System.out.println("result: " + correlationSum/(double)size());
    // }
    return correlationSum / (double) sequences.size();
  }

  /**
   * Calculates the individual information for a sequence with a certain junction position
   * 
   * @param sequenceP
   *            A {@link Sequence}
   * @param junction
   *            The position of the potential splice site on the pattern
   * @param reference
   *            false -> the alternate sequence is used
   * @param limitLength
   *            true -> lengthIntronCryptic in Config defines the maximum intron length
   * @return A summand for each position
   */
  public Result getInformation(Sequence sequenceP, int junction, boolean reference, boolean limitLength) {
    // Log.add("Calculate individual information for the " + (reference ? "reference" : "alternate") + " of "
    // + (sequenceP.isAcceptor() ? "an acceptor" : "a donor") + " sequence with " + (isAcceptor() ? "accptor" : "donor")
    // + " model and " + (limitLength ? "limited" : "unlimited") + " length at position " + junction + ".", 2);
    if (sequenceP.length() < sequences.length()) {
      throw new IllegalArgumentException("The length of the pattern Sequence (" + sequenceP.length()
          + ") has to be equal or bigger than the length of this instance (" + sequences.length() + ").");
    }
    int min = Config.getMinAnalysisPosition(sequenceP.isAcceptor(), isAcceptor());
    int max = Config.getMaxAnalysisPosition(sequenceP.isAcceptor(), isAcceptor());
    if (junction < min || junction > max) {
      throw new IllegalArgumentException("The position (=" + junction + ") has to be a number between " + min + " and " + max);
    }
    int intronLengthMax = Config.getLengthIntronCryptic(isAcceptor());
    int patternStart = junction - sequences.getJunctionPosition();
    int matrixStart = 0;
    int matrixEnd = sequences.length();
    // limit intron length?
    if (limitLength && isAcceptor() && patternStart < junction - intronLengthMax) {
      patternStart = junction - intronLengthMax;
      matrixStart = sequences.getJunctionPosition() - intronLengthMax;
      matrixEnd = matrixStart + Config.getLengthModelExon() + intronLengthMax;
    } else if (limitLength && !isAcceptor() && matrixEnd > Config.getLengthModelExon() + intronLengthMax) {
      matrixEnd = matrixStart + Config.getLengthModelExon() + intronLengthMax;
    }
    // System.out.println("junction: " + junction);
    // System.out.println("getJunctionPosition(): " + sequences.getJunctionPosition());
    // System.out.println("patternStart: " + patternStart);
    // System.out.println("matrixStart: " + matrixStart);
    // System.out.println("matrixEnd: " + matrixEnd);
    double[] individualInformation = Functions.getInitializedDoubleArray(sequences.length());
    int changePos = sequenceP.getPositionChange();
    for (int locM = matrixStart, locP = patternStart; locM < matrixEnd; locM++, locP++) {
      // try {
      int baseNumber = Functions.mapNumber.get(sequenceP.charAt(locP));
      // is alternate sequence and change position?
      if (!reference && locP == changePos) {
        char alt = sequenceP.getAlt().charAt(0);
        baseNumber = Functions.mapNumber.get(alt);
      }
      individualInformation[locM] = weightMatrix[locM][baseNumber];
      // } catch (StringIndexOutOfBoundsException e) {
      // System.out.println(e.getMessage() + "\nl1:" + l1 + "\t l2: " + l2 + "\t patternStart: " + patternStart + "\t seqLen: " +
      // sequenceLength);
      // System.out.println("seq acc: " + sequenceP.isAcceptor() + "\t model acc: " + isAcceptor());
      // System.out.println("exLen: " + sequenceP.lengthExtended() + "\t len: " + sequenceP.length());
      // }
    }
    // Log.add("Information: " + Functions.sum(individualInformation), 2);
    return new Result(individualInformation, junction, sequenceP);
  }

  /**
   * Calculates the individual information for a sequence with the natural junction position
   * 
   * @param sequences
   * @param reference
   *            true -> reverence sequence is used; false -> the alternate sequence is used
   * @param cluster
   * @return sum of Individual Information for every sequence
   */
  public double[] getInformation(Sequences sequences, boolean reference, boolean cluster) {
    double[] indInfo = new double[sequences.size()];
    for (int s = 0; s < sequences.size(); s++) {
      Sequence sequence = sequences.get(s);
      int junction = sequence.getPositionJunction();
      if (cluster) {
        int changeRel = sequence.getPositionChangeRelative();
        indInfo[s] = 0
            + getInformation(sequence, junction, reference, true).getTotalInformation()
            + sequence.getMaxPatternQty(clusterHash, reference)
            - sequence.getMaxPatternQty(clusterHashPatho, reference);
//        Log.add("info: " + getInformation(sequence, junction, reference, true).getTotalInformation());
//        Log.add("+ " + sequence.getMaxPatternQty(clusterHash, reference));
//        Log.add("- " + sequence.getMaxPatternQty(clusterHashPatho, reference));
      } else {
        indInfo[s] = getInformation(sequence, junction, reference, false).getTotalInformation();
      }
    }
    return indInfo;
  }

  public String probabilityToString() {
    System.out.println();
    System.out.println("probability");
    String probStr = 0 + " \t" + Functions.arrayToString(probability[0], 3);
    for (int i = 1; i < probability.length; i++) {
      probStr = probStr + "\n" + i + " \t" + Functions.arrayToString(probability[i], 3);
    }
    return probStr;
  }

  public String correlationToString() {
    System.out.println();
    System.out.println("correlation");
    String corStr = 0 + " \t" + Functions.arrayToString(correlation[0], 4);
    for (int i = 1; i < correlation.length; i++) {
      corStr = corStr + "\n" + i + " \t" + Functions.arrayToString(correlation[i], 4);
    }
    return corStr;
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
   * @param junction
   *            for relative enumeration
   * @return weight matrix as String
   */
  public String matrixToString(String title, int junction) {
    String matrixString = title + "\n" + (0 - junction) + " \t" + Functions.arrayToString(weightMatrix[0], 5);
    for (int i = 1; i < weightMatrix.length; i++) {
      if (i == junction) {
        junction--;
      }
      matrixString = matrixString + "\n" + (i - junction) + " \t" + Functions.arrayToString(weightMatrix[i], 5);
    }
    return matrixString;
  }

  public boolean isAcceptor() {
    return sequences.isAcceptor();
  }

  /**
   * @return
   */
  public int getJunctionPosition() {
    return sequences.getJunctionPosition();
  }

  public Sequences getSequences() {
    return sequences;
  }

  /**
   * Find all containing the Variant position and cluster them. 
   * TODO sub of merged pattern to non sub important pattern when containing it -> -1% 
   * TODO limit cluster merging OR subtract/check quantity of the important cluster from the longer cluster -> ok
   * TODO add pattern only when share > 0.5
   * TODO penalty with cluster out of conditional (cryptic) -> +1% ; much more stable
   * TODO separate cluster for splice site 
   * TODO respect average/median position in cluster
   * TODO respect sub-pattern in countClusterInSequences()
   * TODO cluster merging only for cluster that are one base longer; recursively
   * 
   * important sequence as non-sub-pattern?
   * 
   * @param variantsP
   * @param acceptorP
   */
  static HashMap<String, Cluster> findPattern(Variants variantsP, boolean acceptorP, boolean benign) {
    int lengthIntronMax = Config.lengthIntronPatternMax;
    int distanceJunctionMin = Config.getDistanceClusterMin(acceptorP);
    int distanceJunctionMax = Config.getDistanceClusterMax();
    variantsP = Filter.filterVariantType(variantsP, acceptorP);
    variantsP = Filter.extractVariantsInRelativeRange(variantsP, distanceJunctionMin, distanceJunctionMax);
    variantsP = Filter.deleteDuplicateJunctions(variantsP);
    HashMap<String, Integer> quantityRef = new HashMap<String, Integer>();
    HashMap<String, Integer> quantityAlt = new HashMap<String, Integer>();
    int numOfPattern[] = Functions.getInitializedIntArray(lengthIntronMax + 1);
    for (int i = 0; i < variantsP.size(); i++) {
      Sequence sequence = variantsP.get(i).getSequence();
      int posChangeAbs = sequence.getPositionChange();
      // count pattern on the variant position and in other sequences at the same position
      for (int length = 1; length <= lengthIntronMax; length++) {
        for (int shift = 0; shift < length; shift++) {
          int from = posChangeAbs + shift - length + 1;
          int to = posChangeAbs + shift;
          int junction = sequence.getPositionJunction();
          if (Math.abs(from - junction) > distanceJunctionMin-1 && Math.abs(to - junction) > distanceJunctionMin-1) {
            String patternRef = sequence.substring(from, to + 1, benign);
            if (!quantityRef.containsKey(patternRef)) {
              quantityRef.put(patternRef, 1);
            } else {
              quantityRef.put(patternRef, quantityRef.get(patternRef) + 1);
            }
            // condition
            String patternAlt = sequence.substring(from, to + 1, !benign);
            if (!quantityAlt.containsKey(patternAlt)) {
              quantityAlt.put(patternAlt, 1);
            } else {
              quantityAlt.put(patternAlt, quantityAlt.get(patternAlt) + 1);
            }
            numOfPattern[length]++;
          }
        }
      }
    }

//    Log.add("Ref: " + quantityRef, 2);
//    Log.add("Alt: " + quantityAlt, 2);
    ArrayList<PatternMotif> pattern = Model.createPattern(quantityRef, quantityAlt);
    //		Log.add("pattern: " + pattern, 2);
    ArrayList<Cluster> cluster = createClusterSub(pattern);
    //		Log.add("cluster: " + cluster, 2);
      if (Config.clusteringBasic) {
        cluster = createCluster(cluster);
//       Log.add("core: " + cluster, 2);
      }
      if (!Config.simpleMerging) {
      cluster = mergeCluster(cluster, variantsP);
//      Log.add("merged: " + cluster, 2);
    }
    cluster = countClusterInSequences(cluster, variantsP);
//    Log.add("Ben count: " + cluster, 2);
    HashMap<String, Cluster> clusterHash = createHash(cluster);
    Log.add("Created " + cluster.size() + " cluster with " + variantsP.size() + " variants for " + (benign? "benign " : "pathogene ") + (acceptorP? "acceptor" : "donor") + " sequences", 3);
    Log.add(cluster, 2);
    return clusterHash;
  }

  /**
   * Create pattern containing the quantities
   * 
   * @param quantityAbsolute
   * @param quantityCondition
   * @return
   */
  private static ArrayList<PatternMotif> createPattern(HashMap<String, Integer> quantityAbsolute, HashMap<String, Integer> quantityCondition) {
    ArrayList<PatternMotif> pattern = new ArrayList<PatternMotif>(quantityAbsolute.size());
    Iterator<Entry<String, Integer>> patternIt = quantityAbsolute.entrySet().iterator();
    while (patternIt.hasNext()) {
      Entry<String, Integer> entry = patternIt.next();
      String patternKey = entry.getKey();
      int quantityAbs = entry.getValue();
      int quantityCon = 0;
      if (quantityCondition.containsKey(patternKey)) {
        quantityCon = quantityCondition.get(patternKey);
      }
      PatternMotif patternNew = new PatternMotif(patternKey, quantityAbs, quantityCon, 0);
      pattern.add(patternNew);
    }
    return pattern;
  }

  /**
   * TODO share prüfen
   * Create cluster for the best pattern and add all sub-pattern containing the cluster pattern
   * @param patternCopy
   * @return
   */
  private static ArrayList<Cluster> createClusterSub(ArrayList<PatternMotif> patternP) {
    double limit = Config.quantityRelLimit;
    ArrayList<PatternMotif> patternCopy = PatternMotif.clone(patternP);
    HashMap<PatternMotif, PatternMotif> patternCopyHash = PatternMotif.cloneToMap(patternCopy);
    ArrayList<Cluster> cluster = new ArrayList<Cluster>();
    PatternMotif patternHightest = findHighestPattern(patternCopy);
    while (patternHightest.getQuantityRelative() > limit) {
      PatternMotif patternMain = patternCopyHash.get(patternHightest);
      Cluster clusterNew = new Cluster(patternMain);
      cluster.add(clusterNew);
      patternCopy.remove(patternHightest);
      if (Config.clusteringSub) {
        for (int p = 0; p < patternCopy.size();) {
          PatternMotif patternCurrent = patternCopy.get(p);
          int mainAbs = patternMain.quantityRef;
          int subAbs = patternCurrent.quantityRef;
          if (patternMain.contains(patternCurrent) && 0.5 * subAbs <= mainAbs) {
            // add the unchanged copy of Pattern
            clusterNew.addSub(patternCopyHash.get(patternCurrent));
            patternCopy.remove(patternCurrent);
            //        } else if (0.5 * subAbs <= mainAbs){
            //          patternCurrent.quantityRef -= patternHightest.quantityRef;
            //          p++;
          } else {
            p++;
          }
        }
      }
      patternHightest = findHighestPattern(patternCopy);
    }
    return cluster;
  }

  /**
   * @param patternP
   * @return
   */
  private static PatternMotif findHighestPattern(ArrayList<PatternMotif> patternP) {
    double valueMax = -1000000000;
    PatternMotif patternMax = null;
    // crate a cluster for relevant pattern and remove short ones
    for (PatternMotif patternCurrent : patternP) {
      // Log.add(patternMax + " < " + patternCurrent.getQuantityRelative() + "\t " + (valueMax <
      // patternCurrent.getQuantityRelative()));
      if (valueMax < patternCurrent.getQuantityRelative()) {
        valueMax = patternCurrent.getQuantityRelative();
        patternMax = patternCurrent;
        // Log.add(patternMax);
      }
    }
    return patternMax;
  }
  
  private static ArrayList<Cluster> createClusterSimple(ArrayList<Cluster> clusterP) {
    ArrayList<Cluster> clusterCopy = Cluster.clone(clusterP);
    HashMap<Cluster, Cluster> clusterCopyMap = Cluster.cloneToMap(clusterCopy);
    for (int m = 0; m < clusterCopy.size(); m++) {
        Collections.sort(clusterCopy);
        Cluster patternMain = clusterCopyMap.get(clusterCopy.get(m));
        clusterCopy.set(m, patternMain);
        for (int c = m + 1; c < clusterCopy.size();) {
            Cluster clusterCurrent = clusterCopyMap.get(clusterCopy.get(c));
            if (clusterCurrent.getPatternCore().contains(patternMain.getPatternCore())) {
                clusterCopy.get(m).add(clusterCurrent);
                int mainAbs = clusterCopy.get(m).getPatternCore().quantityRef;
                int subAbs = clusterCopy.get(c).getPatternCore().quantityRef;
                double subRel = clusterCopy.get(m).getPatternCore().getQuantityRelative();
                double limit = Config.quantityRelLimit;
                int lenDif = clusterCopy.get(c).getPatternCore().pattern.length() - clusterCopy.get(m).getPatternCore().pattern.length();
                if (lenDif * subAbs <= mainAbs && subRel > limit) {
                  clusterCopy.remove(c);
                } else {
                  clusterCopy.get(c).getPatternCore().quantityRef -= mainAbs;
                  c++;
                }
            } else {
                c++;
            }
        }
    }
    return clusterCopy;
}

  /**
   * add all cluster as subset of the more important cluster <br>
   * if subset cluster contains the substring of the important cluster
   * 
   * @param clusterCopy
   * @return
   */
  private static ArrayList<Cluster> createCluster(ArrayList<Cluster> clusterP) {
    ArrayList<Cluster> clusterCopy = Cluster.clone(clusterP);
    for (int len = Config.lengthIntronPatternMax - 1; len > 0; len--) {
      for (int m = 0; m < clusterCopy.size(); m++) {
        //				Collections.sort(clusterCopy);
        Cluster clusterMain = clusterCopy.get(m);
        for (int c = m + 1; c < clusterCopy.size();) {
          Cluster clusterCurrent = clusterCopy.get(c);
          boolean contains = clusterCurrent.getPatternCore().contains(clusterMain.getPatternCore());
          int lenMain = clusterMain.getPatternCore().pattern.length();
          int lenSub = clusterCurrent.getPatternCore().pattern.length();
          //          boolean lengthOk = lenMain == len && (lenSub == len + 1 || lenSub == len + 2);
          boolean lengthOk = lenMain == len && lenSub > len;
          int lenDif = lenSub - lenMain;
          int absMain = clusterMain.getPatternCore().quantityRef;
          int absClus = clusterCurrent.getPatternCore().quantityRef;
          boolean quantityOk = (Math.pow(4, lenDif) - 2) * absClus < absMain;
          //          boolean quantityOk = lenDif * absSub < absMain;
          double subRel = clusterCopy.get(m).getPatternCore().getQuantityRelative();
          double limit = Config.quantityRelLimit;
          if (contains && lengthOk && quantityOk) {
            clusterMain.add(clusterCurrent);
            clusterCopy.remove(clusterCurrent);
          } else {
            c++;
          }
        }
      }
    }
    return clusterCopy;
  }

  @SuppressWarnings("unused")
  private static ArrayList<Cluster> mergeCluster(ArrayList<Cluster> clusterP, Variants variants) {
    double limit = Config.mergingCorrelationMin;
    ArrayList<Cluster> cluster = Cluster.clone(clusterP);
    // precalculate first cluster pair
    ClusterCorrelation correlation = new ClusterCorrelation(cluster, variants);
    double[][] correlationCluster = correlation.getOriginal();
    double max = Double.MIN_VALUE;
    int maxPos1 = -1;
    int maxPos2 = -1;
    for (int c = 0; c < correlationCluster.length; c++) {
      double maxTemp = Doubles.max(correlationCluster[c]);
      int maxPosTemp = Doubles.indexOf(correlationCluster[c], maxTemp);
      if (maxTemp > max) {
        max = maxTemp;
        maxPos1 = c;
        maxPos2 = maxPosTemp;
      }
    }
    // merge clusters while limit is not reached
    for (int i = 0; max > limit; i++) {
      // merge
      int newIdx;
      int delIdx;
      if (cluster.get(maxPos1).getPatternCore().quantityRef > cluster.get(maxPos2).getPatternCore().quantityRef) {
        newIdx = maxPos1;
        delIdx = maxPos2;
      } else {
        newIdx = maxPos2;
        delIdx = maxPos1;
      }
      Cluster newCl = cluster.get(newIdx);
      Cluster delCl = cluster.get(delIdx);
      Log.add("Merging\t cl1: " + newCl.getPatternCore().pattern + "\t cl2: " + delCl.getPatternCore().pattern + "\t with correlation " + max, 2);
      boolean success = newCl.add(delCl);
      if (success) {
        cluster.remove(delIdx);
        correlationCluster = correlation.refresh(newIdx, delIdx);
      } else {
        correlationCluster = correlation.notMatchable(newIdx, delIdx);
        Log.add("Clusters can't be aligned:\t " + newCl.getPatternCore().pattern + " and " + delCl.getPatternCore().pattern, 4);
      }
      // reset values
      max = Double.MIN_VALUE;
      maxPos1 = -1;
      maxPos2 = -1;
      // find new pair
      for (int c = 0; c < correlationCluster.length; c++) {
        double maxTemp = Doubles.max(correlationCluster[c]);
        int maxPosTemp = Doubles.indexOf(correlationCluster[c], maxTemp);
        if (maxTemp > max) {
          max = maxTemp;
          maxPos1 = c;
          maxPos2 = maxPosTemp;
        }
      }
    }
    return cluster;
  }

  /**
   * Count only the longest pattern of every cluster that the sequences contain
   * 
   * @param cluster
   * @param variants
   * @return
   */
  private static ArrayList<Cluster> countClusterInSequences(ArrayList<Cluster> cluster, Variants variants) {
    int lengthIntronMax = Config.lengthIntronPatternMax;
    for (int c = 0; c < cluster.size(); c++) {
      cluster.get(c).sortPattern();
      PatternMotif patternCluster = cluster.get(c).getPatternCore();
      for (int v = 0; v < variants.size(); v++) {
        int posChange = variants.get(v).getSequence().getPositionChange();
        int min = posChange - lengthIntronMax + 1;
        int max = posChange + lengthIntronMax - 1;
        String sequenceStr = variants.get(v).getSequence().substring(min, max);
        if (cluster.get(c).isContainedBy(sequenceStr)) {
          int posRel = variants.get(v).getSequence().getPositionChangeRelative();
          if (!cluster.get(c).addQuantityBen(sequenceStr, posRel)) {
            throw new IllegalArgumentException();
          }
        }
        String sequenceAlt = variants.get(v).getSequence().substring(min, max, false);
        if (sequenceAlt.contains(patternCluster.pattern)) {
          if (!cluster.get(c).addQuantityPat(sequenceAlt)) {
            throw new IllegalArgumentException();
          }
        }
      }
    }
    return cluster;
  }

  /**
   * Create hash to cluster by its pattern
   * 
   * @param cluster
   * @return
   */
  private static HashMap<String, Cluster> createHash(ArrayList<Cluster> cluster) {
    Collections.sort(cluster);
    HashMap<String, Cluster> clusterHash = new HashMap<String, Cluster>();
    for (int c = 0; c < cluster.size(); c++) {
      for (int p = 0; p < cluster.get(c).size(); p++) {
        String key = cluster.get(c).getPattern(p).pattern;
        clusterHash.put(key, cluster.get(c));
      }
    }
    return clusterHash;
  }

  /* (non-Javadoc)
   * @see java.lang.Runnable#run()
   */
  @Override
  public void run() {
    weightMatrix = calculateMatrix();
    clusterHash = findPattern(sequences.getVariants(), isAcceptor(), true);
    clusterHashPatho = findPattern(sequences.getVariants(), isAcceptor(), false);
  }
}
