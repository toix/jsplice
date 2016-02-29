/**
 * 
 */
package jsplice.tools;

import java.util.ArrayList;

import jsplice.data.Config;
import jsplice.data.PatternMotif;
import jsplice.io.Variants;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 * TODO differenz der PatternLänge berücksichtigen
 */
public class ClusterCorrelation {

  private ArrayList<Cluster> cluster;
  private Variants variants;
  private int[][] quantityClCl;
  private int[] quantityCl;
  private double[][] correlation;
  private ArrayList<Integer> idx;
  private boolean[][] matchingCluster;
  private boolean[][] notMatchable;

  /**
   * Count only the longest pattern of every cluster that the sequences contains
   * @param variantsP 
   * @param clusterP 
   */
  public ClusterCorrelation(ArrayList<Cluster> clusterP, Variants variantsP) {
    this.cluster = Cluster.cloneSoft(clusterP);
    this.variants = variantsP;
    idx = new ArrayList<Integer>(cluster.size());
    for (int i = 0; i < cluster.size(); i++) {
      idx.add(i);
    }
    notMatchable = new boolean[cluster.size()][cluster.size()];
    this.correlation = createClusterCorrelationMatrix(variants);
  }

  /**
   * TODO slow
   * Count only the longest pattern of every cluster that the sequences contains
   * @param cluster
   * @param variants
   * @return
   */
  private double[][] createClusterCorrelationMatrix(Variants variants) {
    double[][] correlation = new double[cluster.size()][cluster.size()];
    quantityClCl = new int[cluster.size()][cluster.size()];
    quantityCl = new int[cluster.size()];
    int lengthIntronMax = Config.lengthIntronPatternMax;
    matchingCluster = new boolean[cluster.size()][variants.size()];
    for (int c = 0; c < cluster.size(); c++) {
      cluster.get(c).resetQuantityBen();
    }
    for (int v = 0; v < variants.size(); v++) {
      for (int c = 0; c < cluster.size(); c++) {
        Cluster clusterC = cluster.get(c);
        int posChange = variants.get(v).getSequence().getPositionChange();
        int min = posChange - lengthIntronMax + 1;
        int max = posChange + lengthIntronMax - 1;
        String sequenceStr = variants.get(v).getSequence().substring(min, max);
        if (clusterC.isContainedBy(sequenceStr)) {
          matchingCluster[c][v] = true;
          quantityCl[c]++;
          int posRel = variants.get(v).getSequence().getPositionChangeRelative();
          if (!clusterC.addQuantityBen(sequenceStr, posRel)) {
            throw new IllegalArgumentException();
          }
        } else {
          matchingCluster[c][v] = false;
        }
      }
      for (int c1 = 0; c1 < matchingCluster.length; c1++) {
        for (int c2 = 0; c2 < matchingCluster.length; c2++) {
          if (matchingCluster[c1][v] && matchingCluster[c2][v] && c1 != c2) {
            quantityClCl[c1][c2]++;
          }
        }
      }
    }
    for (int c1 = 0; c1 < correlation.length; c1++) {
      for (int c2 = 0; c2 < correlation[c1].length; c2++) {
        int qty1 = quantityCl[c1];
        int qty2 = quantityCl[c2];
        int len1 = cluster.get(c1).getPatternCore().pattern.length();
        int len2 = cluster.get(c2).getPatternCore().pattern.length();
        int lenDif = Math.abs(len1 - len2);
        int lenMin = len1 < len2 ? len1 : len2;
        int qtyMin = (qty1 < qty2 ? qty1 : qty2);
        double divisor = qtyMin / ((3. + lenMin)/(4. + lenMin) + 0.1) + (1. + lenDif) / (2. + lenDif) + (1. + qtyMin) / (2. + qtyMin);
        if (divisor > 0) {
          correlation[c1][c2] = quantityClCl[c1][c2] / divisor;
        } else {
          correlation[c1][c2] = 0;
        }
      }
    }
    return correlation;
  }

  public double[][] refresh(int newIdxP, int delIdxP) {
    int newIdx = idx.get(newIdxP);
    int delIdx = idx.get(delIdxP);
    Cluster newCl = cluster.get(newIdx);
    int lengthIntronMax = Config.lengthIntronPatternMax;
    for (int c = 0; c < matchingCluster.length; c++) {
      quantityClCl[newIdx][c] = 0;
      quantityClCl[c][newIdx] = 0;
    }
    quantityCl[newIdx] = 0;
    newCl.resetQuantityBen();
    for (int v = 0; v < variants.size(); v++) {
      int posChange = variants.get(v).getSequence().getPositionChange();
      int min = posChange - lengthIntronMax + 1;
      int max = posChange + lengthIntronMax - 1;
      String sequenceStr = variants.get(v).getSequence().substring(min, max);
      if (newCl.isContainedBy(sequenceStr)) {
        matchingCluster[newIdx][v] = true;
        quantityCl[newIdx]++;
        int posRel = variants.get(v).getSequence().getPositionChangeRelative();
        if (!newCl.addQuantityBen(sequenceStr, posRel)) {
          throw new IllegalArgumentException();
        }
      } else {
        matchingCluster[newIdx][v] = false;
      }
      for (int c = 0; c < matchingCluster.length; c++) {
        if (matchingCluster[c][v] && matchingCluster[newIdx][v] && c != newIdx) {
          quantityClCl[newIdx][c]++;
          quantityClCl[c][newIdx]++;
        }
      }
    }
    for (int c = 0; c < correlation.length; c++) {
      int qty1 = quantityCl[newIdx];
      int qty2 = quantityCl[c];
      int len1 = cluster.get(newIdx).getPatternCore().pattern.length();
      int len2 = cluster.get(c).getPatternCore().pattern.length();
      int lenDif = Math.abs(len1 - len2);
      int lenMin = len1 < len2 ? len1 : len2;
      int qtyMin = (qty1 < qty2 ? qty1 : qty2);
      double divisor = qtyMin / ((1. + lenMin)/(2. + lenMin)) + (1. + lenDif) / (2. + lenDif) + (1. + qtyMin) / (2. + qtyMin);
      if (divisor > 0) {
        correlation[c][newIdx] = quantityClCl[c][newIdx] / divisor;
        correlation[newIdx][c] = quantityClCl[newIdx][c] / divisor;
      } else {
        correlation[c][newIdx] = 0;
        correlation[newIdx][c] = 0;
      }
    }
    idx.remove(new Integer(delIdx));
    return getCurrentCorrelation();
  }

  /**
   * TODO delete parts
   * @param newIdx
   * @param delIdx
   * @return 
   */
  public double[][] notMatchable(int newIdxP, int delIdxP) {
    int newIdx = idx.get(newIdxP);
    int delIdx = idx.get(delIdxP);
    notMatchable[delIdx][newIdx] = true;
    notMatchable[newIdx][delIdx] = true;
    return getCurrentCorrelation();
  }

  /**
   * @return
   */
  private double[][] getCurrentCorrelation() {
    double[][] correlation = new double[idx.size()][idx.size()];
    for (int c1 = 0; c1 < idx.size(); c1++) {
      for (int c2 = 0; c2 < idx.size(); c2++) {
        int a1 = idx.get(c1);
        int a2 = idx.get(c2);
        if (notMatchable[a1][a2]) {
          correlation[c1][c2] = 0;
        } else {
          correlation[c1][c2] = this.correlation[a1][a2];
        }
      }
    }
    return correlation;
  }

  public double[][] getOriginal() {
    return correlation;
  }

  /**
   * @param clusterP
   * @param variantsP
   * @return
   */
  public void check(ArrayList<Cluster> clusterP, Variants variantsP) {
    double[][] corValid = new ClusterCorrelation(clusterP, variantsP).getOriginal();
    double[][] corTest = getCurrentCorrelation();
    for (int i = 0; i < corTest.length; i++) {
      for (int j = 0; j < corTest[0].length; j++) {
        if (corValid[i][j] != corTest[i][j] && corTest[i][j] != 0) {
          double[] x = corValid[i];
          double[] y = corValid[j];
          double[] a = corTest[i];
          double[] b = corTest[j];
          int ix1 = idx.get(i);
          int ix2 = idx.get(j);
          throw new IllegalArgumentException(i + " " + j);
        }
      }
    }
  }

}
