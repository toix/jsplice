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
  private double[][] quantityClCl;
  private double[] quantityCl;
  private double[][] correlation;
  private ArrayList<Integer> idx;
  private double[][] matchingCluster;
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
   * Count only the longest pattern of every cluster that the sequences contains
   * @param cluster
   * @param variants
   * @return
   */
  private double[][] createClusterCorrelationMatrix(Variants variants) {
    double[][] correlation = new double[cluster.size()][cluster.size()];
    quantityClCl = new double[cluster.size()][cluster.size()];
    quantityCl = new double[cluster.size()];
    int lengthIntronMax = Config.lengthIntronPatternMax;
    matchingCluster = new double[cluster.size()][variants.size()];
    for (int c = 0; c < cluster.size(); c++) {
      Cluster cl = cluster.get(c);
      cl.sortPattern(PatternComp.desc(AttrPat.LENGTH), AttrPat.QTY_BEN_REL); // TODO ben is not there
      cl.resetQuantityBen();
    }
    for (int v = 0; v < variants.size(); v++) {
      for (int c = 0; c < cluster.size(); c++) {
        Cluster clusterC = cluster.get(c);
        int posChange = variants.get(v).getSequence().getPositionChange();
        int min = posChange - lengthIntronMax + 1;
        int max = posChange + lengthIntronMax - 1;
        String sequenceStr = variants.get(v).getSequence().substring(min, max);
        double rating = clusterC.getRatingMatch(sequenceStr);
        if (rating > 0) {
          matchingCluster[c][v] = rating;
          quantityCl[c] += rating;
          int posRel = variants.get(v).getSequence().getPositionChangeRelative();
          if (!clusterC.addQuantityBen(sequenceStr, posRel)) {
            throw new IllegalArgumentException();
          }
        } else {
          matchingCluster[c][v] = 0;
        }
      }
      for (int c1 = 0; c1 < matchingCluster.length; c1++) {
        for (int c2 = 0; c2 < matchingCluster.length; c2++) {
          if (matchingCluster[c2][v] > 0 && matchingCluster[c1][v] > 0 && c2 != c1) {
            quantityClCl[c1][c2] += Math.sqrt(matchingCluster[c1][v] * matchingCluster[c2][v]);
          }
        }
      }
    }
    for (int c1 = 0; c1 < correlation.length; c1++) {
      for (int c2 = 0; c2 < correlation[c1].length; c2++) {
        Cluster cl1 = cluster.get(c1);
        Cluster cl2 = cluster.get(c2);
        correlation[c1][c2] = getCorrelation(c1, c2);
        
//        Cluster cl1 = cluster.get(c1);
//        Cluster cl2 = cluster.get(c2);
//        int matchingCharsCore = Functions.countMatchingChar(cl1.getPatternCore().pattern, cl2.getPatternCore().pattern);
//        int qty1 = quantityCl[c1];
//        int qty2 = quantityCl[c2];
//        int len1 = cl1.getPatternCore().pattern.length();
//        int len2 = cl2.getPatternCore().pattern.length();
//        int lenDif = Math.abs(len1 - len2);
//        int lenMin = len1 < len2 ? len1 : len2;
//        int qtyMin = qty1 < qty2 ? qty1 : qty2;
//        int qtyMax = qty1 > qty2 ? qty1 : qty2;
//        double divisor = (qtyMin + 0.1 * qtyMax) /1.1 / ((3. + lenMin)/(4. + lenMin) + 0.1) + (1. + lenDif) / (2. + lenDif) + (1. + qtyMin) / (2. + qtyMin);
//        if (matchingCharsCore > 0 && divisor > 0) {
//          correlation[c1][c2] = quantityClCl[c1][c2] / divisor;
//        } else {
//          correlation[c1][c2] = 0;
//        }
      }
    }
    return correlation;
  }

  public double[][] refresh(int newIdxP, int delIdxP) {
    int newIdx = idx.get(newIdxP);
    int delIdx = idx.get(delIdxP);
    Cluster clusterNew = cluster.get(newIdx);
    int lengthIntronMax = Config.lengthIntronPatternMax;
    for (int c = 0; c < matchingCluster.length; c++) {
      quantityClCl[newIdx][c] = 0;
      quantityClCl[c][newIdx] = 0;
      cluster.get(c).sortPattern(PatternComp.desc(AttrPat.LENGTH), AttrPat.QTY_BEN_REL);
    }
    quantityCl[newIdx] = 0;
    clusterNew.resetQuantityBen();
    for (int v = 0; v < variants.size(); v++) {
      int posChange = variants.get(v).getSequence().getPositionChange();
      int min = posChange - lengthIntronMax + 1;
      int max = posChange + lengthIntronMax - 1;
      String sequenceStr = variants.get(v).getSequence().substring(min, max);
      double rating = clusterNew.getRatingMatch(sequenceStr);
      if (rating > 0) {
        matchingCluster[newIdx][v] = rating;
        quantityCl[newIdx] += rating;
        int posRel = variants.get(v).getSequence().getPositionChangeRelative();
        if (!clusterNew.addQuantityBen(sequenceStr, posRel)) {
          throw new IllegalArgumentException();
        }
      } else {
        matchingCluster[newIdx][v] = 0;
      }
      for (int c = 0; c < matchingCluster.length; c++) {
        if (matchingCluster[c][v] > 0 && matchingCluster[newIdx][v] > 0 && c != newIdx) {
          quantityClCl[newIdx][c] += Math.sqrt(matchingCluster[newIdx][v] * matchingCluster[c][v]);
          quantityClCl[c][newIdx] += Math.sqrt(matchingCluster[newIdx][v] * matchingCluster[c][v]);
        }
      }
    }
    for (int c = 0; c < correlation.length; c++) {
      correlation[c][newIdx] = getCorrelation(newIdx, c);
      correlation[newIdx][c] = correlation[c][newIdx];
    }
    idx.remove(new Integer(delIdx));
    return getCurrentCorrelation();
  }

  /**
     * @param c1
     * @param c2
     * @return
     */
    private double getCorrelation(int c1, int c2) {
      Cluster cl1 = cluster.get(c1);
      Cluster cl2 = cluster.get(c2);
  //    int lenMin = len1 < len2 ? len1 : len2;
  //    double qtyMax = qty1 > qty2 ? qty1 : qty2;
      int len1 = cl1.getPatternCore().pattern.length();
      int len2 = cl2.getPatternCore().pattern.length();
      int lenMin = len1 > len2 ? len1 : len2;
      double dividend = quantityClCl[c1][c2];
      double qty1 = quantityCl[c1];
      double qty2 = quantityCl[c2];
      double qtyMin = qty1 < qty2 ? qty1 : qty2;
      double divisor = qtyMin * (1.8 - (1. + lenMin)/(3. + lenMin));
      if (divisor > 0 && dividend > 0) {
        return dividend / divisor;
      } else {
        return 0;
      }
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
  public double[][] getCurrentCorrelation() {
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
