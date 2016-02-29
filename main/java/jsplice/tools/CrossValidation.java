/**
 * 
 */
package jsplice.tools;

import jsplice.data.Sequence;
import jsplice.data.Sequences;
import jsplice.data.Variant;
import jsplice.data.VariantsTupel;
import jsplice.exception.Log;
import jsplice.io.Variants;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class CrossValidation implements Runnable{

  private Variants variantsPathogene;
  private Variants variantsBenign;
  private boolean acceptor;
  private String folder;
  private int i;

  /**
   * 
   */
  public CrossValidation(Variants variantsPathogene, Variants variantsBenign, boolean acceptorP, String folder, int i) {
    this.variantsPathogene = variantsPathogene;
    this.variantsBenign = variantsBenign;
    this.acceptor = acceptorP;
    this.folder = folder;
    this.i = i;
  }

  /* (non-Javadoc)
   * @see java.lang.Runnable#run()
   */
  @Override
  public void run() {
    Log.add("\n - - - - - - - - - - - - - - - - - - - - - ", 3);
    Log.add("cross validation run nr. " + i, 3);
    VariantsTupel varBenig = VariantsTupel.createTestPool(variantsBenign);
    VariantsTupel varPatho = VariantsTupel.createTestPool(variantsPathogene);
    
    Log.add(" - - - pathogene standard model - - - ", 3);
    Model modelStd = new Model(varPatho.train, acceptor);
    modelStd.run();
    Functions.writeToFile(modelStd.matrixToString("A\tC\tG\tT", modelStd.getJunctionPosition()), folder + "accMatrix.tsv");
//  Model modelStdOther = new Model(varPatho.train, !acceptorP);
//  Functions.writeToFile(modelStd.matrixToString("A\tC\tG\tT", modelStdOther.getJunctionPosition()), folder + "donMatrix.tsv");
    Log.add(" - - - pathogene test variants - - - ", 3);
    Sequences sequencesPathoTest = new Sequences(varPatho.test, acceptor);
    Log.add("Number of pathogene test sequences for " + (acceptor? "acceptor" : "donor") + " site: " + sequencesPathoTest.size(), 3);
    Log.add(" - - - benign test variants - - - ", 3);
    Sequences sequencesBenTest = new Sequences(varBenig.test, acceptor);
    Log.add("Number of benign test sequences for " + (acceptor? "acceptor" : "donor") + " site: " + sequencesBenTest.size(), 3);
    
    // write benign results to file
    String standardResultsBenign = getResultsTable(sequencesBenTest, modelStd, false);
    Functions.writeToFile(standardResultsBenign, folder + "benAccStd" + i + ".tsv");
    String changeResultsBeingn = getResultsTable(sequencesBenTest, modelStd, true);
    Functions.writeToFile(changeResultsBeingn, folder + "benAccChg" + i + ".tsv");
    // write pathogene results to file
    String standardResultsPathogene = getResultsTable(sequencesPathoTest, modelStd, false);
    Functions.writeToFile(standardResultsPathogene, folder + "pathoAccStd" + i + ".tsv");
    String changeResultsPathogene = getResultsTable(sequencesPathoTest, modelStd, true);
    Functions.writeToFile(changeResultsPathogene, folder + "pathoAccChg" + i + ".tsv");
    Log.writeToFile();
  }

  /**
   * @param testSequences
   * @param cluster 
   * @return
   */
  public static String getResultsTable(Sequences testSequences,  Model modelStdAcc, boolean cluster) {
  	double[] resultsReference = modelStdAcc.getInformation(testSequences, true, cluster);
  	double[] resultsAlternate = modelStdAcc.getInformation(testSequences, false, cluster);
  	String results = "";
  	String line = "pos\t r\t a\t ref\t  alt\t start\t sequence";
  	results += line;
  	for (int i = 0; i < testSequences.size(); i++) {
  		Sequence sequence = testSequences.get(i);
  		Variant variant = sequence.getVariant();
  		String ref = variant.getRef();
  		String alt = variant.getAlt();
  		int pos = sequence.getPositionChangeRelative();
  		line = "\n" + pos + "\t " + ref + "\t " + alt + "\t " + resultsReference[i] + "\t " + resultsAlternate[i] + "\t " + variant.getStart() + "\t " + sequence.getStringExtended();
  		results += line;
  	}
  	return results;
  }

}
