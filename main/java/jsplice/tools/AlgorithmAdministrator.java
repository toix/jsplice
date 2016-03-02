package jsplice.tools;

import jsplice.data.Config;
import jsplice.data.Sequence;
import jsplice.data.Sequences;
import jsplice.data.Variant;
import jsplice.data.VariantsTupel;
import jsplice.exception.Log;
import jsplice.io.Variants;

/**
 * TODO print all variants to analyze them statistically (position, donor vs acceptor)
 * In AlgorithmAdministrator the different algorithms are applied on the data
*/
public class AlgorithmAdministrator {
	
	public AlgorithmAdministrator(Variants variantsPathogene, Variants variantsBenign) {
		
		Log.add("#benign Variants: " + variantsBenign.size());
		Log.add("#pathogene Variants: " + variantsPathogene.size());
		boolean acceptor = true;
		
		String folder = "results/";
		printAllVariants(variantsPathogene, folder);
		
//		Model modelStdAcc = new Model(variantsPathogene, acceptor);
//		Model modelStdDon = new Model(variantsPathogene, !acceptor);
		
		// Filter
//		Variants variantsPathogeneNonCry = Filter.extractCrypticVariants(variantsPathogene, modelStdAcc, modelStdDon, false);
//		variantsPathogene = VariantFile.concat(Filter.filterActivatingVariants(variantsPathogene, modelStdAcc, true), Filter.filterActivatingVariants(variantsPathogene, modelStdDon, true));
		
		folder = Config.folder;
		crossValidate(variantsPathogene, variantsBenign, acceptor, folder);
		
		
		// analyze certain sequence
//		Log.add("\n---- certain sequence ---", 3);
//		Sequence seqTest = findSequencebyPosition(new Sequences(variantsPathogene, acceptor), 66905841);
//		Log.add(seqTest.getVariant()+"", 3);
//		Log.add(seqTest.getStringExtended(), 2);
//		int crypticPos = 0;//seqTest.getPositionChange()+1;
//		double natRef = modelStdAcc.getInformation(seqTest, seqTest.getPositionJunction(), true, false).getTotalInformation();
//		double natAlt = modelStdAcc.getInformation(seqTest, seqTest.getPositionJunction(), false, false).getTotalInformation();
//		double kryRef = modelStdAcc.getInformation(seqTest, crypticPos, true, false).getTotalInformation();
//		double kryAlt = modelStdAcc.getInformation(seqTest, crypticPos, false, false).getTotalInformation();
//		Log.add("nat: " + natRef + ">" + natAlt + "\t kry: " + kryRef + ">" + kryAlt + " at " + (crypticPos), 2);
//		Variants variantsTest = new Variants();
//		variantsTest.add(seqTest.getVariant());
//		Filter.extractCrypticVariants(variantsTest, modelStdAcc, modelStdDon, false);
		
		// separate training and test variants
//		ArrayList<VariantFile> separatedVariants;
//		separatedVariants = separateData(variantsBenign);
//		variantsBenTrain = separatedVariants.get(0);
//		variantsBenTest = separatedVariants.get(1);
//		separatedVariants = separateData(variantsPathogene);
//		variantsPathoTrain = separatedVariants.get(0);
//		variantsPathoTest = separatedVariants.get(1);
		
//		System.out.println("\n\n - - - - - - - - - - - - - - - - - - - - - \n");
//		String folder = "results/nonCryptic/";
//		// standard model pathogene
//		Log.add("\n - - - pathogene standard model - - - ", 3);
//		modelStdAcc = new Sequences(variantsPathoTrain, acceptor);
//		Log.add("Number of pathogene training sequences for acceptor site: " + modelStdAcc.size(), 3);
//		Functions.writeToFile(modelStdAcc.matrixToString("standard matrix", modelStdAcc.getJunctionPosition()),
//				folder + "standardMatrix.tsv");
//		modelStdDon = new Sequences(variantsPathoTrain, !acceptor);
//		Log.add("Number of pathogene training sequences for donor site: " + modelStdDon.size(), 3);
//		// change model pathogene
//		Log.add("\n - - - pathogene change model - - - ", 3);
//		Sequences modelChgAcc = new Sequences(variantsPathoTrain, modelStdAcc, modelStdDon, acceptor);
//		Functions.writeToFile(modelChgAcc.matrixToString("change matrix", modelChgAcc.getJunctionPosition()), folder + "changeMatrix.tsv");
//		// create pathogene test sequences
//		Log.add("\n - - - pathogene test variants - - - ", 3);
//		// variantsPathoTest = VariantFile.extractCrypticVariants(variantsPathoTest, modelStdAcc, modelStdDon, false);
//		Sequences sequencesPathoTestAcc = new Sequences(variantsPathoTest, acceptor);
//		Log.add("Number of pathogene test sequences for acceptor site: " + sequencesPathoTestAcc.size(), 3);
//		// create benign test sequences
//		Log.add("\n - - - benign test variants - - - ", 3);
//		Sequences sequencesBenTestAcc = new Sequences(variantsBenTest, acceptor);
		
		
//		// Algo differences
//		double cutoffStd = Config.getLengthModelIntron()/3.470481 - 1.79548;
//		double cutoffChg = Config.getLengthModelIntron()/4.694624 - 2.14212;
//		for (int i = 0; i < variantsPathoTest.size(); i++) {
//			Sequence sequence = variantsPathoTest.get(i).getSequence();
//			double natRefStd = Functions.sum(modelStdAcc.getIndividualInformation(sequence, sequence.getPositionJunction(), true));
//			double natAltStd = Functions.sum(modelStdAcc.getIndividualInformation(sequence, sequence.getPositionJunction(), false));
//			double natRefChg = Functions.sum(modelChgAcc.getIndividualInformation(sequence, sequence.getPositionJunction(), true));
//			double natAltChg = Functions.sum(modelChgAcc.getIndividualInformation(sequence, sequence.getPositionJunction(), false));
//			if (natRefStd > cutoffStd != natRefChg > cutoffChg) {
//				System.out.println("Std: " + natRefStd + " < " + cutoffStd + "\t Chg: " + natRefChg + " < " + cutoffChg);
//				System.out.println(sequence.getVariant());
//			}
//			if (natAltStd < cutoffStd != natAltChg < cutoffChg) {
//				System.out.println("Std: " + natAltStd + " > " + cutoffStd + "\t Chg: " + natAltChg + " > " + cutoffChg);
//				System.out.println(sequence.getVariant());
//			}
//		}
		
		
		// alle weiter weg kryptisch?
//		System.out.println("\n---- cryptic variants ---");
//		VariantFile variantsNonCrypticAcc = VariantFile.extractCrypticVariants(VariantFile.extractVariantsInRange(variantsPathoTest, -100, -6), modelStdAcc, modelStdDon, true);
//		VariantFile VariantsCrypticAcc = VariantFile.extractCrypticVariants(VariantFile.extractVariantsInRange(variantsPathoTest, -100, -6), modelStdAcc, modelStdDon, false);
//		Sequences sequencesNonCryptic = new Sequences(variantsNonCrypticAcc, acceptor);
//		Sequences sequencesCryptic = new Sequences(VariantsCrypticAcc, acceptor);
//		
//		String crypticStd = getResultsTable(sequencesCryptic, modelStdAcc);
//		Functions.writeToFile(crypticStd, "results/crypticStd.tsv");
//		
//		String crypticCh = getResultsTable(sequencesCryptic, modelChAcc);
//		Functions.writeToFile(crypticCh, "results/crypticCh.tsv");
//		
//		String nonCrypticStd = getResultsTable(sequencesNonCryptic, modelStdAcc);
//		Functions.writeToFile(nonCrypticStd, "results/nonCrypticStd.tsv");
//		
//		String nonCrypticCh = getResultsTable(sequencesNonCryptic, modelChAcc);
//		Functions.writeToFile(nonCrypticCh, "results/nonCrypticCh.tsv");
//		
//		printSomeVariants(variantsNonCrypticAcc);
		
//		// random base change
//		System.out.println("\n---- random changes ---");
//		VariantFile variantsRandomChange = VariantFile.changeRandomBase(variantsPathoTest);
//		Sequences sequencesRandomChg = new Sequences(variantsRandomChange, acceptor);
//		String outRandomChangeChg = getResultsTable(sequencesRandomChg, modelChgAcc);
//		Functions.writeToFile(outRandomChangeChg, "results/randomChangeChg.tsv");
//		String outRandomChangeStd = getResultsTable(sequencesRandomChg, modelStdAcc);
//		Functions.writeToFile(outRandomChangeStd, "results/randomChangeStd.tsv");
		
		
		
		// activating variants
//		System.out.println("\n---- activating variants ---");
//		Sequences modelAcc = new Sequences(variantsPathogene, true);
//		Sequences modelDon = new Sequences(variantsPathogene, false);
//		VariantFile nonCrypticPathogene = VariantFile.extractCrypticVariants(variantsPathogene, modelAcc, modelDon, true);
//		VariantFile activatingVariants = VariantFile.filterActivatingVariants(nonCrypticPathogene, modelAcc, true);
//		System.out.println(activatingVariants);
		
		// cluster activating recursively
//		Sequences strangeSequences = new Sequences(strangeVariants, intronExon);
//		System.out.println(strangeSequences.matrixToString("strange standard matrix 1"));
//		double[] strangeResultsRef = strangeSequences.getIndividualInformation(strangeSequences, true);
//		double[] strangeResultsAlt = strangeSequences.getIndividualInformation(strangeSequences, false);
//		
//		System.out.println("\n\n - results strange - \n");
//		for (int i = 0; i < strangeSequences.size(); i++) {
//			Sequence sequence = strangeSequences.get(i);
//			Variant variant = sequence.getVariant();
//			String ref = variant.getRef();
//			String alt = variant.getAlt();
//			int pos = sequence.getChangePosition() - 30;
//			
//			System.out.println(pos + " " + ref + ">" + alt + "\t " + strangeResultsRef[i] + "\t " + strangeResultsAlt[i]);
//		}
//		
//		for (int j = 2; j < 5; j++) {
//			strangeVariants = strangeVariants.extractStrange(strangeSequences);
//			strangeSequences = new Sequences(strangeVariants, intronExon);
//			System.out.println(strangeSequences.matrixToString("strange standard matrix " + j));
//			strangeResultsRef = strangeSequences.getIndividualInformation(strangeSequences, true);
//			strangeResultsAlt = strangeSequences.getIndividualInformation(strangeSequences, false);
//			
//			System.out.println("\n\n - results strange - \n");
//			for (int i = 0; i < strangeSequences.size(); i++) {
//				Sequence sequence = strangeSequences.get(i);
//				Variant variant = sequence.getVariant();
//				String ref = variant.getRef();
//				String alt = variant.getAlt();
//				int pos = sequence.getChangePosition() - 30;
//				
//				System.out.println(pos + " " + ref + ">" + alt + "\t " + strangeResultsRef[i] + "\t " + strangeResultsAlt[i]);
//			}
//		}
	}
	
	/**
	 * @param variantsP
	 * @param folder
	 */
	private void printAllVariants(Variants variants, String folder) {
		String fileContent = "";
		String line = "pos\t r\t a\t start\t sequence";
		fileContent += line;
		for (int v = 0; v < variants.size(); v++) {
			Variant variant = variants.get(v);
			Sequence sequence = variant.getSequence();
			String ref = variant.getRef();
			String alt = variant.getAlt();
			int posRel = sequence.getPositionChangeRelative();
			int posAbs = variant.getStart();
			line = "\n" + posRel + "\t " + ref + "\t " + alt + "\t " + posAbs + "\t " + sequence.getStringExtended();
			fileContent += line;
		}
		Functions.writeToFile(fileContent, folder + "variants.tsv");
	}

	/**
	 * @param pattern1
	 * @param pattern2
	 * @return
	 */
	private static boolean equalsShifted(String pattern1, String pattern2) {
//		System.out.println("p1:" + pattern1 + "\t p2: " + pattern2);
		if (pattern1.equals(pattern2)) {
			return true;
		}
		if (pattern1.substring(1).equals(pattern2.substring(0, pattern2.length()-1))) {
			return true;
		}
		if (pattern2.substring(1).equals(pattern1.substring(0, pattern1.length()-1))) {
			return true;
		}
		if (pattern1.contains(pattern2)) {
			return true;
		}
		if (pattern2.contains(pattern1)) {
			return true;
		}
//		System.out.println("false");
		return false;
	}

	public static double getRelativePercentage(int quantity, int lengthPattern, int[] numOfPattern) {
		return quantity * ((double)Math.pow(4, lengthPattern) / numOfPattern[lengthPattern/2]);
	}

	/**
	 * @param variantsBenign 
	 * @param variantsPathogene 
	 * @throws UnexpectedException 
	 * 
	 */
	private static void crossValidate(Variants variantsPathogene, Variants variantsBenign, boolean acceptorP, String folder) {
	  int steps = Config.crossValidationSteps + 1;
	  for (int i = 1; i < steps; i++) {
	    CrossValidation crossVal = new CrossValidation(variantsPathogene, variantsBenign, acceptorP, folder, i);
	    Thread thr = new Thread(crossVal);
	    thr.setPriority(Thread.MIN_PRIORITY);
	    thr.start();
	    do {
	      try {
	        Thread.sleep(200);
	      } catch (InterruptedException e) {
	        e.printStackTrace();
	      }
	    } while(Thread.activeCount() >1);// Runtime.getRuntime().availableProcessors());
	  }
	  do {
        try {
          Thread.sleep(1000);
        } catch (InterruptedException e) {
          e.printStackTrace();
        }
      } while(Thread.activeCount() > 1);
	}

	/**
	 * @param sequences
	 * @param pos 
	 * @return
	 */
	private Sequence findSequencebyPosition(Sequences sequences, int pos) {
		Sequence result = null;
		for(int i=0; i<sequences.size() && result == null; i++){
			if (sequences.get(i).getVariant().getStart() == pos) {
				result = sequences.get(i);
			}
		}
		return result;
	}

	/**
	 * @param testSequences
	 * @return
	 */
	private static String getDifferences(Sequences testSequences,  Model modelStd, Model modelCh) {
		double[] resultsReferenceStd = modelStd.getInformation(testSequences, true, false);
		double[] resultsAlternateStd = modelStd.getInformation(testSequences, false, false);
		double[] resultsReferenceCh = modelCh.getInformation(testSequences, true, false);
		double[] resultsAlternateCh = modelCh.getInformation(testSequences, false, false);
		String results = "";
		String line = "pos\t r\t a\t ref\t  alt";
		results += line;
		for (int i = 0; i < testSequences.size(); i++) {
			Sequence sequence = testSequences.get(i);
			Variant variant = sequence.getVariant();
			String ref = variant.getRef();
			String alt = variant.getAlt();
			int pos = sequence.getPositionChangeRelative();
			
			line = "\n" + variant.toString();
			results += line;
		}
		return results;
	}

	private void printVariantsByPosition(Variants variants){
		int sequenceLength = variants.get(0).getSequence().length();
		int[] count = Functions.getInitializedIntArray(sequenceLength);
		
		for (int i = 0; i < variants.size(); i++) {
			count[variants.get(i).getSequence().getPositionJunction()]++;
		}
		
		for (int i = 0; i < count.length; i++) {
			System.out.print(" " + i + " " + count[i] + " ");
			for (int j = 0; j < count[i]; j++) {
				System.out.print(".");
			}
			System.out.println();
		}
	}
	
	public static void printSomeVariants(Variants variants){
		System.out.println("\n");
		System.out.println("0? " + variants.size());
		for (int i = 0, j = 0; i < variants.size(); i++) {
			System.out.println(variants.get(i));
		}
		System.out.println("\n");
	}
}
