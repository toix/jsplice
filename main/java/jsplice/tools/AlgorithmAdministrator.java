package jsplice.tools;

import java.rmi.UnexpectedException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import jsplice.data.Config;
import jsplice.data.Sequence;
import jsplice.data.Sequences;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.io.Variants;

/**
 * In AlgorithmAdministrator the different algorithms are applied on the data
*/
public class AlgorithmAdministrator {
	
	Variants variantsBenTrain;
	Variants variantsBenTest;
	Variants variantsPathoTrain;
	Variants variantsPathoTest;
	
	public AlgorithmAdministrator(Variants variantsPathogene, Variants variantsBenign) throws UnexpectedException {
		
		
		Log.add("#benign Variants: " + variantsBenign.size());
		Log.add("#pathogene Variants: " + variantsPathogene.size());
		boolean acceptor = true;
		
		// standard model pathogene
		Log.add("\n - - - pathogene standard model - - - ", 3);
		Model modelStdAcc = new Model(variantsPathogene, acceptor);
		Model modelStdDon = new Model(variantsPathogene, !acceptor);
		
		// Filter
//		Variants variantsPathogeneNonCry = Filter.extractCrypticVariants(variantsPathogene, modelStdAcc, modelStdDon, false);
//		variantsPathogene = VariantFile.concat(Filter.filterActivatingVariants(variantsPathogene, modelStdAcc, true), Filter.filterActivatingVariants(variantsPathogene, modelStdDon, true));
		
		String folder = "results/patternSep/" + Config.getLengthModelIntron() + "+" + Config.getLengthModelExon() + "/";
		crossValidate(variantsPathogene, variantsBenign, acceptor, folder);
		
		
		
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
		
		// analyze certain sequence
//		Log.add("\n---- certain sequence ---", 3);
//		Sequence seqTest = findSequencebyPosition(new Sequences(variantsPathogene, acceptor), 241810764);
//		Log.add(seqTest.getVariant()+"", 3);
//		Log.add(seqTest.getStringExtended(), 2);
//		int crypticPos = 0;//seqTest.getPositionChange()+1;
//		double natRef = modelStdAcc.getIndividualInformation(seqTest, seqTest.getPositionJunction(), true, false).getTotalInformation();
//		double natAlt = modelStdAcc.getIndividualInformation(seqTest, seqTest.getPositionJunction(), false, false).getTotalInformation();
//		double kryRef = modelStdAcc.getIndividualInformation(seqTest, crypticPos, true, false).getTotalInformation();
//		double kryAlt = modelStdAcc.getIndividualInformation(seqTest, crypticPos, false, false).getTotalInformation();
//		Log.add("nat: " + natRef + ">" + natAlt + "\t kry: " + kryRef + ">" + kryAlt + " at " + (crypticPos), 2);
//		Variants variantsTest = new Variants();
//		variantsTest.add(seqTest.getVariant());
//		Filter.extractCrypticVariants(variantsTest, modelStdAcc, modelStdDon, false);
		
		
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
	private void crossValidate(Variants variantsPathogene, Variants variantsBenign, boolean acceptor, String folder) throws UnexpectedException {
		for (int i = 1; i < 301; i++) {
			System.out.println("\n\n - - - - - - - - - - - - - - - - - - - - - \n");
			// separate training and test variants
			Log.add("cross validation run nr. " + i, 3);
			ArrayList<Variants> separatedVariants;
			separatedVariants = separateData(variantsBenign);
			variantsBenTrain = separatedVariants.get(0);
			variantsBenTest = separatedVariants.get(1);
			separatedVariants = separateData(variantsPathogene);
			variantsPathoTrain = separatedVariants.get(0);
			variantsPathoTest = separatedVariants.get(1);
			
//			System.out.println("\n\n - - - - - - - - - - - - - - - - - - - - - \n");
			// VariantFile trainVariantsBenign = clinVarBenign.filteVariants(intronExon);
			// trainVariantsBenign = trainVariantsBenign.filterNonACGT();
			// Sequences standardModel = new Sequences(trainVariantsBenign, intronExon, false);
			// standard model pathogene
			Log.add("\n - - - pathogene standard model - - - ", 3);
			Model modelStdAcc = new Model(variantsPathoTrain, acceptor);
			Log.add("Number of pathogene training sequences for acceptor site: " + modelStdAcc.getSequences().size(), 3);
			Functions.writeToFile(modelStdAcc.matrixToString("acceptor matrix", modelStdAcc.getJunctionPosition()),
					folder + "accMatrix.tsv");
			Model modelStdDon = new Model(variantsPathoTrain, !acceptor);
			Log.add("Number of pathogene training sequences for donor site: " + modelStdDon.getSequences().size(), 3);
			Functions.writeToFile(modelStdAcc.matrixToString("donor matrix", modelStdDon.getJunctionPosition()),
					folder + "donMatrix.tsv");
			// create pathogene test sequences
			Log.add("\n - - - pathogene test variants - - - ", 3);
			// variantsPathoTest = VariantFile.extractCrypticVariants(variantsPathoTest, modelStdAcc, modelStdDon, false);
			Sequences sequencesPathoTestAcc = new Sequences(variantsPathoTest, acceptor);
			Log.add("Number of pathogene test sequences for acceptor site: " + sequencesPathoTestAcc.size(), 3);
			// create benign test sequences
			Log.add("\n - - - benign test variants - - - ", 3);
			Sequences sequencesBenTestAcc = new Sequences(variantsBenTest, acceptor);
			
			// write benign results to file
			String standardResultsBenign = getResultsTable(sequencesBenTestAcc, modelStdAcc, false);
			Functions.writeToFile(standardResultsBenign, folder + "benAccStd" + i + ".tsv");
			String changeResultsBeingn = getResultsTable(sequencesBenTestAcc, modelStdAcc, true);
			Functions.writeToFile(changeResultsBeingn, folder + "benAccChg" + i + ".tsv");
			// write pathogene results to file
			String standardResultsPathogene = getResultsTable(sequencesPathoTestAcc, modelStdAcc, false);
			Functions.writeToFile(standardResultsPathogene, folder + "pathoAccStd" + i + ".tsv");
			String changeResultsPathogene = getResultsTable(sequencesPathoTestAcc, modelStdAcc, true);
			Functions.writeToFile(changeResultsPathogene, folder + "pathoAccChg" + i + ".tsv");
			Log.writeToFile();
		}
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
	 * @return two Variant Files; The first is the Train Data; The second is the test Data
	 */
	private static ArrayList<Variants> separateData(Variants variantsP) {
		ArrayList<Variants> separatedVariants = new ArrayList<Variants>();
		// create train and test variants
		separatedVariants.add(new Variants(variantsP));		// separatedVariants.get(0) = Train
		separatedVariants.add(new Variants());					// separatedVariants.get(1) = Test
		int testSize = separatedVariants.get(0).size()/15;
		// create array with unique random indices
		ArrayList<Integer> orderedList = new ArrayList<Integer>();
		for (int i = 0; i < separatedVariants.get(0).size(); i++) {
			orderedList.add(new Integer(i));
		}
        Collections.shuffle(orderedList);
        List<Integer> testIndices = orderedList.subList(0, testSize);
        Collections.sort(testIndices);
        // separate test and train variants
        for (int i = testSize-1; i >= 0 ; i--) {
        	separatedVariants.get(1).add(separatedVariants.get(0).get(testIndices.get(i)));
		}
        for (int i = testSize-1; i >= 0 ; i--) {
        	separatedVariants.get(0).remove(testIndices.get(i));
        }
		return separatedVariants;
	}

	/**
	 * @param testSequences
	 * @param cluster 
	 * @return
	 */
	private static String getResultsTable(Sequences testSequences,  Model modelStdAcc, boolean cluster) {
		double[] resultsReference = modelStdAcc.getIndividualInformation(testSequences, true, cluster);
		double[] resultsAlternate = modelStdAcc.getIndividualInformation(testSequences, false, cluster);
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

	/**
	 * @param testSequences
	 * @return
	 */
	private static String getDifferences(Sequences testSequences,  Model modelStd, Model modelCh) {
		double[] resultsReferenceStd = modelStd.getIndividualInformation(testSequences, true, false);
		double[] resultsAlternateStd = modelStd.getIndividualInformation(testSequences, false, false);
		double[] resultsReferenceCh = modelCh.getIndividualInformation(testSequences, true, false);
		double[] resultsAlternateCh = modelCh.getIndividualInformation(testSequences, false, false);
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
