/**
 * 
 */
package jsplice.io;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.regex.Matcher;

import jsplice.data.Config;
import jsplice.data.HgmdVariant;
import jsplice.data.RefGene;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.exception.RefGeneNotFoundException;
import jsplice.exception.VariantToFarAwayException;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 * The class handles the ClinVar variant file
 */
public class HgmdFile extends Variants {

	/**
	 * Read the HGMD file specified in GlobalParameters
	 */
	public HgmdFile() {
		readFile(Config.getClinVarFileName());
		initial();
		Log.add("Number of HGMD variants stored:\t" + variants.size(), 3);
	}

	/**
	 * Read the HGMD and store the variants
	 * @param variantFileName Name and path of the file
	 */
	public HgmdFile(String variantFileName) {
		readFile(variantFileName);
		initial();
		Log.add("Number of HGMD variants stored:\t" + variants.size(), 3);
	}

	/**
	 * Read the variant file
	 * 
	 * @param variantFileName HGMD file name
	 *            Destination of the HGMD file with the variants
	 */
	public void readFile(String variantFileName) {
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
			inputStream = new FileInputStream(variantFileName);
			sc = new Scanner(inputStream, "UTF-8");
//			setHashHeader(sc.nextLine());
			while (sc.hasNextLine()) {
				try {
					addVariant(sc.nextLine());
				} catch (IllegalArgumentException e){
					Log.add(e.getMessage(), 2);
				} catch (Exception e) {
					Log.add(e.getMessage(), 5);
				}
			}
			// Scanner suppresses exceptions
			if (sc.ioException() != null) {
				throw sc.ioException();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (inputStream != null) {
				try {
					inputStream.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			if (sc != null) {
				sc.close();
			}
		}
	}

	/**
	 * Check the line and create a HGMD variant
	 * @param line
	 */
	public void addVariant(String line) {
		String[] variantLine = line.split("!");
		if (variantLine.length < 31) {
			throw new IllegalArgumentException("The following variant line doesn't contain all fields:\n\"" + line + "\"");
		}
		// correct DNA change?
		Matcher m = relevantLines.matcher(variantLine[12]);
		
		// correct evidence?
		final ArrayList<String> checkEvidence = new ArrayList<String>(Arrays.asList("DM","DP","DFP"));
		boolean correctEvidence = false;
		if (checkEvidence.contains(variantLine[18])) {
			correctEvidence = true;
		}
		if (m.matches() && correctEvidence) {
			try {
				variants.add(new HgmdVariant(variantLine));
			} catch (java.lang.NumberFormatException e) {
				Log.add(e.getMessage(), 2);
				missingVariantIDs.add(Integer.parseInt(variantLine[24]));
			} catch (IllegalArgumentException e) {
				if (e.getMessage().contains("Illegal DNA change or non-intronic mutation: ")) {
					Log.add(e.getMessage(), 1);
				} else {
					Log.add(e.getMessage(), 2);
				}
			} catch (VariantToFarAwayException e) {
				Log.add(e.getMessage(), 1);
			}catch (Exception e) {
				Log.add("Unexpected exception in line:\n" + line + "\n" + e.getMessage(), 5);
			}
		} else if(!m.matches()){
			Log.add("The variant with HGMD accession " + variantLine[28] + " will be dropped because the HGVS " + variantLine[12] + " doesn't match the constraints of a splice variant.", 1);
		} else if(!correctEvidence){
			Log.add("The variant with HGMD accession " + variantLine[28] + " will be dropped because the evidence " + variantLine[18] + " doesn't match the constraints.", 1);
		}
	}

	/**
	 * Adds a {@link RefGene} Object and the sequence to every {@link Variant}.
	 * @param refGeneFile 
	 * @param faRefFile 
	 */
	public void addRefSeqData(RefSeqGeneFile refGeneFile, FastaReferenceFiles faRefFile){
		int noRefSeqGeneFound = 0;
		for (int i = 0; i < variants.size(); i++) {
			Variant variant = variants.get(i);
			int distance = variant.getDistanceToExon();
			// Will remove variant if it is farer away from junction than the analyzed sequence is long
			if (Math.abs(distance) < Config.getLengthTrainingSequence()) {
				String vChromosome = variant.getChromosome();
				String vGene = variant.getGene();
				int position = variant.getStart();
				try {
					// Find and set the corresponding RefSeq Gene data 
					RefGene refGene = refGeneFile.findByGeneId(vGene, vChromosome, position);
					variant.setRefGene(refGene);
					// get length parameters
					int extendedSequenceLength = Config.getLengthTrainingSequence();
					int extendedExonLength = Config.getLengthTrainingExon();
					int sequenceLength = Config.getLengthModelSequence();
					int exonLength = Config.getLengthModelExon();
					// Add the corresponding sequence to the variant
					String sequence = faRefFile.getSequence(vChromosome, variant.getJunctionPosition(), variant.isAcceptorSite(), refGene.isPlusStrand(), extendedSequenceLength, extendedExonLength);
					variant.setSequence(sequence, extendedExonLength, sequenceLength, exonLength);
					// Check that the sequence contains data
					if (variant.getSequence() == null) {
						Log.add("The sequence is null\nchr: " + vChromosome + "\tjunct: " + variant.getJunctionPosition() + "\tplusStr: " + refGene.isPlusStrand(), 5);
						// remove variants with unequal reference and reference sequence
						String ref = variant.getSequence().getReference();
						if(!variant.getRef().equals(ref)){
							remove(i);
							i--;
							Log.add("Variant with DNA variant " + variant.getDnaHgvs() + " removed because the reference was unequal.", 4);
							noRefSeqGeneFound++;
						}
					}
				} catch (RefGeneNotFoundException e) {
					Log.add("Variant will be removed: " + e.getMessage(), 2);
					remove(i);
					i--;
					noRefSeqGeneFound++;
				}
			}
		}
		Log.add("Number of HGMD variants with no corresponding RefSeq Gene found:\t" + noRefSeqGeneFound, 3);
	}

}
