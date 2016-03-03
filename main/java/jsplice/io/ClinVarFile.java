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

import jsplice.data.ClinVariant;
import jsplice.data.Config;
import jsplice.data.RefGene;
import jsplice.data.Variant;
import jsplice.exception.Log;
import jsplice.exception.MissingDataException;
import jsplice.exception.RefGeneNotFoundException;
import jsplice.exception.VariantToFarAwayException;



/**
 * @author Tobias Gresser (gresserT@gmail.com)
 * The class handles the ClinVar variant file
 */
public class ClinVarFile extends Variants {
	
	/**
	 * Read in the malign variants of ClinVar file specified in GlobalParameters and store the
	 * information
	 */
	public ClinVarFile() {
		readFile(Config.getClinVarFileName(), false);
		initial();
		Log.add("Number of ClinVar variants stored:\t" + variants.size(), 3);
	}
	/**
	 * Read the ClinVar file and store the information
	 * @param variantFileName Name and path of the file
	 * @param benign whether to read benign or pathogene variants
	 */
	public ClinVarFile(String variantFileName, boolean benign) {
		readFile(variantFileName, benign);
		initial();
		Log.add("Number of ClinVar variants stored:\t" + variants.size(), 3);
	}

	/**
	 * Read the variant file
	 * 
	 * @param variantFileName ClinVar file name
	 *            Destination of the ClinVar file with the variants
	 * @param clinVar True if method reads the ClinVar file. False if method reads the HGVS file.
	 */
	public void readFile(String variantFileName, boolean benign) {
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
			inputStream = new FileInputStream(variantFileName);
			sc = new Scanner(inputStream, "UTF-8");
//			setHashHeader(sc.nextLine());
			while (sc.hasNextLine()) {
				try {
					addVariant(sc.nextLine(), benign);
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
	 * Check the line and create a ClinVar variant
	 * @param line
	 */
	public void addVariant(String line, boolean benign) {
		String[] variantLine = line.split("\t");
		if (variantLine.length < 25) {
			throw new IllegalArgumentException("The following variant line doesn't contain all fields:\n\"" + line + "\"");
		}
		// matches (splice) HGVS constraints?
		Matcher nameMatcher = relevantLines.matcher(variantLine[2]);
		Matcher hgvsMatcher = relevantLines.matcher(variantLine[18]);
		boolean hgvsMatches = nameMatcher.matches() || hgvsMatcher.matches();
		// Variant has chosen assembly version?
		boolean validAssembly = Config.useGrch38() == (variantLine[12].equals("GRCh38"));
		// correct clinical significance?
		final ArrayList<String> checkClinicalSignificance;
		if (benign) {
			checkClinicalSignificance = new ArrayList<String>(Arrays.asList("Benign", "Likely benign"));
		} else {
			checkClinicalSignificance = new ArrayList<String>(Arrays.asList("Pathogenic", "Likely pathogenic"));
		}
		ArrayList<String> clinicalSignificance = new ArrayList<String>(Arrays.asList(variantLine[5].split(";")));
		boolean correctSignificance = false;
		for (int i = 0; i < clinicalSignificance.size(); i++) {
			correctSignificance = correctSignificance || checkClinicalSignificance.contains(clinicalSignificance.get(i));
		}
		if (hgvsMatches && validAssembly && correctSignificance) {
			try {
				variants.add(new ClinVariant(variantLine, Config.useGrch38(), benign));
			} catch (java.lang.NumberFormatException e) {
				Log.add(e.getMessage(), 2);
				missingVariantIDs.add(Integer.parseInt(variantLine[24]));
			} catch (IllegalArgumentException e) {
				if (e.getMessage().contains("or non-intronic mutation: ")) {
					Log.add(e.getMessage(), 1);
				} else {
					Log.add(e.getMessage(), 2);
				}
			} catch (MissingDataException e) {
				Log.add(e.getMessage(), 1);
			} catch (VariantToFarAwayException e) {
				Log.add(e.getMessage(), 1);
			} catch (Exception e) {
				Log.add("Unexpected exception in line:\n" + line + "\n" + e.getMessage(), 5);
			}
		} else {
			if (!validAssembly) {
				Log.add("The variant with allele id " + variantLine[0] + " will be dropped because the assembly " + variantLine[12] + "  doesn't reference the configured assembly.", 1);
			} else if (!hgvsMatches) {
				Log.add("The variant with allele id " + variantLine[0] + " will be dropped because the HGVS " + variantLine[18] + " doesn't match the constraints of a splice variant.", 1);
			} else if (!correctSignificance) {
				Log.add("The variant with allele id " + variantLine[0] + " will be dropped because the clinical Significance " + variantLine[5] + " doesn't match the constraints of a splice variant.", 1);
			} else {
				Log.add("This should be unreachable. Please check your code!", 4);
			}
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
				String vAccession = variant.getRefSeqAccession();
				String vChromosome = variant.getChromosome();
				String vGene = variant.getGene();
				int position = variant.getStart();
				try {
					// Find and set the corresponding RefSeq Gene data 
					RefGene refGene = refGeneFile.find(vAccession, vGene, vChromosome, position);
					variant.setRefGene(refGene);
					// get length parameters
					int extendedSequenceLength = Config.getLengthTrainingSequence();
					int extendedExonLength = Config.getLengthTrainingExon();
					int sequenceLength = Config.getLengthModelSequence();
					int exonLength = Config.getLengthModelExon();
					// Add the corresponding sequence to the variant
					String sequence = faRefFile.getSequence(vChromosome, variant.getJunctionPosition(), variant.isAcceptorSite(), refGene.isPlusStrand(), extendedSequenceLength, extendedExonLength);
					variant.setSequence(sequence, extendedExonLength, sequenceLength, exonLength);
					// check that the sequence contains data
					if (variant.getSequence() == null) {
						Log.add("The sequence is null\nchr: " + vChromosome + "\tjunct: " + variant.getJunctionPosition() + "\tplusStr: " + refGene.isPlusStrand(), 4);
						// remove variants with unequal reference and reference sequence
						String ref = variant.getSequence().getReference();
						if(!variant.getRef().equals(ref)){
							remove(i);
							i--;
							Log.add("Variant with DNA variant " + variant.getDnaHgvs() + " removed because the reference was unequal.", 4);
							noRefSeqGeneFound++;
						}
					}
				// remove variants that have no RefSeq Gene entry
				} catch (RefGeneNotFoundException e) {
					Log.add("Variant will be removed: " + e.getMessage(), 2);
					remove(i);
					i--;
					noRefSeqGeneFound++;
				}
			// remove variants that are to far away from exon junction
			} else {
				remove(i);
				i--;
				
			}
		}
		Log.add("Number of ClinVariants with no corresponding RefSeq Gene found:\t" + noRefSeqGeneFound, 3);
	}

}
