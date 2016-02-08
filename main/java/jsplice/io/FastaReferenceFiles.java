package jsplice.io;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;

import jsplice.data.Config;
import jsplice.tools.Functions;

/**
 * Administrates the 24 chromosome files
 * 
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class FastaReferenceFiles {

	HashMap<String, IndexedFastaSequenceFile> chromosomeFiles = new HashMap<String, IndexedFastaSequenceFile>();

	public FastaReferenceFiles(String chrFileNames) {
		String[] chromosomeNames = Config.getChromosomenames();
		for (int i = 0; i < chromosomeNames.length; i++) {
			String fileName = chrFileNames.replace("*", chromosomeNames[i]);
			File faFile = new File(fileName);
			File faiFile = new File(fileName + ".fai");
			chromosomeFiles.put(chromosomeNames[i],
					new IndexedFastaSequenceFile(faFile,
							new FastaSequenceIndex(faiFile)));
		}
	}

	/**
	 * Return the sequence of the reference file defined by the parameters 
	 * @param chromosome
	 * @param start inclusive
	 * @param stop inclusive
	 * @param isPlusStrand
	 * @return
	 */
	public String getSequence(String chromosome, int start, int stop, boolean isPlusStrand) {
		// Delete the first three character if chromosome parameter starts with "chr"
		if (chromosome.length() > 2 && chromosome.substring(0, 3).equals("chr")) {
			chromosome = chromosome.substring(3, chromosome.length());
		}
		// Throw an Exception if chromosome parameter is not a valid chromosome name
		if (!Arrays.asList(Config.getChromosomenames()).contains(chromosome)) {
			throw new IllegalArgumentException("\"" + chromosome
					+ "\" or \"chr" + chromosome
					+ "\" is not a valid chromosome name.");
		}
		ReferenceSequence refSeq = chromosomeFiles.get(chromosome)
				.getSubsequenceAt("chr" + chromosome, start, stop);
		String sequence = new String(refSeq.getBases());
		if(!isPlusStrand)
			sequence = Functions.reverseAndComplement(sequence);
		return sequence;
	}

	/**
	 * Return the sequence of the reference file surrounding the junction position
	 * @param chromosome
	 * @param junctionPosition
	 * @param intronExonJunction
	 * @param isPlusStrand
	 * @param length
	 * @param exonLength
	 * @return
	 */
	public String getSequence(String chromosome, int junctionPosition, boolean intronExonJunction, boolean isPlusStrand, int sequenceLength, int exonLength) {
		int start;
		int stop;
		if (isPlusStrand == intronExonJunction) {
			start = junctionPosition - sequenceLength + exonLength;
			stop = junctionPosition - 1 + exonLength;
		} else {
			start = junctionPosition + 1 - exonLength;
			stop = junctionPosition + sequenceLength - exonLength;
		}
		return getSequence(chromosome, start, stop, isPlusStrand);
	}

	/*
	 * public void close(){ for (int i = 0; i < chromosomeFiles.size(); i++) {
	 * try { chromosomeFiles.get(i).close(); } catch (IOException e) {
	 * e.printStackTrace(); } catch (NullPointerException e) {
	 * System.out.println(chromosomeNames[i]); e.printStackTrace(); } } }
	 */
}
