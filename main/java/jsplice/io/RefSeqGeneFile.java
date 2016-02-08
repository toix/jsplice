package jsplice.io;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Pattern;

import jsplice.data.RefGene;
import jsplice.exception.Log;
import jsplice.exception.RefGeneNotFoundException;

import com.opencsv.CSVReader;
/**
 * Reads and administrates the metadata about RefSeq from RefGene file. <br/>
 * First line contains identifiers
 * @author Tobias Gresser (gresserT@gmail.com)
 */
//	uint bin;            	"Indexing field to speed chromosome range queries"
//	string name;        	"Name of gene (usually transcript_id from GTF)"
//	string chrom;       	"Chromosome name"
//	char[1] strand;     	"+ or - for strand"
//	uint txStart;       	"Transcription start position"
//	uint txEnd;         	"Transcription end position"
//	uint cdsStart;      	"Coding region start"
//	uint cdsEnd;        	"Coding region end"
//	uint exonCount;     	"Number of exons"
//	uint[exonCount] exonStarts; "Exon start positions"
//	uint[exonCount] exonEnds;   "Exon end positions"
//	uint score;           "score"
//	string name2;       	"Alternate name (e.g. gene_id from GTF)"
//	string cdsStartStat; 	"enum('none','unk','incmpl','cmpl')"
//	string cdsEndStat;   	"enum('none','unk','incmpl','cmpl')"
//	lstring exonFrames; 	"Exon frame {0,1,2}, or -1 if no frame for exon"

public class RefSeqGeneFile {
	
	/**
	 * Names of the RefGene file titles
	*/
	private final String[] identifiers = {"bin", "name", "chrom", "strand",
			"txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount",
			"exonStarts", "exonEnds", "score", "name2", "cdsStartStat",
			"cdsEndStat", "exonFrames" };
	/**
	 * HashMap from title names to the title index of the RefGene file
	*/
	private HashMap<String, Integer> titles = new HashMap<>();
	/**
	 * Tab separated lines of the file as array
	 */
	private ArrayList<String[]> refGeneLines;
	/**
	 * Hashes to a {@link RefGene} by RefSeq accession number
	 */
	private HashMap<String, RefGene> hashByRefSeq = new HashMap<String, RefGene>();
	/**
	 * Hashes to a {@link RefGene} by gene id
	 */
	private ArrayList<RefGene> duplicateRefSeqs = new ArrayList<RefGene>();
	/**
	 * {@link ArrayList} of {@link RefGene} which share the same RefSeq accession with an element in RefGeneMap
	 */
	private HashMap<String, RefGene> hashByGeneId = new HashMap<String, RefGene>();
	/**
	 * {@link ArrayList} of {@link RefGene} which share the same RefSeq accession with an element in RefGeneMap
	 */
	private ArrayList<RefGene> duplicateGeneIds = new ArrayList<RefGene>();
	/**
	 * Tab separated lines of the file as array which don't match the constraints
	 */
	private ArrayList<RefGene> invalidRefGeneLines = new ArrayList<RefGene>();
	
	/**
	 * Creates the HashMap titles with the titles. </br>
	 * Reads the whole RefGene file </br>
	 * Creates the HashMap lines with the RefSeq accession number
	 * @param refSeqGeneFileName 
	 * @throws IOException
	 */
	public RefSeqGeneFile(String refSeqGeneFileName) throws IOException {
		// Create HashMap titles
		for (int i = 0; i < identifiers.length; i++) {
			titles.put(identifiers[i], i);
		}
		// Read the file
		CSVReader reader = new CSVReader(new FileReader(refSeqGeneFileName), '\t');
		refGeneLines = (ArrayList<String[]>) reader.readAll();
		reader.close();
		deleteInvalidData();
		// iterate over all RefSeq Gene lines and create Hash
		for (int i = 1; i < refGeneLines.size(); i++) {
			try {
					
				RefGene refGene = new RefGene(refGeneLines.get(i));
				// add to Hash by RefSeq accession
				if(!hashByRefSeq.containsKey(refGene.getRefSeqAccession())){
					hashByRefSeq.put(refGene.getRefSeqAccession(), refGene);
				} else {
					duplicateRefSeqs.add(refGene);
				}
				// add to Hash by Gene ID
				if(!hashByGeneId.containsKey(refGene.getGene())){
					hashByGeneId.put(refGene.getGene(), refGene);
				} else {
					duplicateGeneIds.add(refGene);
				}
			} catch (IllegalArgumentException e) {
				// Try to add the line to invalidRefGeneLines accepting uncommon chromosome names
				try {
					invalidRefGeneLines.add(new RefGene(refGeneLines.get(i), false));
				} catch (IllegalArgumentException e2) {
					Log.add(e2.getStackTrace().toString(), 2);
				}
			}
		}
		Log.add("Number of RefSeq Gene stored in hash:\t" + size(), 3);
		Log.add("Number of duplicate RefSeq accession found:\t" + duplicateRefSeqs.size(), 3);
		Log.add("Number of duplicate gene name found:\t" + duplicateGeneIds.size(), 3);
		Log.add("Number of invalid RefSeq Gene found:\t" + invalidRefGeneLines.size(), 3);
	}

	/**
	 * Validate data and move invalid lines in InvalidRefGeneLines
	 */
	private void deleteInvalidData() {
		for (int i = refGeneLines.size()-1; i >= 0 ; i--) {
			// validate chromosome
			String chromosome =refGeneLines.get(i)[titles.get("chrom")];
			String accession =refGeneLines.get(i)[titles.get("name")];
			boolean invalidAccession = accession.subSequence(0, 3).equals("NR_");
			boolean invalidChromosome = !Pattern.matches("^(chr)?([1-9]|1[0-9]|2[0-2]|[XY])?.*$", chromosome);
			if(invalidAccession || invalidChromosome){
				refGeneLines.remove(i);
			}
		}
	}
	
//	/**
//	 * Get a the Entry of the line containing the accession number
//	 * 
//	 * @param accession RefSeq accession number
//	 * @param entry name of the entry in the line
//	*/
//	public String get(String entry, String accession) {
//		return RefGeneLines.get(lines.get(accession))[titles.get(entry)];
//	}
//	
//	/**
//	 * Get a the Entry of the line containing the accession number
//	 * 
//	 * @param accession RefSeq accession number
//	 * @param entry number of the entry in the line
//	 */
//	public String get(int index, String accession){
//		return RefGeneLines.get(lines.get(accession))[index];
//	}
	
//	/**
//	 * Returns a RefSeqInformation Object which contains the accession number
//	 * 
//	 * @param accession
//	 *            RefSeq accession number
//	 * @return The first {@link RefSeq} matching the accession Number
//	 * @throws RefSeqNotFoundException
//	 */
//	public RefSeq get(String accession) throws RefSeqNotFoundException {
//		if (hashByRefSeq.containsKey(accession)) {
//			return hashByRefSeq.get(accession);
//		} else {
//			throw new RefSeqNotFoundException(
//					"No RefSeqInformation found for RefSeq accession "
//							+ accession);
//		}
//	}
	
	/**
	 * Returns a RefSeq Gene Object which contains the accession number
	 * and where the position lies in the transcribed region of the RefSeq <br/>
	 * Also searches for other RefSeqInformation containing the same RefSeq
	 * accession number
	 * 
	 * @param geneId RefSeq accession number
	 * @param geneId Gene ID
	 * @param position Absolute variant position on the chromosome
	 * @param chromosome Chromosome of the RefSeq
	 * @return The first {@link RefGene} matching the accession Number
	 *         and containing the Position
	 * @throws RefGeneNotFoundException
	 */
	public RefGene find(String accession, String geneId, String chromosome, int position) throws RefGeneNotFoundException{
		RefGene GeneByName = null;
		RefGene GeneByRefSeq = null;
		try {
			GeneByRefSeq = findByRefSeq(accession, chromosome, position);
		} catch (Exception e) {
			Log.add("No RefSeq gene fond by RefSeq: " + accession, 1);
		}
		
		try {
			GeneByName = findByGeneId(geneId, chromosome, position);
		} catch (Exception e) {
			Log.add("No RefSeq gene fond by gene id: " + geneId, 1);
		}
		if(GeneByRefSeq != null || GeneByName != null){
			if(GeneByRefSeq != null){
				if(GeneByName != null && GeneByName != GeneByRefSeq){
					Log.add("A different RefSeq Gene was found by gene name\t RefSeq accession: " + accession + "\t gene ID: " + geneId + "\t position: " + position, 1);
				}
				return GeneByRefSeq;
			} else {
				return GeneByName;
			}
		} else {
			throw new RefGeneNotFoundException("No RefSeq Gene entery found for RefSeq accession " + accession + ", gene ID " + geneId + " and position " + position);
		}
	}
	
	/**
	 * Returns a RefSeqInformation Object which contains the accession number
	 * and where the position lies in the transcribed region of the RefSeq <br/>
	 * Also searches for other RefSeqInformation containing the same RefSeq
	 * accession number
	 * 
	 * @param accession RefSeq accession number
	 * @param position Absolute variant position on the chromosome
	 * @param chromosome Chromosome of the RefSeq
	 * @return The first {@link RefGene} matching the accession Number
	 *         and containing the Position
	 * @throws RefGeneNotFoundException
	 */
	public RefGene findByRefSeq(String accession, String chromosome, int position) throws RefGeneNotFoundException{
		if(hashByRefSeq.containsKey(accession)){
			RefGene refSeqInfo = hashByRefSeq.get(accession);
			if(refSeqInfo.hasJunctionNextTo(position) && refSeqInfo.getChromosome().equals(chromosome)){
				return refSeqInfo;
			} else {
				for (int i = 0; i < duplicateRefSeqs.size(); i++) {
					if (	duplicateRefSeqs.get(i).getRefSeqAccession().equals(accession)
							&& duplicateRefSeqs.get(i).hasJunctionNextTo(position)
							&& duplicateRefSeqs.get(i).equals(chromosome)) {
						return duplicateRefSeqs.get(i);
					}
				}
			}
		}
		throw new RefGeneNotFoundException("No RefSeq Gene entery found for RefSeq accession " + accession + " and position " + position);
	}

	/**
	 * Returns a RefSeqInformation Object which contains the accession number
	 * and where the position lies in the transcribed region of the RefSeq <br/>
	 * Also searches for other RefSeqInformation containing the same RefSeq
	 * accession number
	 * 
	 * @param geneId Gene ID
	 * @param position Absolute variant position on the chromosome
	 * @param chromosome Chromosome of the RefSeq
	 * @return The first {@link RefGene} matching the accession Number
	 *         and containing the Position
	 * @throws RefGeneNotFoundException
	 */
	public RefGene findByGeneId(String geneId, String chromosome, int position) throws RefGeneNotFoundException{
		if(hashByGeneId.containsKey(geneId)){
			RefGene refSeqInfo = hashByGeneId.get(geneId);
			if(refSeqInfo.hasJunctionNextTo(position) && refSeqInfo.getChromosome().equals(chromosome)){
				return refSeqInfo;
			} else {
				for (int i = 0; i < duplicateGeneIds.size(); i++) {
					if (	duplicateGeneIds.get(i).getRefSeqAccession().equals(geneId)
							&& duplicateGeneIds.get(i).hasJunctionNextTo(position)
							&& duplicateGeneIds.get(i).equals(chromosome)) {
						return duplicateGeneIds.get(i);
					}
				}
			}
		}
		throw new RefGeneNotFoundException("No RefSeq Gene entery found for gene ID " + geneId + " and position " + position);
	}

	/**
	 * Returns a list of RefSeqInformation Objects which contain the accession
	 * number and where the position lies in the transcribed region of the
	 * chromosome RefSeq
	 * 
	 * @param accession
	 *            RefSeq accession number
	 * @param chromosome
	 * @param position
	 * @return
	 * @throws Exception
	 */
	public ArrayList<RefGene> getInvalidRefGeneLines(String accession, String chromosome, int position) throws Exception {
		ArrayList<RefGene> invLines = new ArrayList<RefGene>();
		for (int i = 0; i < invalidRefGeneLines.size(); i++) {
			String invAccession = invalidRefGeneLines.get(i).getRefSeqAccession();
			String invChromosome = invalidRefGeneLines.get(i).getChromosome();
			int invTranscriptionStart = invalidRefGeneLines.get(i).getTranscriptionStart();
			int invTranscriptionEnd = invalidRefGeneLines.get(i).getTranscriptionEnd();
			if(invTranscriptionStart<position && invTranscriptionEnd>position && invAccession.equals(accession) && invChromosome.contains(chromosome)){
				invLines.add(invalidRefGeneLines.get(i));
			}
		}
		if (invLines.size() == 0)
			throw new RefGeneNotFoundException("No invalid line matches the paramters.");
		return invLines;
	}
	
	public int size(){
		return hashByRefSeq.size();
	}

//	public String getChromosome(String accession){
//		return get("chrom", accession);
//	}
//	public boolean isPlusStrand(String accession) throws Exception{
//		String strand = get("strand", accession);
//		if (strand.equals("+")){
//			return true;
//		} else if (strand.equals("-")) {
//			return false;
//		} else {
//			throw new Exception("Invalid Strand orientation:\t" + strand);
//		}
//	}
//	public int getTranscriptionStart(String accession){
//		return Integer.parseInt(get("txStart", accession));
//	}
//	public int getTranscriptionEnd(String accession){
//		return Integer.parseInt(get("txEnd", accession));
//	}
//	public int getExonCount(String accession){
//		return Integer.parseInt(get("exonCount", accession));
//	}
//	public ArrayList<Integer> getExonStarts(String accession){
//		String[] exonStartStrings = get("exonStarts", accession).split(",");
//		ArrayList<Integer> exonStarts = new ArrayList<Integer>();
//		for (int i = 0; i < exonStarts.size(); i++) {
//			exonStarts.add(Integer.parseInt(exonStartStrings[i]));
//		}
//		return exonStarts;
//	}
//	public ArrayList<Integer> getExonEnds(String accession){
//		String[] exonEndStrings = get("exonStarts", accession).split(",");
//		ArrayList<Integer> exonEnds = new ArrayList<Integer>();
//		for (int i = 0; i < exonEnds.size(); i++) {
//			exonEnds.add(Integer.parseInt(exonEndStrings[i]));
//		}
//		return exonEnds;
//	}
}
