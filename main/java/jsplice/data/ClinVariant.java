/**
 * 
 */
package jsplice.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Pattern;

import jsplice.exception.DebuggingException;
import jsplice.exception.Log;
import jsplice.exception.MissingDataException;
import jsplice.exception.VariantToFarAwayException;
import jsplice.tools.Functions;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 */
public class ClinVariant extends Variant {

	/**
	 * ID of the Variant; alleleID field in ClinVa
	 */
	private int alleleId;
	/**
	 * type of variation
	 */
	private String type;
	/**
	 * preferred name for the variation 
	 */
	private String name;
	/**
	 * Gene ID in NCBI's Gene database
	 */
	private int geneID;
	/**
	 * list of values of clinical significance reported for this variation
	 */
	private ArrayList<String> clinicalSignificance = new ArrayList<String>();
	/**
	 * index of region in dbVar; stored as int without nsv-prefix, -1 if not available
	 */
	private int nsv;
	/**
	 * accession number of the submission reference (version) in ClinVar;
	 * stored as int without "RCV"-prefix
	 */
	private ArrayList<Integer> rcv = new ArrayList<Integer>();
	/**
	 * if there is a test registered as specific to this variation in the NIH 
	 * genetic Testing Registry
	 */
	private boolean testedInGTR;
	/**
	 * list of database names and identifiers for phenotype(s) reported for this
	 * variant
	 */
	private ArrayList<String> phenotypeIDs;
	/**
	 * list of all allelic origins for this variation
	 */
	private ArrayList<String> origin;
	/**
	 * highest review status for reporting this measure
	 */
	private ArrayList<String> reviewStatus = new ArrayList<String>();
	/**
	 * protein-based HGVS expression with RefSeq as reference. 
	 * Ignored: There is no protein-based HGVS for splice mutation. 
	 * private String pHGVS;
	 */
	// number of submissions with this variant
	private int numberSubmitters;
	/**
	 * the latest time any submitter reported clinical significance
	 */
	private String lastEvaluated;
	/**
	 * list of other identifiers or sources of information about this variant
	 */
	private ArrayList<String> otherIDs;
	/**
	 * the value used to build the URL for the current default report
	 */
	private int variantID;
	/**
	 * Pattern for Substitutions
	 * NM_ - Prefix of RefSeq mRNA accession number
	 */
	final Pattern checkSubstitution = Pattern
			.compile("^NM_[0-9]{6,9}\\.[1-9][0-9]?:c\\."
					+ "[1-9][0-9]*[+-][1-9][0-9]*"
					+ "([ACGT]>[ACGT])$");
	/**
	 * Pattern for Deletions, Duplications and Insertions. </br>
	 * NM_ is the Prefix of RefSeq mRNA accession number </br>
	 * (no Inversions and no Translocations) </br>
	 */
	final Pattern checkOtherHGVS = Pattern
	.compile("^NM_[0-9]{6,9}\\.[1-9][0-9]?:c\\."
			+ "([1-9][0-9]*[+-][1-9][0-9]*(_[1-9][0-9]*([+-][1-9][0-9]*)?)?"		// 0+0(_0(+0))
			+ "|[1-9][0-9]*_[1-9][0-9]*[+-][1-9][0-9]*)"							// 0_0+0
			+ "(del|dup|(del[ACGT]*)?ins[ACGT])[ACGT]*$");
	/**
	 * HGVS expression with LRG as reference. Only if used for cHGVS
	 */
	String cLRG;
	/**
	 * @param lines
	 * @param annotationTitles
	 * @param irrelevantAnnotations
	 * @param useGrch38
	 * @throws NumberFormatException
	 * @throws IllegalArgumentException
	 */
	public ClinVariant(String[] lines, boolean useGrch38, boolean benign) throws NumberFormatException, IllegalArgumentException, DebuggingException, VariantToFarAwayException {
		super(lines, useGrch38);
		// if a relevant annotation is empty or contains no data, throw IllegalArgumentException
		List<Integer> relevantAnnotations = Arrays.asList(0, 2, 4, 5, 12, 13, 14, 15, 18, 24);
		for (int i = 0; i < lines.length; i++) {
			if (relevantAnnotations.contains(i) && Pattern.matches("^[\\s\\-]*$", lines[i])){
				throw new MissingDataException("The following row contains no data: " + i + "\t AlleleID: " + lines[0]);
			}
		}
		// if cHGVS contains LRG reference, it will be replaced by name field
		if(lines[18].length()>2 && lines[18].substring(0, 3).equals("LRG")){
			cLRG =lines[18];
			lines[18] = lines[2].replaceAll("\\(.*\\)", "");
		}
		// delete the number in the end for deletions and duplications
		if(Pattern.matches(".*(del|dup)[0-9]*$", lines[18]))
			lines[18].replaceAll("[0-9]*$", "");
		// check for single nucleotide splice variant or throw exception
		if(!checkSubstitution.matcher(lines[18]).matches())// && !checkOtherHGVS.matcher(annotationLine[18]).matches())
			throw new IllegalArgumentException("Illegal HGVS or non-intronic mutation: " + lines[18] + "\t AlleleID: " + lines[0]);
		// allele ID
		alleleId = Integer.parseInt(lines[0]);
		// type
		final ArrayList<String> checkType = new ArrayList<String>(Arrays.asList("copy number gain", "copy number loss",
				"indel", "duplication", "fusion", "deletion", "insertion", "inversion", "NT expansion", "protein only",
				"short repeat", "single nucleotide variant", "structural variant", "undetermined variant"));
		if (checkType.contains(lines[1]))
			type = lines[1];
		else
			throw new IllegalArgumentException("Illegal Type: " + lines[1] + "\t AlleleID: " + lines[0]);
		// name
		name = lines[2];
		// gene ID
		geneID = Integer.parseInt(lines[3]);
		// gene
		gene = lines[4].split(";")[0];
		// clinical significance
//		final ArrayList<String> checkClinicalSignificance = new ArrayList<String>(Arrays.asList("association",
//				"Benign", "Likely benign", "Pathogenic", "Likely pathogenic", "risk factor", "drug response",
//				"protective", "Uncertain significance", "other", "confers sensitivity", "not provided"));
		final ArrayList<String> checkClinicalSignificance;
		if (benign) {
			checkClinicalSignificance = new ArrayList<String>(Arrays.asList("Benign", "Likely benign"));
		} else {
			checkClinicalSignificance = new ArrayList<String>(Arrays.asList("Pathogenic", "Likely pathogenic", "risk factor"));
		}
		String[] clinicalSignificanceTemp = lines[5].split(";");
		for (int i = 0; i < clinicalSignificanceTemp.length; i++) {
			if (checkClinicalSignificance.contains(clinicalSignificanceTemp[i])) {
				clinicalSignificance.add(clinicalSignificanceTemp[i]);
			}
		}
		if (clinicalSignificance.size() < 1) {
			throw new IllegalArgumentException("Illegal clinical significance: " + lines[5] + "\t AlleleID: " + lines[0]);
		}
		// snp
		snp = lines[6];
		// nsv
		if(lines[7].equals("-"))
			nsv = -1;
		else
			nsv = Integer.parseInt(lines[7].substring(3));
		// RCV
		String[] rcvTemp = lines[8].split(";");
		for (int i = 0; i < rcvTemp.length; i++) {
			rcv.add(Integer.parseInt(rcvTemp[i].substring(3)));
		}
		// Tested in GTR
		if (Pattern.matches("[Yy](es)?", lines[9]))
			testedInGTR = true;
		else if (Pattern.matches("[Nn]o?", lines[9]))
			testedInGTR = false;
		else
			throw new IllegalArgumentException("Illegal GTR test state: " + lines[9] + "\t AlleleID: " + lines[0]);
		// phenotype IDs
		phenotypeIDs = new ArrayList<String>(Arrays.asList(lines[10].split(",")));
		// Origin
		final ArrayList<String> checkOrigin = new ArrayList<String>(Arrays.asList("germline", "biparental", "de novo",
				"maternal", "paternal", "somatic", "inherited", "unknown", "not provided", "tested-inconclusive"));
		origin = new ArrayList<String>(Arrays.asList(lines[11].split(";")));
		for (int i = 0; i < origin.size(); i++) {
			if (!checkOrigin.contains(origin.get(i)))
				throw new IllegalArgumentException("Illegal origin: " + lines[11] + "\t AlleleID: " + lines[0]);
		}
		// assembly
		assembly = "GRCh37";
		if(useGrch38)
			assembly = "GRCh38";
		if (!assembly.equals(lines[12]))
			throw new IllegalArgumentException("Illegal reference genome assembly: " + lines[12] + "\t AlleleID: " + lines[0]);
		// chromosome
		setChromosome(lines[13]);
		// start
		start =Integer.parseInt(lines[14]);
		// end
		stop = Integer.parseInt(lines[15]);
		// chromosome band
		try {
			setCytocinetic(lines[16]);
		} catch (IllegalArgumentException e) {
			Log.add(e.getMessage(), 1);
		}
		// review status
		final ArrayList<String> checkReviewStatus = new ArrayList<String>(Arrays.asList("expert panel", "reviewed by expert panel",
				"professional society", "multiple submitters", "single submitter", "criteria provided", "no conflicts",
				"practice guideline", "conflicting interpretations", "no assertion criteria provided"));
		String[] reviewStatusTemp = lines[17].split(", ");
		for (int i = 0; i < reviewStatusTemp.length; i++) {
			if (checkReviewStatus.contains(reviewStatusTemp[i])) {
				reviewStatus.add(reviewStatusTemp[i]);
			}
		}
		if(reviewStatus.size() < 1)
			throw new IllegalArgumentException("Illegal review status: " + lines[17] + "\t AlleleID: " + lines[0]);
		// HGVS for DNA
		dnaHgvs = lines[18];
		// HGVS for Protein is ignored
		// number of submitters
		numberSubmitters =Integer.parseInt(lines[20]);
		// date of last evaluation
		lastEvaluated =lines[21];
		// Guideline is ignored
		// Other IDs
		otherIDs =new ArrayList<String>(Arrays.asList(lines[23].split(",")));
		variantID =Integer.parseInt(lines[24]);
		// Parse additional Information from File
		parseInformation();
		Log.add("ClinVariant with AlleleID " + alleleId + " added.", 1);
	}

	/**
	 * Validate and set cytogenetic
	 * @param chrBand chromosome number and band on the chromosome
	 */
	private void setCytocinetic(String chrBand) throws IllegalArgumentException {
		if (!Pattern.matches("^"+ chromosome +"[pq][1-4][0-9](\\.[1-9][0-9]?(\\.[1-9])?)?"
				+ "(-(("+ chromosome +")?[pq])?[1-4][0-9](\\.[1-9][0-9]?(\\.[1-9])?)?)?$",chrBand))
			throw new IllegalArgumentException("Illegal chromosome band:\t" + chrBand + "\t AlleleID: " + alleleId);
		cytogenetic =chrBand;
	}

	/**
	 * Parse the Information of the Variant in other fields.
	 */
	private void parseInformation() throws DebuggingException, VariantToFarAwayException{
		// refSeq
		String refSeqWithVersion = dnaHgvs.split(":", 2)[0];
		setRefSeqAccession(refSeqWithVersion);
		// refSeqVersion
		String refSeqVersion = refSeqWithVersion.split("\\.", 2)[1];
		setRefSeqVersion(Integer.parseInt(refSeqVersion));
		// beforeExon
		acceptor = dnaHgvs.contains("-");
		// distanceToExon
		int end = dnaHgvs.length()-3;
		int begin;
		if(acceptor){
			begin = dnaHgvs.indexOf('-');
		}
		else{
			begin = dnaHgvs.indexOf('+');
		}
		if(begin == -1)
			throw new DebuggingException("Can't find + or - in DNA HGVS\t AlleleID: " + alleleId);
		if(start != stop)
			throw new DebuggingException("start != stop\t AlleleID: " + alleleId);
		distanceToExon = Integer.parseInt(dnaHgvs.substring(begin, end));
		if(Math.abs(distanceToExon) > Config.getLengthTrainingIntron())
			throw new VariantToFarAwayException("The variant with HGVS " + dnaHgvs + " is to far away from junction: " + Math.abs(distanceToExon));
		// alt
		if(!dnaHgvs.contains(">")){
			throw new DebuggingException("DNA change contains no \">\"\t AlleleID: " + alleleId);
		}
		String chgvs[] = dnaHgvs.split(">");
		alt = chgvs[1];
		if (!Functions.isValidDNA(alt)) {
			throw new DebuggingException("DNA change contains invalid DNA:" + alt + "\t AlleleID: " + alleleId);
		}
		// ref
		chgvs = chgvs[0].split("[0-9]");
		ref = chgvs[chgvs.length-1];
		if (!Functions.isValidDNA(ref)) {
			throw new DebuggingException("DNA change contains invalid DNA:" + ref + "\t AlleleID: " + alleleId);
		}
	}

//	/**
//	 * returns all Annotations as tab-separated String
//	 */
//	// @override
//	public String toString() {
//		return alleleID + '\t' + type + '\t' + name + '\t' + geneID + '\t' + gene + '\t' + clinicalSignificance
//				+ '\t' + snp + '\t' + nsv + '\t' + rcv + '\t' + testedInGTR + '\t' + phenotypeIDs + '\t' + origin + '\t'
//				+ assembly + '\t' + chromosome + '\t' + start + '\t' + stop + '\t' + cytogenetic + '\t' + reviewStatus
//				+ '\t' + dnaHgvs + '\t' + cLRG + '\t' + numberSubmitters + '\t' + lastEvaluated + '\t' + otherIDs + '\t'
//				+ variantID + '\t' + refSeqAccession + '\t' + refSeqVersion + '\t' + ref + '\t' + alt;
//	}

	public int getAlleleId() {
		return alleleId;
	}

	public String getName() {
		return name;
	}

	/**
	 * List of values of clinical significance reported for this variation
	 */
	public ArrayList<String> getClinicalSignificance() {
		return clinicalSignificance;
	}

}
