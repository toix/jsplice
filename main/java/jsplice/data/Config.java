package jsplice.data;

import jsplice.tools.Functions;

import com.beust.jcommander.ParameterException;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class Config {

  private static Config instance;
  private static String clinVarFileName = "../data/variant_summary.txt";
  private static String hgmdFileName = "../data/hgmd.txt";
  private static String vcfFileName = "../data/bowtie.vcf";
  private static String chrFileNames = "../data/chromFa/chr*.fa";
  private static String refGeneFileName = "../data/RefSeqGenesGRCh37.tsv";
  /**
   * Length of the whole sequence used for the algorithm. <br/>
   * 2 * sequenceLength defines the length of the training data <br/>
   * Has to be even. Is the sum of modelIntronLength and modelExonLength
   */
  private static int lengthModel = 20;
  /**
   * Length of the exonic part of the sequence. <br/>
   * sequenceLength - exonLength defines the intron length
   */
  private static int lengthModelExon = 2;
  /**
   * Multiplier by that the training data is longer than the analysis length (sequenceLength) <br/>
   * Has to be at least 3
   */
  private static int factor = 5;
  /**
   * false -> GRCh37 will be used <br/>
   * true -> GRCh38 will be used
   */
  private static boolean useGrch38 = false;
  public static final String[] chromosomeNames = {"1", "2", "3", "4", "5", "6", "7", "8", "9",
      "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"};
  public static final double quantityRelLimit = 1;
  public static final double ClusterCorrelationLimit = 0.5;
  /**
   * The maximum length of the pattern for clustering
   */
  public static final int lengthIntronPatternMax = 12;
  public static final int distanceClusterMinAcceptor = 4;
  public static final int distanceClusterMinDonor = 6;
  public static final int distanceClusterMax = getLengthTrainingIntron() - lengthIntronPatternMax;
  public static final boolean multiClusterRel = true;
  public static final boolean simpleMerging = false;
  public static final double mergingCorrelationMin = 2.0;
  public static final boolean clusteringSub = true;
  public static final boolean clusteringBasic = false;
  public static final String folder = "results/cluster"
      + (multiClusterRel ? "Rel" : "Std")
      + (clusteringSub ? "Csub" : "")
      + (clusteringBasic ? "Cbasic" : "")
      + (simpleMerging ? "Spl" : ("Cpx" + mergingCorrelationMin))
      + "/" + Config.getLengthModelIntron() + "+" + Config.getLengthModelExon() + "/";
  private static String logFile = folder + "jsplice.log";
  public static final int crossValidationSteps = 1000;
	
	/**
	 * 
	 * @param clinVarFileName
	 * @param vcfFileName
	 * @param chrFileNames
	 * @param refGeneFileName
	 * @param sequenceLength Length of the whole Sequence. Has to be even.
	 * @param useGrch38
	 * @throws ParameterException
	 */
	public Config(String clinVarFileName, String hgmdFileName, String vcfFileName,
			String chrFileNames, String refGeneFileName, String logFile, int sequenceLength,
			int exonLength, boolean useGrch38) throws ParameterException {
		if (instance == null) {
			if (exonLength/2 >= sequenceLength) {
				throw new ParameterException("\""+ exonLength +"\" is an invalid exon length. It has to shorter than sequence length.");
			}
			Config.clinVarFileName = clinVarFileName;
			Config.hgmdFileName = hgmdFileName;
			Config.vcfFileName = vcfFileName;
			Config.chrFileNames = chrFileNames;
			Config.refGeneFileName = refGeneFileName;
			Config.logFile = logFile;
			Config.lengthModel = sequenceLength;
			Config.lengthModelExon = exonLength;
			Config.useGrch38 = useGrch38;
			Config.instance = this;
		} else
			throw new ParameterException(
					"An Instance of class GlobalParameters already exists");
	}
	public Config() throws ParameterException {
		if (instance == null) {
			Config.instance = this;
		} else
			throw new ParameterException(
					"An Instance of class GlobalParameters already exists");
	}

	//	/**
	//	 * @return 
	//	 */
	//	public static double getThresholdAltRef() {
	//		return thresholdAltRef;
	//	}
		/**
		 * @param lengthIntron
		 */
		public static void setLengthModelIntron(int lengthIntron) {
			lengthModel = lengthModel - getLengthModelIntron() + lengthIntron;
			lengthModelExon = lengthModel - lengthIntron;
		}
	public static String[] getChromosomenames() {
		return chromosomeNames;
	}

	public static Config getInstance() {
		return instance;
	}

	public static boolean isInitialized() {
	
		return false;
	}
	public static String getClinVarFileName() {
		return clinVarFileName;
	}
	
	public static String getHgmdFileName() {
		return hgmdFileName;
	}

	public static String getVcfFileName() {
		return vcfFileName;
	}

	public static String getChrFileNames() {
		return chrFileNames;
	}

	public static String getRefGeneFileName() {
		return refGeneFileName;
	}

	public static String getLogFile() {
		return logFile;
	}
	
	public static int getLengthModelSequence() {
		return lengthModel;
	}
	public static int getLengthModelExon() {
		return lengthModelExon;
	}
	public static int getLengthModelIntron(){
		return getLengthModelSequence()-getLengthModelExon();
	}
	public static int getLengthTrainingSequence(){
		return (int) Math.round(lengthModel*factor + getLengthModelIntron());
	}
	public static int getLengthTrainingExon(){
		return (int) Math.round(lengthModelExon*factor + getLengthModelIntron());
	}
	public static int getLengthTrainingIntron(){
		return Math.round(getLengthTrainingSequence()-getLengthTrainingExon());
	}

//	/**
//	 * Calculates the start position of the analysis by junction type and junction position
//	 * @param junctionPosition Position of the exon junction
//	 * @param intronExonJunction Defines whether the sequence belongs to an intron exon junction
//	 * @return Position in the sequence
//	 */
//	public static int getAnalysisStartPosition(int junctionPosition, boolean intronExonJunction){
//		return junctionPosition - getMinAnalysisPosition(intronExonJunction);
//	}
//	/**
//	 * Calculates the end position of the analysis by junction type and junction position
//	 * @param junctionPosition Position of the exon junction
//	 * @param intronExonJunction Defines whether the sequence belongs to an intron exon junction
//	 * @return Position in the sequence
//	 */
//	public static int getAnalysisEndPosition(int junctionPosition, boolean intronExonJunction){
//		return junctionPosition - getMaxAnalysisPosition(intronExonJunction);
//	}
	
	/**
	 * Calculates the minimum possible position that can be analyzed
	 * @param acceptorSequence Defines whether the sequence belongs to an acceptor junction
	 * @return Position in the sequence
	 */
	public static int getMinAnalysisPosition(boolean acceptorSequence, boolean acceptorModel){
		if (acceptorSequence && acceptorModel) {
//			System.out.println("-  (" + getTrainingIntronLength()+ " - " + getModelIntronLength() + ") + " + getModelIntronLength());
			return - (getLengthTrainingIntron() - getLengthModelIntron()) + getLengthModelIntron();
//			return - (factor - 1) * getModelIntronLength() + getModelIntronLength();
		} else if (acceptorSequence && !acceptorModel) {
//			return - (factor - 1) * getModelIntronLength() + modelExonLength;
			return - (getLengthTrainingIntron() - getLengthModelIntron()) + getLengthModelExon();
		} else if (!acceptorSequence && acceptorModel) {
//			return - (factor - 1) * modelExonLength;
			return - (getLengthTrainingExon() - getLengthModelExon()) + getLengthModelIntron();
		} else {
//			return - (factor - 2) * modelExonLength - getModelIntronLength();
			return - (getLengthTrainingExon() - getLengthModelExon()) + getLengthModelExon();
		}
	}
	/**
	 * Calculates the maximum possible position that can be analyzed in a pattern sequence
	 * @param acceptorSequence Has to be true if the sequence belongs to an intron exon junction
	 * @param sequenceLength Length of the pattern sequence
	 * @return Position in the sequence
	 */
	public static int getMaxAnalysisPosition(boolean acceptorSequence, boolean acceptorModel){
		if (acceptorSequence && acceptorModel) {
//			return modelSequenceLength-1 + (factor - 2) * modelExonLength + getModelIntronLength();
			return getLengthTrainingSequence()-1 - (getLengthTrainingIntron() - getLengthModelIntron()) - getLengthModelExon();
		} else if (acceptorSequence && !acceptorModel) {
//			System.out.println(getTrainingSequenceLength() + "-1 - (" + getTrainingIntronLength() + " - " + getModelIntronLength() + ") - " + getModelIntronLength());
//			return modelSequenceLength-1 + (factor - 2) * modelExonLength;
			return getLengthTrainingSequence()-1 - (getLengthTrainingIntron() - getLengthModelIntron()) - getLengthModelIntron();
		} else if (!acceptorSequence && acceptorModel) {
//			return modelSequenceLength-1 + (factor - 2) * getModelIntronLength();
			return getLengthTrainingSequence()-1 - (getLengthTrainingExon() - getLengthModelExon()) - getLengthModelExon();
		} else {
//			return modelSequenceLength-1 + (factor - 2) * getModelIntronLength();
			return getLengthTrainingSequence()-1 - (getLengthTrainingExon() - getLengthModelExon()) - getLengthModelIntron();
		}
	}
	public static boolean useGrch38() {
		return useGrch38;
	}
	
	/**
	 * @param acceptorP
	 * @return The length of the Intron for the cryptic site prediction
	 */
	public static int getLengthIntronCryptic(boolean acceptorP) {
		return (int) Functions.round((getDistanceClusterMin(acceptorP) + 0), 0);
	}
	/**
	 * @param acceptorP 
	 * @return
	 */
	public static int getDistanceClusterMin(boolean acceptorP) {
		return acceptorP? distanceClusterMinAcceptor : distanceClusterMinDonor;
	}
	/**
	 * @param acceptorP 
	 * @return
	 */
	public static int getDistanceClusterMax() {
		return distanceClusterMax;
	}
}
