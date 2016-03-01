package jsplice.tools;

/**
 * TODO initalization for static functions
 * General Array and Sequence Functions <br/>
 */

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import jsplice.data.Sequence;
import jsplice.data.Variant;

import com.beust.jcommander.ParameterException;

public class Functions {
	
//	public final static char[] bases = {'A','T','G','C', 'a', 'c', 'g', 't', 'N'};#
	public static final String fastaBases = new String("ACGTacgtN");
	public static final String bases = new String("ACGT");
	public static final int[] baseNumber ={0,1,2,3,0,1,2,3,4};
	/**
	 * Hashes from Base to its corresponding number <br/>
	 * 0=A, 1=T, 2=G, 3=C, 4=N
	 */
	public static HashMap<Character, Integer> mapNumber = new HashMap<Character, Integer>();
	/**
	 * Flag whether the Class has been initialized.
	 */
	private static boolean initialized;

	public Functions() {
		for (int i = 0; i < fastaBases.length(); i++) {
			mapNumber.put(fastaBases.charAt(i), baseNumber[i]);
		}
		initialized = true;
	}
	// TODO useless?
	public static double[] getInitializedDoubleArray(int length) {
		double[] array = new double[length];
		for (int i = 0; i < length; i++) {
			array[i] = 0;
		}
		return array;
	}
	
	public static int[] getInitializedIntArray(int length) {
		int[] array = new int[length];
		for (int i = 0; i < length; i++) {
			array[i] = 0;
		}
		return array;
	}
	
	/**
	 * Will throws an {@link ParameterException} if the sequence contains invalid chars
	 * @param sequence May contain the char A, C, G, T, a, c, g, t and N
	 * @return 
	 */
	public static boolean isValidDNA(String sequence){
		if (java.util.regex.Pattern.matches("[" + fastaBases + "]*", sequence) && sequence != null) {
			return true;
		} else {
			return false;
		}
//		for (int i = 0; i < sequence.length(); i++) {
//			if(!fastaBases.contains(""+ sequence.charAt(i))){
//				return false;
//			}
//		}
//		return true;
	}

	/**
	 * Complements and reverses DNA sequence
	 * @param sequence May contain the Character A, C, G, T, a, c, g, t and N
	 * @return complement, reverse sequence
	 */
	public static String reverseAndComplement(String sequence){
		return reverse(complement(sequence));
	}
	
	/**
	 * Complements the DNA sequence
	 * @param sequence May contain the Character A, C, G, T, a, c, g, t and N
	 * @return complement sequence
	 */
	public static String complement(String sequence) {
		if (!isValidDNA(sequence)) {
			throw new IllegalArgumentException("The sequence contains invalid DNA:\n" + sequence);
		}
		String complement = "";
		for (int i = 0; i < sequence.length(); i++) {
			switch (sequence.charAt(i)) {
			case 'A':
				complement += "T";
				break;
			case 'C':
				complement += "G";
				break;
			case 'G':
				complement += "C";
				break;
			case 'T':
				complement += "A";
				break;
			case 'a':
				complement += "t";
				break;
			case 'c':
				complement += "g";
				break;
			case 'g':
				complement += "c";
				break;
			case 't':
				complement += "a";
				break;
			case 'N':
				complement += "N";
				break;
			default:
				throw new ParameterException("This should be impossible");
			}
		}
		return complement;
	}

	/**
	 * Reverses the DNA sequence
	 * @param sequence May contain the Character A, C, G, T, a, c, g, t and N
	 * @return reverse sequence
	 */
	public static String reverse(String sequence){
		if (!isValidDNA(sequence)) {
			throw new IllegalArgumentException("The sequence contains invalid DNA:\n" + sequence);
		}
		return new StringBuilder(sequence).reverse().toString();
	}

	public static String[] transpose(String[] array) {
		int numberOfLines = array.length;
		int lineLength = array[0].length();
		String[] transposed = new String[lineLength];
		for (int i = 0; i < lineLength; i++) {
			transposed[i] = "";
			for (int j = 0; j < numberOfLines; j++) {
				transposed[i] += array[j].charAt(i);
			}
		}
		return transposed;
	}
	
	/**
	 * @param sequences
	 * @return String array where every String represents a sequence position
	 */
	public static String[] transpose(ArrayList<Sequence> sequences, boolean reference) {
		int numberOfLines = sequences.size();
		int lineLength = sequences.get(0).length();
		String[] transposed = new String[lineLength];
		for (int i = 0; i < lineLength; i++) {
			transposed[i] = "";
			for (int j = 0; j < numberOfLines; j++) {
				Sequence sequence = sequences.get(j);
				if(!reference && i == sequence.getPositionChange()){
					transposed[i] += sequence.getAlt();
				} else {
					transposed[i] += sequence.charAt(i);
				}
			}
		}
		return transposed;
	}

	public static String arrayToString(String[] array){
		String arrayString = new String(array[0]);
		for (int i = 1; i < array.length; i++) {
			arrayString += ", " + array[i];
		}
		return arrayString;
	}
	
	public static String arrayToString(double[] array, int decimalNumbers){
		double rounded = round(array[0], 4);
		String arrayString = new String(limitNumber(rounded) +"");
		for (int i = 1; i < array.length; i++) {
			rounded = round(array[i], 4);
			arrayString += "\t " + limitNumber(rounded);
		}
		return arrayString;
	}
	
	public static String arrayToString(int[] correlationCluster){
      double rounded = round(correlationCluster[0], 4);
      String arrayString = new String(limitNumber(rounded) +"");
      for (int i = 1; i < correlationCluster.length; i++) {
          rounded = round(correlationCluster[i], 4);
          arrayString += "\t " + limitNumber(rounded);
      }
      return arrayString;
  }
	
	/**
	 * @param number
	 * @return Number between -999 and 999
	 */
	public static String limitNumber(double number){
		String smallNumber;
		if (number > 999) {
			smallNumber = 999 + "";
		} else if (number < -999) {
			smallNumber = -999 + "";
		} else {
			smallNumber = number + "";
		}
		return smallNumber;
	}

	/**
	 * Converts the sequence to a sequence containing solely big DNA letters.
	 * @param sequence A DNA sequence containing only "aAtTcCgG".
	 * @return
	 */
	public static String convertToUpperBases(String sequence){
		Functions.isValidDNA(sequence);
		sequence = sequence.toUpperCase();
		for (int i = 0; i < sequence.length(); i++) {
			if (!Functions.bases.contains(sequence.substring(i, i+1))) {
				throw new IllegalArgumentException("The sequence contains an invalid letter: " + sequence.charAt(i));
			}
		}
		return sequence;
	}

	/**
	 * The sum of an array
	 * @param summmands
	 * @return sum
	 */
	public static double sum(double[] summmands) {
		double sum = 0;
		for (int i = 0; i < summmands.length; i++) {
			sum += summmands[i];
		}
		return sum;
	}
	
	public static int sum(ArrayList<Integer> list) {
	  int sum = 0;
	  for (int i = 0; i < list.size(); i++) {
	    sum += list.get(i);
	  }
	  return sum;
	}

	/**
	 * Probability that a base occurs at a certain position
	 * @return
	 */
	public static double[] getFrequencies(String sequence){
		if (!initialized) {
			for (int i = 0; i < fastaBases.length(); i++) {
				mapNumber.put(fastaBases.charAt(i), baseNumber[i]);
			}
			initialized = true;
		}
		double[] count = getInitializedDoubleArray(bases.length());
		double[] prob = new double[bases.length()];
		for (int i = 0; i < sequence.length(); i++) {
			char base =sequence.charAt(i);
			Integer baseIdx = mapNumber.get(base);
			try {
				count[baseIdx]++;
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new IllegalArgumentException();
			}
		}
		for (int i = 0; i < bases.length(); i++) {
			if (count[i] != 0) {
				prob[i] = count[i] / sequence.length();
			} else {
				prob[i] = 1 / (sequence.length() + 2);
			}
		}
		return prob;
	}
	
	/**
	 * Probability that a base occurs at a certain position on multiple positions
	 * @return
	 */
	public static double[][] getFrequencies(ArrayList<Sequence> sequences, boolean reference) {
		int sequenceLength = sequences.get(0).length();
		double[][] probability = new double[sequenceLength][bases.length()];
		String[] sequenceRows = Functions.transpose(sequences, reference);
		for (int i = 0; i < sequenceLength; i++) {
			probability[i] = getFrequencies(sequenceRows[i]);
		}
		return probability;
	}

	/**
	 * round function
	 */
	public static double round(double number, int decimalPos){
		double multiplier = Math.pow(10, decimalPos);
		return (Math.round(number*multiplier)/multiplier);
	}
	
	public static boolean equals(double[] array1, double[] array2){
		boolean equals = true;
		if(array1.length != array2.length){
			throw new IllegalArgumentException("Both arrays must have equal lenght.");
		}
		for (int i = 0; i < array1.length; i++) {
			equals = equals && array1[i] == array2[i];
		}
		return equals;
	}
	
	/**
	 * Write String to File
	 * @param text desired file content
	 * @param fileName name, file extension and directory of the file
	 */
	public static void writeToFile(String text, String fileName) {
		File file = new File(fileName);
		FileWriter fWriter = null;
		try {
			new File(fileName).mkdirs();
			// delete old file and create file
			if (file.exists())
				file.delete();
			if (file.createNewFile()) {
				file.setReadable(true);
				fWriter = new FileWriter(file);
				// write variants lines to file
				fWriter.write(text);
			} else
				System.out.println("Unable to create file.");
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (fWriter != null){
				try {
					fWriter.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * @param from included
	 * @param to included
	 * @return
	 */
	public static int random(int from, int to) {
		return (int) (Functions.round(Math.random() * (to - from + 1) + from -0.5,0));
	}
	
	public static int countMatchingChar(String str1, String str2) {
		if (str1.length() != str2.length()) {
			throw new IllegalArgumentException("Both strings must have same length.");
		}
		int count = 0;
		for (int p = 0; p < str1.length(); p++) {
			if (str1.charAt(p) == str2.charAt(p)) {
				count++;
			}
		}
		return count;
	}

  /**
   * @param matrix1
   * @param matrix2
   * @return
   */
  public static double probabilityMatrixSimilarity(double[][] matrix1, double[][] matrix2) {
    if (matrix1.length != matrix2.length) {
      throw new IllegalArgumentException("Both strings must have same length.");
    }
    double similarity = 0;
    for (int l = 0; l < matrix1.length; l++) {
      double similaritySum = 0;
        for (int b = 0; b < matrix1[l].length; b++) {
          double sumMatrix = matrix1[l][b] + matrix2[l][b];
          if (sumMatrix > 0) {
            similaritySum += sumMatrix;
          }
        }
      // System.out.println(Functions.arrayToString(matrix1[p], 1));
//      System.out.println(Functions.arrayToString(matrix2[p], 1));
//      System.out.println("add: " + distanceSum / matrix1[p].length);
      similarity += similaritySum / matrix1[l].length;
    }
    return similarity;
  }
  /**
   * @param list
   */
  public static double mean(ArrayList<Integer> list) {
    if (list.size() == 0) {
      throw new IllegalArgumentException();
    }
    double count = sum(list);
    return count / list.size();
  }
}

