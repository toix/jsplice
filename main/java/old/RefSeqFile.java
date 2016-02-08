package old;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

import jsplice.data.Variant;

public class RefSeqFile {
	
	/**
	 * File path of the fasta-file where the RefSeq sequences are stored
	*/
	private String file;
	private String[] header;
	private HashMap<String, Character> strand = new HashMap<String, Character>();
	private HashMap<String, Integer> chromosome = new HashMap<String, Integer>();
	private HashMap<String, Integer> start = new HashMap<String, Integer>();
	private HashMap<String, Integer> stop = new HashMap<String, Integer>();
	private ArrayList<Variant> variants;

	public RefSeqFile(String file) {
		this.file = file;
	}
	
//	public void extendVariantByStrand(ArrayList<Variant> variants) {
//		this.variants = variants;
//		readFile();
//		addInformation();
//	}
	
	private void readFile(){
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
			inputStream = new FileInputStream(file);
			sc = new Scanner(inputStream, "UTF-8");
			while (sc.hasNextLine()) {
				String line = sc.nextLine();
				if(line.charAt(0) == '>'){
					parseInformation(line);
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
	 * Create Hash with RefSec mRNA accession numbers for for relevant Information
	 * @param line is the meta information line of Sequence in the fa-file
	 * @param lineNr is the line number in the RefSeq sequence file
	*/
	private void parseInformation(String line){
		String[] lineArray =line.split(" ");
		// RefSeqID
		String refSeqID = lineArray[0].split("_", 3)[2];
		if (refSeqID.matches("^NM_[0-9]{6,9}$")) {
			// chromosome
			String[] position = lineArray[1].split(":");
			chromosome.put(refSeqID, Character.getNumericValue((position[0].charAt(position[0].length() - 1))));
			// start
			String[] range = position[1].split("-");
			start.put(refSeqID, Integer.parseInt(range[0]));
			// stop
			stop.put(refSeqID, Integer.parseInt(range[1]));
			// strand
			strand.put(refSeqID, lineArray[4].charAt(lineArray[4].length() - 1));
			// file path
		}
	}

//	private void addInformation() {
//		for(int i=0; i<variants.size(); i++){
//			String refSeqID = variants.get(i).getRefSeqAccession();
//			// add Strand orientation to variant
//			if(strand.containsKey(refSeqID)){
//				variants.get(i).setPlusStrand(strand.get(refSeqID));
//				variants.get(i).setTranscriptStart(start.get(refSeqID));
//				variants.get(i).setTranscriptStop(stop.get(refSeqID));
//				variants.get(i).setRefSeqFileName(file);
//			} 
//			// or delete variant
//			else{
//				if(variants.remove(variants.get(i)))
//					System.out.println("Removed variant:\n"+ variants.get(i));
//				else
//					System.out.println("Can't remove variant:\n"+ variants.get(i));
//			}
//		}
//	}
	
	public String getTranscriptSequence(String refSeqID){
		String transcript = null;
		FileInputStream inputStream = null;
		Scanner sc = null;
		try {
			inputStream = new FileInputStream(file);
			sc = new Scanner(inputStream, "UTF-8");
			sc.useDelimiter(">");
			int lineNr=0;
			boolean finished = false;
			while (sc.hasNextLine() && !finished) {
				String line = sc.nextLine();
				// Sequence to RefSeqID found?
				if (line.matches("^>hg[1-9][0-9]_refGene_" + refSeqID + ".*")) {
					transcript = "";
//					System.out.println("found in line "+i+":\t"+ line);
				} else if (transcript != null) {
					if (line.charAt(0) == '>')
						finished = true;
					else{
						transcript = transcript.concat(line);
					}
				}
				lineNr++;
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
		return transcript;
	}
}
