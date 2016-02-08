package jsplice.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import jsplice.data.Variant;

public class VcfFile {
	private String file;
	private String[] header;
	private ArrayList<Variant> variants;

	public VcfFile(String file) {
		this.file =file;
	}
	
	public void writeVariantFile(ArrayList<Variant> variants) throws IOException{
		this.variants = variants;
		createFile();
		BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		writeHeader(writer);
		
		for(Variant variant : variants){
			writer.write(variant.getLine());
			writer.newLine();
		}
	}
	
	private void createFile() throws IOException{
		File osFile = new File(file);
		if(osFile.exists()){
			osFile.delete();
		}
		osFile.createNewFile();
	}
	
	private void writeHeader(BufferedWriter writer) throws IOException{
		writer.write("#CHRO\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1");
		writer.newLine();
	}
}