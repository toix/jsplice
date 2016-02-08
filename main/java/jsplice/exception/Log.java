/**
 * 
 */
package jsplice.exception;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import jsplice.data.Config;

/**
 * Stores all log entries <br/> 
 * <br/>
 * - 1 Trace - Only when I would be "tracing" the code and trying to find one
 * part of a function specifically <br/>
 * - 2 Debug - Information that is diagnostically helpful to people more than
 * just developers (IT, sysadmins, etc)<br/>
 * - 3 Info - Generally useful information to log (service start/stop,
 * configuration assumptions, etc). Info I want to always have available but
 * usually dont care about under normal circumstances. This is my
 * out-of-the-box config level<br/>
 * - 4 Warn - Anything that can potentially cause application oddities, but for
 * which I am automatically recoverring (such as switching from a primary to
 * backup server, retrying an operation, missing secondary data, etc)<br/>
 * - 5 Error - Any error which is fatal to the operation but not the service or
 * application (cant open a required file, missing data, etc). These errors
 * will force user (administrator, or direct user) intervention. These are
 * usually reserved (in my apps) for incorrect connection strings, missing
 * services, etc.<br/>
 * - 6 Fatal - Any error that is forcing a shutdown of the service or
 * application to prevent data loss (or further data loss). I reserve these
 * only for the most heinous errors and situations where there is guaranteed
 * to have been data corruption or loss.<br/>
 * @author Tobias Gresser (gresserT@gmail.com)
 */
public class Log {
	
	/**
	 * The main entries of the Log
	 */
	private static ArrayList<String> logEntries = new ArrayList<String>();
	/**
	 * Every log entry has a priority
	 */
	private static ArrayList<Integer> entryPriorities = new ArrayList<Integer>();
	/**
	 * All entries higher or equal to the priority level will be printed
	 */
	private static int printLogLevel = 3;
	/**
	 * All entries higher or equal to the priority level will be written to the log file
	 */
	private static int fileLogLevel = 2;

	/**
	 * 
	 */
	public Log() {
		super();
	}
	
	/**
	 * Add an error string with priority 6
	 * @param message 
	 */
	public static void add(String message){
		logEntries.add(message);
		entryPriorities.add(6);
	}
	
	/**
	 * Add an error string
	 * @param message 
	 * @param priority int between 1 and 6: <br/>
	 * 1 Trace, 2 Debug, 3 Info, 4 Warn, 5 Error, 6 Fatal
	 */
	public static void add(String message, int priority){
		if(priority > 6 || priority < 1)
			throw new IllegalArgumentException("Priority level has to be between 1 and 6: " + priority);
		logEntries.add(message);
		entryPriorities.add(priority);
	}

	public static String toStringStatic(){
		String str = "";
		for (int i = 0; i < logEntries.size(); i++) {
			if(entryPriorities.get(i) >= printLogLevel)
			str += logEntries.get(i) + "\n";
		}
		return str;
	}
	
	@Override
	public String toString(){
		return Log.toStringStatic();
	}
	
	/**
	 * Write log to File
	 */
	public static void writeToFile() {
		File file = new File(Config.getLogFile());
		FileWriter fWriter = null;
		try {
			// delete old file and create file
			if (file.exists())
				file.delete();
			if (file.createNewFile()) {
				file.setReadable(true);
				fWriter = new FileWriter(file);
				String title = "\t\t---   jsplice log File   ---";
				fWriter.write(title);
				fWriter.append("\n");
				// write variants lines to file
				for (int i=0; i < logEntries.size(); i++) {
					if (entryPriorities.get(i) >= fileLogLevel) {
						fWriter.write(logEntries.get(i));
						fWriter.append("\n");
					}
				}
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
	
	public static int getPrintLogLevel() {
		return printLogLevel;
	}
	
	public static int getFileLogLevel() {
		return fileLogLevel;
	}

	/**
	 * 
	 * @param logLevel int between 1 and 6: Trace, Debug, Info, Warn, Error, Fatal
	 */
	public static void setPrintLogLevel(int logLevel) {
		if(logLevel > 6 || logLevel < 1)
			throw new IllegalArgumentException("Priority level has to be between 1 and 6: " + logLevel);
		Log.printLogLevel = logLevel;
	}
	
	/**
	 * 
	 * @param logLevel int between 1 and 6: Trace, Debug, Info, Warn, Error, Fatal
	 */
	public static void setFileLogLevel(int logLevel) {
		if(logLevel > 6 || logLevel < 1)
			throw new IllegalArgumentException("Priority level has to be between 1 and 6: " + logLevel);
		Log.fileLogLevel = logLevel;
	}
	
	public static int size() {
		return logEntries.size();
	}
	
	/**
	 * reset the static object for test cases
	 */
	public static void reset() {
		logEntries = new ArrayList<String>();
		entryPriorities = new ArrayList<Integer>();
		printLogLevel = 3;
		fileLogLevel = 2;
	}
}
