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
public class LogThread {

  private static boolean writing = false;
  /**
   * The main entries of the Log
   */
  private ArrayList<String> logEntries = new ArrayList<String>();
  /**
   * Every log entry has a priority
   */
  private ArrayList<Integer> entryPriorities = new ArrayList<Integer>();
  private static boolean fileExists = false;
  private static FileWriter fWriter = null;

  /**
   * 
   */
  public LogThread() {
    super();
  }

  /**
   * Write log to File and to console
   */
  public void writeToFile() {
    while (writing) {
      try {
        Thread.sleep(100);
      } catch (InterruptedException e) {
        e.printStackTrace();
      }
    }
    writing = true;
    try {
      String title = "";
      if (!fileExists) {
        File file = new File(Config.getLogFile());
        // delete old file and create file
        if (file.exists()){
          file.delete();
        }
        fileExists = file.createNewFile();
        fileExists = true;
        title = "\t\t---   jsplice log File   ---";
        fWriter = new FileWriter(file);
        file.setReadable(true);
      }
      if (fileExists) {
        fWriter.write(title);
        // write variants lines to file
        for (int i=0; i < logEntries.size(); i++) {
          if (entryPriorities.get(i) >= Log.getFileLogLevel()) {
            fWriter.append("\n");
            fWriter.write("[" + entryPriorities.get(i) + "] " + logEntries.get(i));
          }
        }
      } else
        System.out.println("Unable to create file.");
    } catch (IOException e) {
      e.printStackTrace();
    }
    System.out.println(this);
    logEntries = new ArrayList<String>();
    entryPriorities = new ArrayList<Integer>();
    writing = false;
  }

  /**
   * reset the static object for test cases
   */
  public void reset() {
    logEntries = new ArrayList<String>();
    entryPriorities = new ArrayList<Integer>();
  }

  public static void close() {
    if (fWriter != null){
      try {
        fWriter.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
  }

  /**
   * Add an error string
   * @param message 
   * @param priority int between 1 and 6: <br/>
   * 1 Trace, 2 Debug, 3 Info, 4 Warn, 5 Error, 6 Fatal
   */
  public void add(Object message, int priority){
    if(priority > 6 || priority < 1)
      throw new IllegalArgumentException("Priority level has to be between 1 and 6: " + priority);
    logEntries.add(message.toString());
    entryPriorities.add(priority);
  }

  @Override
  public String toString(){
    String str = "";
    for (int i = 0; i < logEntries.size(); i++) {
      if(entryPriorities.get(i) >= Log.getPrintLogLevel())
        str += logEntries.get(i) + "\n";
    }
    return str;
  }

  public int size() {
    return logEntries.size();
  }
}
