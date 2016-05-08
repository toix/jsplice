/**
 * 
 */
package jsplice.exception;

import java.util.HashMap;

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

  private static HashMap<Thread, LogThread> logThread = new HashMap<Thread, LogThread>();
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
    getLog().add("Log file: " + Config.getLogFile(), 3);
  }

  /**
   * @return
   */
  private static LogThread getLog() {
    if (!logThread.containsKey(Thread.currentThread())) {
      logThread.put(Thread.currentThread(), new LogThread());
    }
    return logThread.get(Thread.currentThread());
  }

  /**
   * Write log to File and to console
   */
  public static void writeToFile() {
    getLog().writeToFile();
    logThread.remove(getLog());
  }
  
  /**
   * reset log for thread
   */
  public static void remove() {
    logThread.remove(getLog());
  }
  
  /**
   * reset the static object for test cases
   */
  public void reset() {
    getLog().reset();
    printLogLevel = 3;
    fileLogLevel = 2;
  }

  /**
   * 
   */
  public static void close() {
    LogThread.close();
  }

  /**
   * Add an error string with priority 6
   * @param message 
   */
  public static void add(Object message){
    getLog().add(message.toString(), 6);
  }

  /**
   * Add an error string
   * @param message 
   * @param priority int between 1 and 6: <br/>
   * 1 Trace, 2 Debug, 3 Info, 4 Warn, 5 Error, 6 Fatal
   */
  public static void add(Object message, int i) {
    getLog().add(message.toString(), i);
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

}
