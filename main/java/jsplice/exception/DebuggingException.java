/**
 * 
 */
package jsplice.exception;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 *
 */
public class DebuggingException extends RuntimeException {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * 
	 */
	public DebuggingException() {
		super();
	}

	/**
	 * @param message
	 */
	public DebuggingException(String message) {
		super(message);
	}

	/**
	 * @param cause
	 */
	public DebuggingException(Throwable cause) {
		super(cause);
	}

	/**
	 * @param message
	 * @param cause
	 */
	public DebuggingException(String message, Throwable cause) {
		super(message, cause);
	}

	/**
	 * @param message
	 * @param cause
	 * @param enableSuppression
	 * @param writableStackTrace
	 */
	public DebuggingException(String message, Throwable cause,
			boolean enableSuppression, boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}

}
