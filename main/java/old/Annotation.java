package old;

public class Annotation<T> {
	public T annotation;

	public Annotation(T annotation) {
		this.annotation = annotation;
	}

	/**
	 * @return the Annotation
	 */
	public T get() {
		return annotation;
	}

	/**
	 * @param value
	 *            the Annotation to set
	 */
	public void set(T annotation) {
		this.annotation = annotation;
	}

	public String getString() {
		return String.valueOf(annotation);
	}
}