/**
 * 
 */
package jsplice.tools;

import jsplice.data.Sequence;

/**
 * @author Tobias Gresser (gresserT@gmail.com)
 * Results containing information, position and sequence
 */
public class Result {

	public double[] information;
	public int position;
	public Sequence sequence;
	/**
	 * 
	 */
	public Result(double[] informationP, int positionP, Sequence sequenceP) {
		this.information = informationP;
		this.position = positionP;
		this.sequence = sequenceP;
	}

	public double getTotalInformation() {
		return Functions.sum(information);
	}
}
