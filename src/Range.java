package srma;

public class Range {

	public int referenceIndex = -1;
	// one based
	public int startPosition = -1;
	public int endPosition = -1;

	public Range(int referenceIndex, int startPosition, int endPosition)
	{
		this.referenceIndex = referenceIndex;
		this.startPosition = startPosition;
		this.endPosition = endPosition;
	}
}

