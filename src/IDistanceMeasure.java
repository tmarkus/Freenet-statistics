import java.util.BitSet;

public interface IDistanceMeasure {

	public double getDistance(BitSet a, BitSet p);
	public double getMaxError(long size, double epsilon);
	public boolean acceptable(double distance, double max_error);
	public String getName();
	
}
