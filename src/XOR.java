import java.util.BitSet;


public class XOR implements IDistanceMeasure {

	@Override
	public double getDistance(BitSet a, BitSet b) {
		BitSet n_clone = a.get(0, a.size());
		n_clone.xor(b);
		return n_clone.cardinality();
	}

	@Override
	public double getMaxError(long size, double epsilon) {
		return Math.round(size * epsilon);
	}

	@Override
	public boolean acceptable(double distance, double maxError) {
		if (distance < maxError) 	return true;
		else 						return false;
	}

	@Override
	public String getName() {
		return "XOR";
	}

	
	
}
