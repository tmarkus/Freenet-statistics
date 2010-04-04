import java.util.BitSet;


public class Jaccard implements IDistanceMeasure{

	@Override
	public boolean acceptable(double distance, double maxError) {
		if (distance >= maxError) 			return true;
		else								return false;
	}

	@Override
	public double getDistance(BitSet a, BitSet b) {
		BitSet clone_a = (BitSet) a.clone();
		BitSet clone_b = (BitSet) a.clone();
	
		clone_a.and(b);
		clone_b.or(a);
		
		if (clone_b.cardinality() == 0) return 0.0;
		else {
			return new Double(clone_a.cardinality()) / clone_b.cardinality();
		}
	}

	@Override
	public double getMaxError(long size, double epsilon) {
		return 1.0 - epsilon;
	}

	@Override
	public String getName() {
		return "Jaccard";
	}

}
