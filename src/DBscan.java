/*
 * Freenet statistics analysis
 * TODO:
 * - auto upload to freenet
 * - 
 */

import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.optimization.fitting.PolynomialFitter;
import org.apache.commons.math.optimization.general.AbstractLeastSquaresOptimizer;
import org.apache.commons.math.optimization.general.LevenbergMarquardtOptimizer;

public class DBscan {

	private static final int extrapolate_time = 0;
	private static final long cutoff_time = 0;
	private static final int max_nodes = 1000000;
	private static final int min_nodes_for_cluster = 100;
	private static final float eps = 0.1f;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {

		final String filename = args[0];
		final String outputPath = args[1];
		
		final File aFile = new File(filename);
		final int max_time = getMaxTime(aFile);
		final long bitset_size = getDimensions(max_time);
		
		final List<BitSet> bitsets = getBitSets(aFile, max_time, max_nodes);
		final long total = bitsets.size();

		final List<IDistanceMeasure> measures = new LinkedList<IDistanceMeasure>();
		measures.add(new XOR());
		measures.add(new Jaccard());
		
		
		for(IDistanceMeasure measure : measures)
		{
			LinkedList<BitSet> compare_bitsets = new LinkedList<BitSet>(bitsets);
			final List<ArrayList<BitSet>> clusters = Dbscan(measure, bitsets, compare_bitsets, eps, min_nodes_for_cluster, bitset_size);
			
			//generate a graph for each cluster
			int total_clustered_nodes = 0;
			for(List<BitSet> cluster : clusters)
			{
				System.out.println("cluster count = " + cluster.size());	
				generateGraph(measure, cluster, bitset_size, outputPath, "Cluster");
				total_clustered_nodes += cluster.size();
			}
		
			//generate a graph for the unclustered nodes
			generateGraph(measure, compare_bitsets, bitset_size, outputPath, "Graph of nodes which have not been clustered");
			
			//store the number of nodes not present in any cluster
			writeString(new File(outputPath+measure.getName()+"_unclustered.size"), Long.toString(total-total_clustered_nodes));
			
			//store the total number of nodes
			writeString(new File(outputPath+measure.getName()+"_total.size"), Long.toString(total));
		}
	}

	/**
	 * Generate a gnuplot input file for graph generation
	 * @param cluster
	 */
	
	public static void generateGraph(IDistanceMeasure measure, List<BitSet> cluster, long size, String outputPath, String title)
	{
		final long start_seconds= 1258884000; //start time of evanbd's dataset
		AbstractLeastSquaresOptimizer optimizer = new LevenbergMarquardtOptimizer();
		optimizer.setMaxIterations(Integer.MAX_VALUE);
		
		PolynomialFitter fitter = new PolynomialFitter(2, optimizer);
		
		//calculate the total number of nodes at each time index
		String data = "";
		for(int i=0; i < size; i++)
		{
			int count = 0;
			for(BitSet node : cluster)
			{
				if (node.get((int)i)) count++;
			}

			if (count >= 0)
			{
				if (i > (size-10)) fitter.addObservedPoint(1.0, i, count);
				data += (start_seconds + i*5*60*60) + "\t " + count + "\n";
			}
		}
		
		//calculate mean samples
		List<Integer> samples = new LinkedList<Integer>();
		long nodes_mean_sample = 0; //In how many samples are the nodes by average?
		for(BitSet node : cluster)
		{
			nodes_mean_sample += node.cardinality();
			samples.add(node.cardinality());
		}
		nodes_mean_sample = nodes_mean_sample / cluster.size();
		
		Collections.sort(samples);
		int minimum = samples.get(0);
		int maximum = samples.get(samples.size()-1);
		int median = samples.get(Math.round(samples.size() / 2));
		
		
		
		
		PolynomialFunction function = null;
		
		try {
			function = fitter.fit();
		} catch (MathException e1) {
			e1.printStackTrace();
		}
		
		for(long i=size; i < size+extrapolate_time; i++)
		{
			//System.out.println(i);
			double value = i;
				data += i + "\t " + function.value(value) + "\n";
		}
		
		
		String name = Integer.toString(Math.abs(cluster.hashCode()));
		
		//write data to file
		File outputFileData = new File(outputPath + measure.getName() + "_" + name + ".data");
		writeString(outputFileData, data);
	
		String gnuplot = "set autoscale\n" + 
		"set terminal png medium\n" +
		"set output \"" + measure.getName() + "_" + name+".png\"\n" +
		"unset log\n" +
		"unset label\n" +
		"set title \""+title+" containing "+ cluster.size() +" nodes.\\nMean number of samples that nodes appear in: "+nodes_mean_sample + "\\n" +
			"minimum: "+minimum+", maximum: "+maximum+", median: "+median+"\"\n" +
		"set xtic auto rotate \n" +
		"set ytic auto\n" +
		"set key right bottom\n" +
		"set ylabel \"Nodes in probe request sample\"\n" +
		"set xlabel \"Time\"\n" +
		"set xdata time\n" +
	    "set timefmt \"%s\"\n" +
		"plot \""+ measure.getName() + "_" + Math.abs(cluster.hashCode())+".data\" using 1:2 with lines\n";
		
		File outputFileDataGnuPLot = new File(outputPath + measure.getName() + "_" + Math.abs(cluster.hashCode()) + ".p");
		writeString(outputFileDataGnuPLot, gnuplot);
	}

	private static void writeString(File outputFileData, String data)
	{
		Writer output = null;
		try {
			output = new BufferedWriter(new FileWriter(outputFileData));
			output.write(data);
		} catch (IOException e) {
			e.printStackTrace();
		}
		finally {
			try {
				output.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	public static List<ArrayList<BitSet>> Dbscan(IDistanceMeasure measure, List<BitSet> bitsets, List<BitSet> compare_bitsets, double eps, int minPTS, long size)
	{
		TIntHashSet visited = new TIntHashSet(); 
		List<ArrayList<BitSet>> clusters = new ArrayList<ArrayList<BitSet>>(); 
		
		
		final double max_error = measure.getMaxError(size, eps); 
		System.out.println("Max error = " + max_error);
		
		for (BitSet p : bitsets)
		{
			if (!visited.contains(p.hashCode()))
			{
				visited.add(p.hashCode());
				List<BitSet> neighborhood = getNeighbours(measure, compare_bitsets, p, max_error, visited);
				
				if (neighborhood.size() >= minPTS)
				{
					ArrayList<BitSet> cluster = new ArrayList<BitSet>();
					expandCluster(measure, p, neighborhood, cluster, clusters, minPTS, compare_bitsets, max_error, visited);
					clusters.add(cluster);
					System.out.println("Found new cluster of size: " + cluster.size());
				}
			}
		}
		
		return clusters;
	}
	
	public static List<BitSet> getNeighbours(IDistanceMeasure measure, List<BitSet> bitsets, BitSet p, double max_error, TIntHashSet visited)
	{
		//System.out.println("Runnng getNeighbours with bitsets size: " + bitsets.size());
		List<BitSet> neighborhood = new LinkedList<BitSet>();
		
		for(BitSet n : bitsets)
		{
			if (!visited.contains(n.hashCode()))
			{
				if (measure.acceptable(measure.getDistance(n, p), max_error)) neighborhood.add(n);
			}
		}
		
		//System.out.println("done (got: " + neighborhood.size() + ")");
		
		return neighborhood;
	}
	
	public static void expandCluster(IDistanceMeasure measure, BitSet p, List<BitSet> N, ArrayList<BitSet> C, List<ArrayList<BitSet>> clusters, int minPTS, List<BitSet> bitsets, double max_error, TIntHashSet visited)
	{
		C.add(p); //add p to cluster c
		final ListIterator<BitSet> N_iterator = N.listIterator();
		HashSet<BitSet> local_visited = new HashSet<BitSet>();
		
		local_visited.addAll(N);
		
		int last_size = 0;
		
		while(N_iterator.hasNext())
		{
			final BitSet p_prime = N_iterator.next();
			if (!visited.contains(p_prime.hashCode()))
			{
				visited.add(p_prime.hashCode());
				final List<BitSet> N_prime = getNeighbours(measure, bitsets, p_prime, max_error, visited);
				if (N_prime.size() >= minPTS)
				{
					for(BitSet N_prime_element : N_prime)
					{
						if (!local_visited.contains(N_prime_element)) {
							N_iterator.add(N_prime_element);
						}
					}
				
					for(BitSet N_prime_element : N_prime)
					{
						local_visited.add(N_prime_element);
					}
					
				}
			}				
		
			if (N.size() > last_size + 100) {
				System.out.println("Growing cluster... current size: " + N.size());
				last_size = N.size();
			}
		}	
		
		for(BitSet p_prime : N)
		{
			C.add(p_prime);
			bitsets.remove(p_prime);
		}
	}
	
	public static List<BitSet> getBitSets(File aFile, int max_time, int max_nodes)
	{
		System.out.println("Creating bitsets...");
		
		/* Load a dataset */
		Map<Long, BitSet> IDtoInstance = new HashMap<Long, BitSet>(100000); 
		
		final int hourly_time_interval = 5;

		//System.out.println("Creating instances...");

		//walk through the file and set the instance values
		try {
			BufferedReader input =  new BufferedReader(new FileReader(aFile));
			try {
				String line = null; //not declared within while loop
				while (( line = input.readLine()) != null){
					final String[] splitted = line.split("\\s+");

					if (splitted.length == 2)
					{
						final long time = new Long(splitted[0]);

						if (time > cutoff_time)
						{
							final long id = new Long(splitted[1]);

							if (IDtoInstance.keySet().size() < max_nodes || IDtoInstance.containsKey(id))
							{
								if (!IDtoInstance.containsKey(id))
								{
									//System.out.println("creating new bitset: " + getDimensions(max_time));
									
									//create a new bitset
									BitSet instance = new BitSet((int) getDimensions(max_time));
									IDtoInstance.put(id, instance);
								}	

								//init a certain bit, at a certain location
								IDtoInstance.get(id).set((int)time/hourly_time_interval);
							}
						}
					}
					else
					{
						System.out.println("Bad unsplittable line: " + line);
					}
				}
			}
			finally {
				input.close();
			}
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		System.out.println("Finished creating bitsets!");
		
		return new LinkedList(IDtoInstance.values());
	}

	public static int getMaxTime(File aFile)
	{
		//determine maximum time index
		int max_time = 0;

		try {
			final BufferedReader input =  new BufferedReader(new FileReader(aFile),300000);
			try {
				String line = null; //not declared within while loop
				int current_time = 0;
				while (( line = input.readLine()) != null){
					final String[] splitted = line.split("\\s+");
					current_time = new Integer(splitted[0]);
					if (current_time > max_time) {
						System.out.println(current_time);
						max_time = current_time;
					}
				}
			}
			finally {
				input.close();
			}
		}
		catch (IOException ex){
			ex.printStackTrace();
		}

		return max_time;
	}

	private static long getDimensions(int max_time)
	{
		return ((max_time - cutoff_time) / 5); 
	}

	
}
