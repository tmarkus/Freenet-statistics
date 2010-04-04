/*
 * Freenet statistics analysis
 * TODO:
 * - auto upload to freenet
 * - 
 */

import gnu.trove.map.hash.THashMap;
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
		
		File aFile = new File(filename);
		final int max_time = getMaxTime(aFile);
		final long size = getDimensions(max_time);
		
		final List<BitSet> bitsets = getBitSets(aFile, max_time, max_nodes);
		final long total = bitsets.size();

		final List<ArrayList<BitSet>> clusters = Dbscan(bitsets, eps, min_nodes_for_cluster, size);
	
		//generate a graph for each cluster
		int total_clustered_nodes = 0;
		for(List<BitSet> cluster : clusters)
		{
			System.out.println("cluster count = " + cluster.size());	
			generateGraph(cluster, size, outputPath);
			total_clustered_nodes += cluster.size();
		}
	
		//store the number of nodes not present in any cluster
		writeString(new File(outputPath+"unclustered.size"), Long.toString(total-total_clustered_nodes));
		
		//store the total number of nodes
		writeString(new File(outputPath+"total.size"), Long.toString(total));
	}

	/**
	 * Generate a gnuplot input file for graph generation
	 * @param cluster
	 */
	
	public static void generateGraph(List<BitSet> cluster, long size, String outputPath)
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
		long nodes_mean_sample = 0; //In how many samples are the nodes by average?
		for(BitSet node : cluster)
		{
			nodes_mean_sample += node.cardinality();
		}
		nodes_mean_sample = nodes_mean_sample / cluster.size();
		
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
		File outputFileData = new File(outputPath + name + ".data");
		writeString(outputFileData, data);
	
		String gnuplot = "set autoscale\n" + 
		"set terminal png medium\n" +
		"set output \""+name+".png\"\n" +
		"unset log\n" +
		"unset label\n" +
		"set title \"Cluster containing "+ cluster.size() +" nodes. Mean number of samples that nodes appear in: "+nodes_mean_sample+"\"\n" +
		"set xtic auto rotate \n" +
		"set ytic auto\n" +
		"set key right bottom\n" +
		"set ylabel \"Nodes in probe request sample\"\n" +
		"set xlabel \"Time\"\n" +
		"set xdata time\n" +
	    "set timefmt \"%s\"\n" +
		"plot \""+Math.abs(cluster.hashCode())+".data\" using 1:2 with lines\n";
		
		File outputFileDataGnuPLot = new File(outputPath + Math.abs(cluster.hashCode()) + ".p");
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
	
	public static List<ArrayList<BitSet>> Dbscan(List<BitSet> bitsets, double eps, int minPTS, long size)
	{
		TIntHashSet visited = new TIntHashSet(); 
		List<ArrayList<BitSet>> clusters = new ArrayList<ArrayList<BitSet>>(); 
		
		final long max_error = Math.round(size * eps);
		int noise = 0;
		System.out.println("Max error = " + max_error);
		
		for(int i=0; i < bitsets.size(); i++)
		{
			BitSet p = bitsets.get(i);
			if (!visited.contains(p.hashCode()))
			{
				visited.add(p.hashCode());
				//bitsets.remove(i); //remove p from the set
				
				final List<BitSet> neighborhood = getNeighbours(bitsets, p, max_error, visited);
				
				//System.out.println(neighborhood.size());
				
				if (neighborhood.size() < minPTS)
				{
					// p is noise, so don't put it in any cluster
				}
				else //neighborhood is large enough, so expand it
				{
					ArrayList<BitSet> cluster = new ArrayList<BitSet>();
					expandCluster(p, neighborhood, cluster, clusters, minPTS, bitsets, max_error,visited);
					if (cluster.size() >= minPTS)
					{
						System.out.println("Adding a new cluster...");
						clusters.add(cluster);
					}
				}
			}
		}
		
		System.out.println("NOISE: " + noise);
		return clusters;
	}
	
	public static List<BitSet> getNeighbours(List<BitSet> bitsets, BitSet p, long max_error, TIntHashSet visited)
	{
		List<BitSet> neighborhood = new LinkedList<BitSet>();
		
		//System.out.println("bitset contains this many elements: " + bitsets.size());
		
		for(BitSet n : bitsets)
		{
				BitSet n_clone = n.get(0, n.size());
				n_clone.xor(p);
				//System.out.println("CARDINALITY: " + n.cardinality());
				
				if (n_clone.cardinality() < max_error) neighborhood.add(n);
		}
		
		return neighborhood;
	}
	
	public static void expandCluster(BitSet p, List<BitSet> N, ArrayList<BitSet> C, List<ArrayList<BitSet>> clusters, int minPTS, List<BitSet> bitsets, long max_error, TIntHashSet visited)
	{
		C.add(p); //add p to cluster c
		final ListIterator<BitSet> N_iterator = N.listIterator();
		
		while(N_iterator.hasNext())
		{
			final BitSet p_prime = N_iterator.next();
			if (!visited.contains(p_prime.hashCode()))
			{
				visited.add(p_prime.hashCode());
				final List<BitSet> N_prime = getNeighbours(bitsets, p_prime, max_error, visited);
				if (N_prime.size() >= minPTS)
				{
					for(BitSet N_prime_element : N_prime)
					{
						if (!visited.contains(N_prime_element.hashCode())) {
							N_iterator.add(N_prime_element);
							visited.add(N_prime_element.hashCode());
						}
					}
				}
			}				
			boolean p_prime_in_cluster = false;
			for(List<BitSet> cluster : clusters)
			{
				if (cluster.contains(p_prime)) p_prime_in_cluster = true;
			}
			if (!p_prime_in_cluster) C.add(p_prime);
		}
		
		//System.out.println("visited size: " + visited.size());
		
		/*
		expandCluster(P, N, C, eps, MinPts)
		   add P to cluster C
		   for each point P' in N 
		      if P' is not visited
		         mark P' as visited
		         N' = getNeighbors(P', eps)
		         if sizeof(N') >= MinPts
		            N = N joined with N'
		      if P' is not yet member of any cluster
		         add P' to cluster C
		*/
	}
	
	public static List<BitSet> getBitSets(File aFile, int max_time, int max_nodes)
	{
		System.out.println("Creating bitsets...");
		
		/* Load a dataset */
		List<BitSet> data = new LinkedList<BitSet>();
		Map<Long, Integer> IDtoIndex = new THashMap<Long, Integer>(100000); 
		
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

							if (IDtoIndex.keySet().size() < max_nodes || IDtoIndex.containsKey(id))
							{
								if (!IDtoIndex.containsKey(id))
								{
									//System.out.println("creating new bitset: " + getDimensions(max_time));
									
									//create a new bitset
									BitSet instance = new BitSet((int) getDimensions(max_time));
									
									data.add(instance);
									IDtoIndex.put(id, data.indexOf(instance));
								}	

								//init a certain bit, at a certain location
								data.get(IDtoIndex.get(id)).set((int)time/hourly_time_interval);
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
		
		return data;
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
