import java.util.ArrayList;


public class Stopwatch {

	private long startTime;
    private long stopTime;

    public static final double NANOS_PER_SEC = 1000000000.0;

	/**
	 start the stop watch.
	*/
	public void start()
	{	System.gc();
		startTime = System.nanoTime();
	}

	/**
	 stop the stop watch.
	*/
	public void stop()
	{	stopTime = System.nanoTime();	}

	/**
	elapsed time in secods.
	@return the time recorded on the stopwatch in seconds
	*/
	public double time()
	{	return (stopTime - startTime) / NANOS_PER_SEC;	}

	public String toString()
	{   return "elapsed time: " + time() + " seconds.";
	}

	/**
	elapsed time in nanoseconds.
	@return the time recorded on the stopwatch in nanoseconds
	*/
	public long timeInNanoseconds()
	{	return (stopTime - startTime);	}

	/**
	 * Prints information for the image. For example : 
	 * 
	 * img in apply for loop is: 
	 * C:\Users\Nona\workspace\Tracking\Original\Mouses_1-1_size_x_2.png
	 * Type of immage: 5
	 * Center of mass: (103.28125, 87.203125)
	 * Type of immage: 5
	 * Time for images: elapsed time: 0.292680465 seconds.
	 * ***********************************************
	 * 
	 * @param params
	 * @param filteredFilesInFolder
	 * @param j
	 * @param imageEffect TODO
	 */
	public void printImageInformation(ArrayList<ObjectParams> params, ArrayList<String> filteredFilesInFolder, int j, BoundaryFinder imageEffect) {
		System.out.println("img in apply for loop is: \n" + imageEffect.data.path + "\\"
				+ filteredFilesInFolder.get(j));
		
		for (int k = 0; k < params.size(); k++) 
		System.out.println("Center of mass: (" + params.get(k).xcoordinate + ", "
				+ params.get(k).ycoordinate + ")");
	
		
		System.out.println("Time for images: " + toString());
		System.out.println("\n"
				+ "***********************************************");
	}
}
