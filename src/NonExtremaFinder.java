import java.util.ArrayList;


public class NonExtremaFinder {

	/** Removes all points which are not extrema from the set of points. 
	 * In other words, only leave cases where localMax - localMin > thresh.
	 * @param firstDeriv - the set of forward derivatives.
	 * @param thresh - the threshold as described above.
	 * @param xs - the x coordinates of the points.
	 * @param ys - the y coordinates of the points.
	 * @return the points which express a sharp enough extrema.
	 */
	private int numCyclesForMax;
	private int numCyclesForMin;
	private boolean negToPos; //indicates if we have switched from - to +

	int[][] removeNonExtremas(double[] firstDeriv, double thresh,
			int[] xs, int[] ys) {

		ArrayList <Integer> indexes = new ArrayList <Integer> ();
		int maxPrev = -1; //the previous point of a local max 
		numCyclesForMax = 0; 
		numCyclesForMin = 0;
		negToPos = false;  //have we seen switch from - to +
		int maxNext = -1;

		while (numCyclesForMax <= 1 || !negToPos) {
			negToPos = false; //have we seen switch from - to +
			maxPrev = getPrevPt(firstDeriv.length, maxPrev, maxNext, true);
			maxNext = getNextPt(firstDeriv.length, maxPrev, true);
			
			while(firstDeriv[maxNext] == 0) {//if + 0...0 -
				maxNext = getNextPt(firstDeriv.length, maxNext, true);
			}
			
			//Search for switch from + to - in the derivative
			//when found, look for next switch from - to +
			if (firstDeriv[maxPrev] > 0 && firstDeriv[maxNext] < 0){
				indexes = findLocalMin(firstDeriv, thresh, indexes, maxPrev, maxNext);	
			}
		} //end of outer while
		
		int[][] newPoints = reformatPoints(xs, ys, indexes);
		return newPoints;
	}


	/**
	 * If more than 1 consecutive 0's are found, find the middle point between them and leave
	 * it only.
	 * @param firstDeriv - the unaltered set of first derivatives
	 * @param xs the x coordinates of the points of the countour
	 * @param ys the y coordinates of the points of the countour
	 * @return List containing the altered first derivative, x and y cooridnates
	 */
	private double[] removeConsecutiveZeros(double[] firstDeriv, int[] xs,
			int[] ys) {
		
		int numIterations = 0; 
		int prevI = -1;
		int nextI = -1;
		
		while(numIterations <=1) {
			//increment nextI and prevI
			prevI+=1; 
			nextI = prevI+1;
			if (prevI == firstDeriv.length -1) {
			   nextI = 0; 
			}
			if (prevI == firstDeriv.length) {
				prevI = 0; 
				nextI = 1; 
				numIterations +=1;
			}
			//done incrementing
			
			
			if (firstDeriv[prevI] == 0 && firstDeriv[nextI]==0) {
				nextI = getLastRunning0(firstDeriv, prevI, nextI);
			}
			
			
		}
		
		return null;
	}


	public int getLastRunning0(double[] firstDeriv, int prevI, int nextI) {
		while (firstDeriv[nextI] == 0 && nextI != prevI) {
			nextI+=1;
			if (nextI == firstDeriv.length)
				nextI = 0;
		}
		nextI = nextI -1;
		if (nextI < 0 ) nextI = firstDeriv.length -1;
		
		return nextI;
	}


	/** Reformat the contour points into one 2D array.
	 * Takes the points indicated in "indexes" and put them
	 * from the xs and ys sets into one 2D array where each
	 * row contains 2 columns. Each row represents a point
	 * where newPoints[n][0] = the x coordinate of the n-th point
	 * and newPoints[n][1] = the y coordinate of the n-th point.
	 * 
	 * @param xs - set of x coordinates.
	 * @param ys - set of y coordinates.
	 * @param indexes - set of indexes to keep.
	 * @return reformatted x and y point.
	 */
	private int[][] reformatPoints(int[] xs, int[] ys, ArrayList<Integer> indexes) {
		int [][] newPoints = new int [indexes.size()][2];
		Integer index = -1;
		for (int i = 0; i < indexes.size(); i++) {
			index = indexes.get(i);
			
			newPoints[i][0] = xs[index];
			newPoints[i][1] = ys[index];
			//System.out.println ("i = " + i + " index = " + 
			//		(index + 1) + " pt : " + xs[index] + 
			//		" " + ys[index]);
		}
		return newPoints;
	}


	private ArrayList<Integer> findLocalMin(double[] firstDeriv, double thresh,
			ArrayList<Integer> indexes, int maxPrev, int maxNext) {
		int minPrev = maxNext-1;//the -1 is just for set up.
		int minNext = maxNext-1;
		numCyclesForMin = 0;
		while (numCyclesForMin <=1) {	
			minPrev =getPrevPt(firstDeriv.length, minPrev, minNext, false);
			minNext = getNextPt(firstDeriv.length, minPrev, false);
			
			while (firstDeriv[minNext] == 0) { //if - 0...0 + 
				minNext = getNextPt(firstDeriv.length, minNext, false);
				System.out.println("minNext in while: " + (minNext +1));
			}
			
			boolean localMin = firstDeriv[minPrev] < 0 && firstDeriv[minNext] >0 ;
			//boolean signifficant = isSignifficant(maxPrev, minPrev, thresh, firstDeriv);
			//take the two points and "-" them and check if diff is signifficant
			if (localMin) { //&& signifficant) {
				System.out.println("minPrev = " + (minPrev +1) + " minNext " + (minNext +1));
				negToPos = true;
				if (!indexes.contains(minPrev)) {
					indexes.add(minPrev);
				}
				break; //break out of the inner for loop
			}
			
		}//end of inner while
		return indexes;
	}
	

	/**
	 * Updates the value of the next point for either a local max or a local min.
	 * Also, updates the number of times we have cycled through the list of
	 * points. If we are cycling for a local max, updates numCyclesForMax.
	 * If we are cycling for local min, updateds numCyclesForMin.
	 * 
	 * @param numDerivs - the number of derivatives.
	 * @param maxPrev - the previous point of the previous cycle.
	 * @param max - true if we are looking for a local max. False if
	 * we are looking for a local min
	 * @return
	 */
	private int getNextPt(int numDerivs, int maxPrev, boolean max) {
		int maxNext = maxPrev + 1; //the next point of a local max 
		if (maxNext > numDerivs-1) {
			maxNext = 0;
		}
		return maxNext;
	}

	/**
	 * Choose tne new previous point. If for the last cycle we were at the end of
	 * the list of points, no we cicrcle back to point 0.
	 * 
	 * @param prev - the previous point for the previous cycle.
	 * @param next - the next point for the previous cycle.
	 * @return - the new "previous" point.
	 */
	private int getPrevPt(int numPts, int prev, int next, boolean max) {
		prev +=1;
		if (prev >= numPts) {
			prev = 0;
			if (max) {
				numCyclesForMax += 1;
			} else if (!max) {
				numCyclesForMin+=1;
			}
		}
		return prev;
	}

	/** The difference is signifficant if the difference of the previous points
	 * of the two peeks is larger than the given threshold.
	 * @param maxPrev - the previous point of the local max
	 * @param minPrev - the previous point of the local min
	 * @param thresh - the threshold which determines if the difference
	 * is significant. This must be a positive number or 0.
	 * @param firstDeriv - the list of first derivatives.
	 * @return True if the difference is signifficant and false otherwise.
	 */
	private boolean isSignifficant(int maxPrev, int minPrev, double thresh,
			double[] firstDeriv) {
		//System.out.println("threshold:" + thresh + "pts: " + maxPrev + " and " + minPrev +
		//		" diff = " + Math.abs(firstDeriv[maxPrev] - firstDeriv[minPrev]) + " signifficant? " + 
		//		(Math.abs(firstDeriv[maxPrev] - firstDeriv[minPrev]) >= thresh) ); 
		if (Math.abs(firstDeriv[maxPrev] - firstDeriv[minPrev]) >= thresh) {
		
			return true;
		}
		return false;
	}

}
