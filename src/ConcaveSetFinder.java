import java.awt.image.BufferedImage;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedList;


public class ConcaveSetFinder {

	/**
	 * Finds the concave set and returns the points which make it up.
	 * @param comX center of mass x
	 * @param comY center of mass y
	 * @param xs   all x vals
	 * @param ys   all y vals
	 * @return 2D array of integers. Column 0 is x vals. Column 1 is y vals. 
	 */
	public int [][] findConcaveSet(int comX, int comY, 
			int[] xs, int[] ys) {

		// checks that distances will not be empty
		if(xs != null && xs.length >0 && ys != null && ys.length >0) {
			// the distance between each point and the COM
			double[] distances = findDistances(comX, comY, xs, ys);
			double thresh = findThreshold(distances);
			double[] firstDeriv = takeDeriv(distances);
			
			NonExtremaFinder extrema = new NonExtremaFinder();
			int [][] newPoints = extrema.removeNonExtremas(firstDeriv, thresh, xs, ys);

			if (newPoints.length >3) newPoints = keepConcavePoints(newPoints);
			
			
			System.out.println("num points:" + firstDeriv.length);
			for (int i = 0 ; i< firstDeriv.length; i++)
				System.out.println("i = " + (i+1) +" : dist = " + distances[i] + " 1st deriv= " + 
						firstDeriv[i] + "pt: (" + xs[i] + ", " + ys[i] + ")");
			return newPoints;
		}
		System.out.println("NULL RETURNED! THIS SHOULD NOT HAPPEN!!!!");
		return null;

		/*double[] secondDeriv = takeDeriv(firstDeriv);
		// get the change from - to +
		int [][] newPoints = getMins(secondDeriv, xs, ys);
		// find concave point
		System.out.println("num points:" + secondDeriv.length);
		for (int i = 0 ; i< secondDeriv.length; i++)
			System.out.println("i = " + (i+1) +" : dist = " + distances[i] + " 1st deriv= " + 
					firstDeriv[i]+ " 2nd deriv= " +  secondDeriv[i]);//+ " " + newPoints[i][1]);


		newPoints = keepConcavePoints(newPoints);
		return newPoints;
		 */

		/*
		ArrayList<Integer> localMinIndexes = findMinDistanceIndexes(distances);
		System.out.println("local min distances: " + localMinIndexes.size());
		int [][] newPoints = new int [localMinIndexes.size()][2];

		//reshape the array 
		for (int i = 0 ; i< newPoints.length; i++) {
			newPoints [i][0] = xs[localMinIndexes.get(i)];
			newPoints [i][1] = ys[localMinIndexes.get(i)];
		}

		//find concave point
		newPoints = keepConcavePoints(newPoints);

		return newPoints;
		 */

	}


	/**
	 * Given a non-empty set of distances find the threshold which we will use
	 * for the derivatives. Threshold = floor(abs(max - min)/10). This means that
	 * the threshold should be a percentage of the difference of the largest and the 
	 * smallest distance. 
	 * @param distances - A non empty array of distances. 
	 * @return the threshold.
	 */
	private double findThreshold(double[] distances) {
		double min = distances[0];
		double max = distances[0];

		for (int i = 0; i< distances.length; i++) {	
			if (distances[i] < min) {
				min = distances[i];
			}
			if (distances[i] > max) {
				max = distances[i];
			}
		}//end of for loop

		double threshold =  Math.ceil((Math.abs((max - min))/10));
		return threshold;
	}

	/**
	 * Get all positive values the list containing second derivatives. Then,
	 * get the xs and ys points associated with each such value.
	 *
	 * @param derivs the second derivatives found using forward differences.
	 * @param xs   all original x vals
	 * @param ys   all original y vals
	 * @return The new points of the concave set
	 */
	private int[][] getMins(double[] derivs, int[] xs, int[] ys) {
		ArrayList <Integer> indexes = new ArrayList <Integer>();
		for (int i = 0; i < derivs.length; i++) {
			if (derivs[i] > 0)
				indexes.add(i);
		}

		int [][] newPoints = new int [indexes.size()][2];
		for (int i = 0; i < indexes.size(); i++) {
			newPoints[i][0] = xs[i];
			newPoints[i][1] = ys[i];
		}


		return newPoints;
	}

	/**
	 * Given a set of distances, takes the forward derivative of those distances.
	 * @param distances distances from center of mass to a given point. 
	 * for example: distances[0] = distance from point 0 to center of mass.
	 * @return forward derivative of the points given in ds.
	 * The forward derivative is found as follows: d[i+1] - d[i]= deriv
	 */
	private double[] takeDeriv(double[] distances) {
		

		double[] deriv = new double[distances.length];
		for (int i = 0; i < distances.length; i++) {
			if (i < (distances.length - 1)) {
				deriv[i+1] = roundToTenth(distances[i+1] - distances[i]);
			}
			else {// last point
				deriv[0] = roundToTenth(distances[0] - distances[i]);
			}
		}
		return deriv;
	}
	
	private double roundToTenth (double value){
        //DecimalFormat df = new DecimalFormat("0.0"); //this creates a new instance of decimalformat and tell it that you want the formatted
        //System.out.println(df.format(value)); // this would print out 54.94
        //String s = df.format(value);
        //value = Double.parseDouble(s);
        value = Math.round(value);
        return  value;
  }

	/**
	 * Finds the points which form a concave set
	 * @param points all points in the minimum set. May not form a concave set.
	 * @return the set of convex points
	 */
	private int[][] keepConcavePoints(int[][] points) {

		ArrayList<Integer[]> pts = new ArrayList<Integer[]>();
		// add to linked list
		for (int i = 0; i <points.length; i++) {
			Integer [] curArr = {points[i][0], points[i][1]};
			pts.add(curArr);
		}

		int prev = -1; // previous i
		int next = -1; // next i  
		for (int cur = 0; cur < pts.size(); cur++) {
			if (cur == pts.size()-1) {
				prev = cur - 1; 
				next = 0;
			} 

			else if (cur == 0) {
				prev = pts.size()-1; 
				next = cur+1;
			}

			else if (cur >= 1) { 
				prev = cur - 1; 
				next = cur + 1; 
			}

			boolean convex = isConvex(next, cur, prev, pts);
			if (convex) {
				//remove cur
				//System.out.println("points: " + (prev +1) + ", " + (cur +1) + ", " + (next +1) + " removing: " + (next+1));
				//System.out.println("point vals: " + pts.get(next)[0] + " " + pts.get(next)[1]);
				pts.remove(next);
				cur = 0;
				//cur = cur-1;
				//if (cur < 0) cur = 0;
			}
		}


		System.out.println("remaining points: " );
		int [][] newPoints = new int[pts.size()][2];
		for (int i = 0; i <newPoints.length; i++) {
			newPoints[i][0] = pts.get(i)[0];
			newPoints[i][1] = pts.get(i)[1];
			System.out.println("i = " + (i+1) +" : val = " + newPoints[i][0] + " " + newPoints[i][1] );
		}
			
		return newPoints;
	}

	/**
	 * Finds out whether the three points are convex.
	 * @param iNext current (largest index)
	 * @param iCur previous
	 * @param iPrev previous before previous
	 * @param pts all points
	 * @return true if convex, false otherwise
	 */
	private boolean isConvex(int iNext, int iCur, int iPrev, ArrayList<Integer[]> pts) {
		double xNext = pts.get(iNext)[0];
		double yNext = pts.get(iNext)[1];
		double xCur = pts.get(iCur)[0];
		double yCur = pts.get(iCur)[1];
		double xPrev = pts.get(iPrev)[0];
		double yPrev = pts.get(iPrev)[1];

		double leftHandSide = (xCur - xPrev)*(yNext - yCur);//-30
		double rightHandSide =(xNext - xCur)*(yCur - yPrev);//-13
		//System.out.println("next pt coordinates: " + xNext + " "+ yNext);
		//System.out.println("point: (45 28?)" + xNext + " "+ yNext);
		//System.out.println("l = " + leftHandSide + " r " + rightHandSide + " <= " + (leftHandSide <= rightHandSide));
		if (leftHandSide <= rightHandSide) //dad said it has to be: <
			//I think it has to be: >=	
			return true;//false

		else return false; //true
	}

	private ArrayList<Integer> findMinDistanceIndexes(double[] distances) {

		ArrayList <Integer> localMinIndexes = new ArrayList<Integer>();
		//for each point compare it to the one before it and the one after it.
		int prevIndex = distances.length-1; 
		int nextIndex = 1;
		for (int cur = 0; cur < distances.length; cur++) {

			prevIndex = cur-1;
			nextIndex = cur+1;

			if (cur == distances.length-1) {
				nextIndex = 0; 
			}
			if (cur == 0) {
				prevIndex = distances.length-1; 
				nextIndex = 1;
			}
			//System.out.println("prev/cur/next: " + prevIndex());
			if (isLocalMin(distances[prevIndex], distances[cur], distances[nextIndex])) {
				localMinIndexes.add(cur);
			}	
		}
		return localMinIndexes;
	}

	/**
	 * For cur to be a local min, it must be smaller then both prev and next.
	 * This method finds if cur is a local min.
	 * 
	 * @param prev a distance between com and a point on the shape
	 * @param cur a distance between com and a point on the shape
	 * @param next a distance between com and a point on the shape
	 * @return true if cur is a local min and false if it is not.
	 */
	private boolean isLocalMin(double prev, double cur, double next) {
		if (prev >= cur && next >= cur)
			return true;
		return false;
	}

	/**
	 * Finds the distance between each point of the contour and
	 * the center of mass.
	 * @param comX = x value of the center of mass
	 * @param comY = y value of the center of mass
	 * @param xs = the x values of all points in the contour
	 * @param ys = the y values of all points in the contour
	 * @return a 1D array of doubles containing distances.
	 */
	private double [] findDistances(int comX, int comY, int[] xs, int[] ys) {
		double [] distances = new double [xs.length]; 
		int x ; int y ; double dSquared;

		for (int i = 0; i < distances.length; i++) {
			//use the distance formula:
			dSquared = Math.pow(xs[i] - comX, 2) + Math.pow(ys[i] - comY, 2);
			distances[i]= Math.sqrt(dSquared);
		}
		return distances;
	}


}
