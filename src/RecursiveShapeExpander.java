import java.awt.Polygon;
import java.util.ArrayList;


public class RecursiveShapeExpander extends ShapeExpander {
	/**
	 * Performs the expansion of the shape using recursion.
	 */
	public RecursiveShapeExpander() {}


	/**
	 * Finds shape by expanding recursively
	 * @param initialShape - the initial shape we are given.
	 * @param difference - the differences between pixel values we must get.
	 * difference is a array of length 3.
	 * @param pixels - the original image
	 * @param threshold - the distance we want to achieve between each set of points.
	 * @return
	 */
	public int[][] runK1(int [][] initialShape, double[] difference, int[][] pixels,
			double threshold, int xStart, int yStart) { 
		ArrayList<Integer[]> prevShape = new ArrayList<Integer[]> (); //points of the shape we have found so far
		ArrayList<Integer[]> nextShape = new ArrayList<Integer[]> ();

		//put origShape into prevShape
		for (int i = 0; i < initialShape.length; i++) {
			Integer []sh = {initialShape[i][0], initialShape[i][1]};
			nextShape.add(sh);
		}

		ArrayList<Integer> largeDistanceIndexes = getIndexesOfLargeDistances(nextShape, threshold);
		while (largeDistanceIndexes.size() > 0) {
			//move prevShape in nextShape. Also add point and move it to the end. 
			prevShape = moveNextToPrev(nextShape);
			nextShape = addAndMovePts(prevShape, largeDistanceIndexes, difference,
					pixels, threshold, xStart, yStart);

			largeDistanceIndexes = getIndexesOfLargeDistances(nextShape, threshold);
			if (nextShape.size() > 1000) break;
			//System.out.println ("shape size:   " + nextShape.size());
		}
		//add the starting point: 
		
		Integer[] initialPoint = {xStart, yStart};
		nextShape.add(initialPoint);
		initialShape = reformatShape(nextShape);
		return initialShape;
	}


	private ArrayList<Integer> getIndexesOfLargeDistances(ArrayList<Integer[]> origShape, 
			double threshold) {
		//holds the index in prevShape at which a large distance occured.
		ArrayList<Integer> largeDistanceIndexes = new ArrayList <Integer>();
		for (int i = 0; i < origShape.size(); i++) {
			int x1 = origShape.get(i)[0]; 
			int y1 = origShape.get(i)[1];

			int next = i+1; 
			if (i == origShape.size()-1) {
				next = 0;
			}
			int x2 = origShape.get(next)[0]; 
			int y2 = origShape.get(next)[1];
			double xdist = Math.pow((x2-x1), 2);
			double ydist = Math.pow((y2-y1), 2);
			double sum = xdist + ydist; 
			double distance = Math.sqrt(sum);
			//System.out.println ("0. prev: " + x1 + " " + y1 + " next: " + x2 + " " + y2 + 
			//		"distance: " + distance + " threshold: " + threshold);
			if ((distance-1) > threshold) {
				largeDistanceIndexes.add(i);
			}
		}
		return largeDistanceIndexes;
	}


	private ArrayList<Integer[]> moveNextToPrev(ArrayList<Integer[]> nextShape) {
		ArrayList<Integer[]> prevShape = new ArrayList<Integer[]>();
		for (int i = 0 ; i< nextShape.size(); i++) {
			Integer[] pt = {nextShape.get(i)[0], nextShape.get(i)[1]};
			prevShape.add(pt);
		}
		return prevShape;
	}

	/**
	 * 
	 * @param shape - points of the shape where each element of type Integer[]
	 * is a point.
	 * @return the points reformatted into a 2D array, where each row is a point. 
	 */
	private int[][] reformatShape(ArrayList<Integer[]> shape) {
		int [][] newShape = new int [shape.size()][2];
		for (int i = 0 ; i < shape.size();i++) {
			newShape [i] = new int []{shape.get(i)[0], shape.get(i)[1]} ;
		}
		return newShape;
	}

	private boolean pointInImage(int[] point, int[][] pixels) {
		return (point[0] < pixels[0].length-1 &&
				point[1] < pixels.length-1 && point[0] >= 1 && point[1] >= 1);
	}

	/** @param x coordinate
	 * @param b y intercept
	 * @param slope
	 * @return the y coordinate corresponding to the given x 
	 * cooridnate.
	 */
	private int findY(int x, double b, double slope) {
		return  (int) (slope*x +b);
	}

	/** @param y coordinate
	 * @param b y intercept
	 * @param slope
	 * @return the x coordinate corresponding to the given y
	 * cooridnate.
	 */
	private int findX(int y, double b, double slope) {
		return  (int) ((y-b)/slope);
	}

	/**
	 * Adds necessary points and moves them forward.
	 * @param shape - the points of the previous shape. Format {x, y}
	 * @param largeDistanceIndexes - indexes of large distances in increasing order.
	 * @param threshold - the distance at which we want to position all points.
	 * @return the new whape after the added points are moved forward
	 */
	private ArrayList<Integer[]> addAndMovePts(ArrayList<Integer[]> shape,
			ArrayList<Integer> largeDistanceIndexes, double[] difference,
			int[][] pixels, double threshold, int xStart, int yStart) {
		int prevIndex = -1; 
		int nextIndex = -1;
		int maxIndex = shape.size()-1;
		//go over indexes in largeDistanceIndexes and add point on the line.
		for (int i = largeDistanceIndexes.size()-1; i >= 0; i--) {
			//get indexes starting with the highest. this way, the lower ones don't change
			prevIndex = largeDistanceIndexes.get(i);
			nextIndex = prevIndex +1;
			if (prevIndex == maxIndex) nextIndex = 0;
			ArrayList <Integer[]> newPts = addPoint(shape.get(prevIndex), shape.get(nextIndex),
					shape, difference, pixels, threshold, xStart, yStart);

			//System.out.println("START points added from method  addAndMovePts");
			//for (int j = 0; j<newPts.size(); j++) {
				//System.out.println("( " + newPts.get(j)[0] + ", " + newPts.get(j)[1] + ")");
			//}
			//System.out.println("END points added from method  addAndMovePts");
			
			
			/*Integer [] newPt = addPoint(shape.get(prevIndex), shape.get(nextIndex),
					shape, difference, pixels, threshold);
			 */
			//shape.add(prevIndex+1, newPt);
			shape = addPtsToShape (shape, prevIndex, newPts);
			
			//System.out.println("ALL points now:");
			//for (int j = 0; j<shape.size(); j++) {
			//	System.out.println("( " + shape.get(j)[0] + ", " + shape.get(j)[1] + ")");
			//}
			//System.out.println("END all points.");
			
		}
		//move point along the perpendicular of the line until it stops
		return shape;
	}

	/**
	 * Adds the new points to the shape starting at index prevIndex + 1
	 * @param shape
	 * @param prevIndex
	 * @param newPts
	 * @return the shape containing the new points.
	 */
	private ArrayList<Integer[]> addPtsToShape(ArrayList<Integer[]> shape,
			int prevIndex, ArrayList<Integer[]> newPts) {

		if (prevIndex+1 < shape.size()) {
			for (int i = 0 ; i < newPts.size(); i++) {
				shape.add((prevIndex+1 + i), newPts.get(i));
			}
		} else {
			for (int i = 0 ; i < newPts.size(); i++) {
				shape.add(newPts.get(i));
			}
		}
		return shape;
	}


	/**
	 * Move the mid point along the line specified by b, slope and the
	 * direction to move in (dirMove). Find the difference among the pixels
	 * and decide where to stop the point
	 * @param b = the y intercept of the line along which we will move the point.
	 * @param slope = the slope of the line along which we will move the point
	 * @param midPt point to move. midPoint[0]= x coordinate of mid point and 
	 * midPoint[1] = y coordinate of mid point.
	 * @param dirMove = booleans indicating direction in which to move.
	 * {increaseX, decreaseX, increaseY, decreaseY}
	 * @param difference 
	 * @param pixels
	 * @return
	 */
	private Integer[] movePointOut(double b, double slope, Integer[] midPt,
			boolean[] dirMove, double[] difference, int[][] pixels,
			Integer[] prevPt, int xStart, int yStart) {
		//System.err.println("CHANGE UNMOVEDMIDPT TO {prevPt[0], prevPt[1]}");
		int[] unmovedMidPt = {xStart, yStart};//{prevPt[0], prevPt[1]};//{midPt[0], midPt[1]}; //
		int[] nextMidPt = {midPt[0], midPt[1]};

		if (dirMove[0]) { // increaseX
			// p2[0] with pixels[0].length || p2[1] and pixels.length
			while (pointInImage(nextMidPt, pixels)) { //point is in img
				nextMidPt[0] += 1;
				nextMidPt[1] = findY(nextMidPt[0], b, slope);
				if (satisfiesPenalty(unmovedMidPt, nextMidPt, pixels, difference) == 1) {
					return midPt;
				} else {
					midPt = new Integer[]{nextMidPt[0], nextMidPt[1]};
				}
			} 
			return midPt;

		} else if (dirMove[1]) { //decreaseX
			while (pointInImage(nextMidPt, pixels)) {
				nextMidPt[0] -= 1;
				nextMidPt[1] = findY(nextMidPt[0], b, slope);
				if (satisfiesPenalty(unmovedMidPt, nextMidPt, pixels, difference) == 1) {
					return midPt;
				} else {
					midPt = new Integer[]{nextMidPt[0], nextMidPt[1]};
				}
			}
			return midPt;

		} else if (dirMove[2]) { // increaseY
			while (pointInImage(nextMidPt, pixels)) {
				nextMidPt[1] +=1;
				nextMidPt[0] = findX(nextMidPt[1], b, slope);
				if (satisfiesPenalty(unmovedMidPt, nextMidPt, pixels, difference) == 1) {
					return midPt;
				} else {
					midPt = new Integer[]{nextMidPt[0], nextMidPt[1]};
				}
			}
			return midPt;

		} else if (dirMove[3]){ // decreaseY
			while (pointInImage(nextMidPt, pixels)) {
				nextMidPt[1] -=1;
				nextMidPt[0] = findX(nextMidPt[1], b, slope);
				if (satisfiesPenalty(unmovedMidPt, nextMidPt, pixels, difference) == 1) {
					return midPt;
				} else {
					midPt = new Integer[]{nextMidPt[0], nextMidPt[1]};
				}
			}
			return midPt;
		}
		
		System.err.println("from movePointOut - no direction for point move was found" + 
				" so point will stay in the same place = " + midPt[0]+ " " + midPt[1] );
		//return new Integer[]{-1, -1};
		return midPt;


	}

	/**
	 * Move the mid point along the line specified by b, slope and the
	 * direction to move in (dirMove). Find the difference among the pixels
	 * and decide where to stop the point
	 * @param b = the y intercept of the line along which we will move the point.
	 * @param slope = the slope of the line along which we will move the point
	 * @param midPoint point to move. midPoint[0]= x coordinate of mid point and 
	 * midPoint[1] = y coordinate of mid point.
	 * @param dirMove = booleans indicating direction in which to move.
	 * {increaseX, decreaseX, increaseY, decreaseY}
	 * @param difference 
	 * @param pixels
	 * @return
	 */
	private Integer[] moveVertPointOut(Integer[] midPoint, boolean[] dirMove, 
			double[] difference, int[][] pixels, Integer[] prevPoint, 
			int xStart, int yStart) {
		//System.err.println("CHANGE UNMOVEDMIDPT TO {prevPt[0], prevPt[1]}");
		int[] unmovedMidPt = {xStart, yStart};//{prevPoint[0], prevPoint[1]};//{midPt[0], midPt[1]}; //
		//int[] unmovedMidPt = {prevPoint[0], prevPoint[1]};//{midPoint[0], midPoint[1]};
		int[] nextMidPoint = {midPoint[0], midPoint[1]};
	
		if (dirMove[2]) { // increaseY 
			while (pointInImage(nextMidPoint, pixels)) {
				nextMidPoint[1] += 1;
				if (satisfiesPenalty(unmovedMidPt, nextMidPoint, pixels, difference) == 1) {
					return midPoint;
				} else {
					midPoint = new Integer[]{nextMidPoint[0], nextMidPoint[1]};
				}
			}
			return midPoint;
	
		} else if (dirMove[3]){ // decreaseY
			while (pointInImage(nextMidPoint, pixels)) {
				/*System.out.println("START From Y DEcrease for vertical line");
				System.out.println("nextMidPoint: " + nextMidPoint[0] +
						", " + nextMidPoint[1]);
				System.out.println("midPoint: " + midPoint[0] + " " +  midPoint[1]);
				System.out.println("END From Y DEcreas for vertical line");
				 */
				nextMidPoint[1] -=1;
				if (satisfiesPenalty(unmovedMidPt, nextMidPoint, pixels, difference) == 1) {
					return midPoint;
				} else {
					midPoint = new Integer[]{nextMidPoint[0], nextMidPoint[1]};
				}
			}
			return midPoint;
		}
		
		//System.out.println("from moveVertPointOut - no direction for point move was found" + 
		//" so point will stay in the same place.");
		return midPoint;
		//return new Integer[]{-1, -1};
	
	}


	private Integer[] findPointOnLine(Integer [] prevPt, Integer[] nextPt) {

		double slope = prevPt[1] - nextPt[1]/prevPt[0] - nextPt[0];
		//y = mx + b ==> b = y - mx
		double b = prevPt[1] - slope*prevPt[0];
		Integer x = Math.min(prevPt[0], nextPt[0]) + 3;
		double y = slope*x + b;

		Integer[] pt = {x, (int) y};
		return pt;
	}


	/**
	 * Adds a new point bisecting the two given points. Moves the new
	 * point pependicular to the line formed by the given two points. 
	 * Stops moving the new point when it hits the border of the object
	 * or the border of the image.
	 * 
	 * @param prevPt - the x and y coordinates of the first point.
	 * @param nextPt - the x and y coordinates of the second point.
	 * @param shape - each index into the ArrayList is a point. Each point 
	 * is connected to the 2 points whose indexes are on the left and right
	 * of the index of this point.
	 * @param difference - the threshold for red, blue and green which stop
	 * the points from moving.
	 * @param pixels - the image pixels.
	 * @param threshold - the threshold at which we want to position points.
	 * @return the new point which we add between prevPt and nextPt. The point
	 * we return should be moved to a border of the object.
	 */
	private ArrayList <Integer[]> addPoint(Integer[] prevPt, Integer[] nextPt,
			ArrayList<Integer[]> shape, double[] difference,
			int[][] pixels, double threshold, int xStart, int yStart) {

		//add all points in between
		ArrayList <Integer[]> allPtsInBetween = getBetweenPoint(prevPt.clone(), nextPt, threshold);

		//System.out.println ("Num points added after adding pts(from addPoint): " + 
		//		allPtsInBetween.size() + " prev: (" 
		//		+ prevPt[0] + ", " + prevPt[1] + 
		//		") next pt: (" + nextPt[0] + ", " + nextPt[1] + ")" );
		
		//for (int k = 0; k< allPtsInBetween.size(); k++) {
		//	System.out.println("point added: ("+ allPtsInBetween.get(k)[0] + 
		//			", " + allPtsInBetween.get(k)[1] + ") ");
		//}
		//Integer[] midPoint = getBetweenPoint(prevPt, nextPt, threshold);
		ArrayList <Integer[]> movedPts = new ArrayList <Integer[]>();

		//move mid point
		//find the perpendicular
		double slope, bPerp = -1;
		double negRecip = -1;
		boolean[] xyIncrease = new boolean[5];
		boolean[] largexyIncrease = new boolean[5];
		Integer [] endPoint;
		
		int [] xyVals = new int [5];
		
		for (int i = 0; i< allPtsInBetween.size(); i++) {
			Integer[] curPoint = allPtsInBetween.get(i);	
			if (getSlope(prevPt, nextPt)[1] != 0 && getSlope(prevPt, nextPt)[0] != 0) {
				//System.out.println(" 1 ");
				slope = getSlope(prevPt, nextPt)[0];
				negRecip = -1/slope;
				bPerp = curPoint[1] - negRecip * curPoint[0];
				largexyIncrease = findSideToIncrease(shape, bPerp, negRecip, curPoint, prevPt, nextPt, pixels, difference);
			} else if ((prevPt[0] - nextPt[0]) == 0 && prevPt[1] - nextPt[1] != 0 ){ //original line is vertical so line we build is horiozntal.
				//System.out.println(" 2 ");
				bPerp = curPoint[1]; // mx+bPerp = y -> 0*m+bPerp = y -> bPerp = y
				negRecip = 0; //slope of line we build is 0
				largexyIncrease = findSideToIncrease(shape, bPerp, negRecip, curPoint,
						prevPt, nextPt, pixels, difference);
				largexyIncrease[2] = false;
				largexyIncrease[3] = false;
				//System.out.println("move: x+ or x- " + xyIncrease[0] + " and " + xyIncrease[1] );
			} else if (getSlope(prevPt, nextPt)[0] == 0 ){ //original line is horizontal so line we build is vertical
				//System.out.println(" 3 ");
				largexyIncrease = findVerticalSideToIncrease(shape, curPoint, prevPt, nextPt, pixels, difference);	
			}
			
			//xyVals
			//System.out.print("largexyIncrease = ");
			for (int j = 0; j<largexyIncrease.length; j++) {
				//System.out.print(largexyIncrease[j] + " ");
				if(largexyIncrease[j]) { 
					xyVals[j]+=1;
					largexyIncrease[j] = false;
				}
			}
			//System.out.println();
			
		}
		
		//find the max index
		int maxIndex = -1;
		int max = Integer.MIN_VALUE;
		//System.out.println("xyVals value: ");
		for(int k = 0; k < xyVals.length-1; k++) {
			//System.out.print(xyVals[k] + " ");
			largexyIncrease[k] = false;
		      if(xyVals[k] > max) {
		         max = xyVals[k];
		         maxIndex = k;
		      }
		}
		//System.out.print(xyVals[4] + " " + " max index = " + maxIndex);
		//System.out.println("\n");
		
		largexyIncrease [maxIndex] = true;
		//System.out.println("comparison: " + (allPtsInBetween.size()==1) + 
		//		" and " + (xyVals[4] ==1));
		if (xyVals[4] >= allPtsInBetween.size()/2 || 
				(allPtsInBetween.size()== 1 && xyVals[4] == 1))
			largexyIncrease[4] = true;
			
		if (allPtsInBetween.size()==1 && xyVals[4] == 0) {
			largexyIncrease[4] = false;
		}
		
		//System.out.println("largexyIncrease final: " + 
		//		largexyIncrease[0] + " " +
		//		largexyIncrease[1] + " " +
		//		largexyIncrease[2] + " " +
		//		largexyIncrease[3] + " " +
		//		largexyIncrease[4] + " " );
		
		for (int i = 0; i< allPtsInBetween.size(); i++) {
			Integer[] curPoint = allPtsInBetween.get(i);	
			
			if (getSlope(prevPt, nextPt)[1] != 0 && getSlope(prevPt, nextPt)[0] != 0) {
				//System.out.println ("!!!!!REGULAR line!");
				slope = getSlope(prevPt, nextPt)[0];
				negRecip = -1/slope;
				//line equation is y = slope*x + b
				bPerp = curPoint[1] - negRecip * curPoint[0];
				//perpendicular equation = y = negRecip*x +bPrep
				xyIncrease = findSideToIncrease(shape, bPerp, negRecip, curPoint, prevPt, nextPt, pixels, difference);
				
				if ((xyIncrease[0]== true && xyIncrease[1]== true) || 
						(xyIncrease[2]== true && xyIncrease[3] == true))
					xyIncrease = largexyIncrease;
				
				
				if (xyIncrease[4]) { //if true, we are going out
					endPoint = movePointOut(bPerp, negRecip, curPoint, xyIncrease, 
							difference, pixels, prevPt, xStart, yStart);
					movedPts.add(endPoint);
				}
				else { 
					endPoint = curPoint;
					movedPts.add(endPoint);
					//System.out.println("moving inward. BAAAD! point is: " + endPoint[0] + " " + endPoint[1]);
				}

			} else if ((prevPt[0] - nextPt[0]) == 0 ){//original line is vertical so line we build is horiozntal.
				//System.out.println("!!!!!HORIZONTAL perp line!!!! SLOPE = 0");
				bPerp = curPoint[1]; // mx+bPerp = y -> 0*m+bPerp = y -> bPerp = y
				negRecip = 0; //slope of line we build is 0
				xyIncrease = findSideToIncrease(shape, bPerp, negRecip, curPoint, prevPt, nextPt, pixels, difference);
				xyIncrease[2] = false;
				xyIncrease[3] = false;
				//answer = {!xUp, !xDn, !yUp, !yDn, true};
				
				if ((xyIncrease[0]== true && xyIncrease[1]== true) || 
						(xyIncrease[2]== true && xyIncrease[3] == true))
					xyIncrease = largexyIncrease;
				

				if (xyIncrease[4]) { //if true, we are going out
					endPoint = movePointOut(bPerp, negRecip, curPoint, xyIncrease, 
							difference, pixels, prevPt, xStart, yStart);
					movedPts.add(endPoint);
				}
				else { 
					endPoint = curPoint;
					movedPts.add(endPoint);
					//System.out.println("moving inward. BAAAD! point is: " + endPoint[0] + " " + endPoint[1]);
				}

			} else if (getSlope(prevPt, nextPt)[0] == 0 ){ //original line is horizontal so line we build is vertical
				//System.out.println("!!!!!VERTICAL perp line!!!! SLOPE IS UNDEFINED");
				// no matter what y is, x stays the same
				xyIncrease = findVerticalSideToIncrease(shape, curPoint, prevPt, nextPt, pixels, difference);
				
				if ((xyIncrease[0]== true && xyIncrease[1]== true) || 
						(xyIncrease[2]== true && xyIncrease[3] == true))
					xyIncrease = largexyIncrease;
				
				
				if (xyIncrease[4] ) { //if true, we are going out
					endPoint = moveVertPointOut(curPoint, xyIncrease, 
							difference, pixels, prevPt, xStart, yStart);
					movedPts.add(endPoint);	
				}
				else {
					endPoint = curPoint;
					movedPts.add(endPoint);
					//System.out.println("moving inward. BAAAD! point is: " + endPoint[0] + " " + endPoint[1]);
				}
				//return endPoint;	
				
			}
		}//end of for loop.
		
		return movedPts; //return endPoint;

	}


	/**
	 * Get a point that falls between prevPt and nextPt. If the distance
	 * between prevPt and nextPt is small enough - get the mid point. Otherwise
	 * get a point that is close to one of the end points. 
	 * @param prevPt - one of the end points
	 * @param nextPt - the other of the end points.
	 * @return arraylist where each element is an array with 2 elements. 
	 * The first is the x coordinate and the second is the y coordinate.
	 */
	private ArrayList<Integer[]> getBetweenPoint(Integer[] prevPt, Integer[] nextPt,
			double threshold) {
		ArrayList <Integer[]> pts = new ArrayList <Integer[]> ();

		//find the equation of the poins
		double[] lineSlope = getSlope(prevPt, nextPt);
		pts = getPoints(lineSlope, prevPt, nextPt, threshold);
		return pts;
	}


	/**
	 * Finds the starting coordinates of the points which we insert
	 * betwen prevPt and nextPt. 
	 * @param lineSlope - slope of the line
	 * @param prevPt - the first of the two points
	 * @param nextPt - the next of the two points
	 * @param threshold - the distance at which we will insert points.
	 * @return starting positions of the inserted points.
	 */
	private ArrayList<Integer[]> getPoints(double[] lineSlope,
			Integer[] prevPt, Integer[] nextPt, double threshold) {

		ArrayList <Integer[]> pts = new ArrayList <Integer[]> ();

		if (lineSlope[1] == 0.0 ) { //vertical line
			pts = addVerticalPoints(prevPt, nextPt, threshold, pts);
		}

		else { // normal line
			pts = addNonVerticalPoints(prevPt, nextPt, threshold);
		}
		return pts;
	}


	/**
	 * Finding the correct location on which to put the new points on the line:
	 * http://stackoverflow.com/questions/12550365/calculate-a-point-along-the-line-a-b-at-a-given-distance-from-a
	 * 
	 * @param prevPt = A
	 * @param nextPt = B
	 * @param threshold
	 * @return
	 */
	private ArrayList<Integer[]> addNonVerticalPoints(Integer[] prevPt,
			Integer[] nextPt, double threshold) {
		ArrayList<Integer[]> pts = new ArrayList<Integer[]>();

		//calculate the vector of prevPt and nextPt:
		Integer [] vab = new Integer[] {nextPt[0]- prevPt[0], nextPt[1] - prevPt[1]};

		//Calculate the length of vector AB
		double vabLength = Math.sqrt(Math.pow(vab[0],2) + Math.pow(vab[1],2));

		// Calculate the unit vector
		Double [] unitVec = new Double [] {vab[0]/vabLength, vab[1]/vabLength};

		// Calculate the vector with length "threshold"
		Double [] vecWithLen = new Double[] { (unitVec[0]*threshold),  (unitVec[1]*threshold)};
		
		// How many points will be inserted between A and B?
		int numPts = (int) ((vabLength -1)/threshold); 
		/*System.out.println("points are: " + prevPt[0] + ", " + prevPt[1] + " and " + 
				nextPt[0] + ", " +  nextPt[1]);
		System.out.println("vabLength = " + vabLength + "  threshold: " + threshold + 
				 " numPts: " + (vabLength-1.0)/threshold);
		System.out.println("num points to add " + numPts);
		System.out.println("from addNonVerticalPoints: " +
				" unit vec = " + unitVec[0] + " " + unitVec[1] +
				" with length: " + vabLength + 
				" num pts to insert: " + numPts + " ");
		System.out.println("vec with len:" + vecWithLen[0] + " " + vecWithLen[1]);
		*/
		for (int i = 0; i<numPts; i++) {//

			//prevPt[0] = (int) (prevPt[0] + vecWithLen[0]);
			//prevPt[1] = (int) (prevPt[1] + vecWithLen[1]);
			pts.add(new Integer []{
					(int) Math.round(prevPt[0] + (i+1)*vecWithLen[0]), 
					(int) Math.round(prevPt[1]+(i+1)*vecWithLen[1])});
		}
		return pts;
		/*if (prevPt[0] > nextPt[0] ) { //subtract from prev
			while(prevPt[0] >= nextPt[0]) {
				prevPt[0] = (int) (prevPt[0] - vecWithLen[0]);
				prevPt[1] = (int) (prevPt[1] - vecWithLen[1]);
				pts.add(new Integer []{prevPt[0], prevPt[1]});
			}
		} else if (prevPt[0] < nextPt[0]) { //add to prev 
			while (prevPt[0] <= nextPt[0] ) {
				prevPt[0] = (int) (prevPt[0] + vecWithLen[0]);
				prevPt[1] = (int) (prevPt[1] + vecWithLen[1]);
				pts.add(new Integer []{prevPt[0], prevPt[1]});
				System.out.println("adds point:" + prevPt[0] + ", " + prevPt[1]);
			}
		}*/

	}


	public ArrayList <Integer[]> addVerticalPoints(Integer[] prevPt, Integer[] nextPt,
			double threshold, ArrayList<Integer[]> pts) {
		//go from previous to next point
		if ((prevPt[1] - threshold) > nextPt[1] ) { //subtract from prev
			while(prevPt[1] > nextPt[1]) {
				prevPt[1] = (int) (prevPt[1] - threshold);
				pts.add(new Integer []{prevPt[0], prevPt[1]});	
			}
		} else { //add to prev 
			while ((prevPt[1] + threshold) < nextPt[1]) {
				prevPt[1] = (int) (prevPt[1] + threshold);
				pts.add(new Integer []{prevPt[0], prevPt[1]});	
			}
		}
		return pts;
	}


	/**
	 * Finds the slope of the line that is formed by the points prev and next.
	 * If the line is vertical, the slope is null.
	 * 
	 * @param prev - point where the 1st coordinate is the x coordinate and the 
	 * 2nd coordinate is the y coordinate.
	 * @param next - point where the 1st coordinate is the x coordinate and the 
	 * 2nd coordinate is the y coordinate.
	 * @return a tuple. The first element is the slope, the second is 0 if the slope
	 * is undefined and 1 if it is defined.
	 */
	private double[] getSlope (Integer[] prev, Integer [] next) {
		double yPart = (double)(prev[1] - next[1]);
		double xPart = (double)(prev[0] - next[0]);

		if (xPart == 0.0) { //vertical line has undefined slope.
			return new double [] {0.0, 0.0};
		}
		return new double [] {yPart/xPart, 1.0};
	}

	/**
	 * Uses the function shown here: 
	 * http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html#The Method
	 * to find a point on the line which is not inside the polygon we have so far. 
	 * @param shape - current shape
	 * @param bPerp - b value of the perpendicular line
	 * @param slope - slope of the perpendicular line
	 * @param midPoint - point to move.
	 * @return
	 */
	private boolean[] findSideToIncrease(ArrayList<Integer[]> shape,
			double bPerp, double slope, Integer[] midPoint, 
			Integer[] prevPt, Integer[] nextPt, int[][] pixels, 
			double[]difference) {

		int increaseX = midPoint[0]+1;
		int decreaseX = midPoint[0]-1;
		int increaseY = midPoint[1]+1;
		int decreaseY = midPoint[1]-1;

		int increaseXY = findY(increaseX, bPerp, slope);
		int decreaseXY = findY(decreaseX, bPerp, slope);
		int increaseYX = findX(increaseY, bPerp, slope);
		int decreaseYX = findX(decreaseY, bPerp, slope);

		int [] a = {increaseX, increaseXY};
		int [] b = {decreaseX, decreaseXY};
		int [] c = {increaseYX, increaseY};
		int [] d = {decreaseYX, decreaseY};

		boolean xUp = polygonContains(a,shape);
		boolean xDn = polygonContains(b,shape);
		boolean yUp = polygonContains(c,shape);
		boolean yDn = polygonContains(d,shape);

		boolean [] answer = {!xUp, !xDn, !yUp, !yDn, true};
		return answer;

		/*
		//if midPoint is in the shape 
		int[] prevInt = {prevPt[0], prevPt[1]};
		int[] midInt = {midPoint[0], midPoint[1]};

		 if (satisfiesPenalty(prevInt, midInt, pixels, difference) == 0) {			
			//{increaseX, decreaseX, increaseY, decreaseY, takePrev}
			boolean [] answer = {!xUp, !xDn, !yUp, !yDn, true};
			return answer;
		}

		//if midPoint is not in the shape
		else {	
			//{increaseX, decreaseX, increaseY, decreaseY, takePrev}
			boolean [] answer = {xUp, xDn, yUp, yDn, false};
			return answer;
		}
		 */

	}


	/**
	 * If the contour line is horizontal, we have to move our points 
	 * vertically. This method finds out in which vertical direction
	 * we need to move the points.
	 * In a vertical line, no matter what y is, x stays the same
	 * @param shape
	 * @param midPoint
	 * @return
	 */
	private boolean[] findVerticalSideToIncrease(ArrayList<Integer[]> shape,
			Integer[] midPoint, Integer[] prevPt, Integer[] nextPt, 
			int[][] pixels, 
			double[]difference) {

		int increaseY = midPoint[1]+1;
		int decreaseY = midPoint[1]-1;

		int [] c = {midPoint[0], increaseY};
		int [] d = {midPoint[0], decreaseY};

		boolean yUp = polygonContains(c, shape);
		boolean yDn = polygonContains(d, shape);

		boolean [] answer = {false, false, !yUp, !yDn, true};
		return answer;

		/*
		int[] prevInt = {prevPt[0], prevPt[1]};
		int[] midInt = {midPoint[0], midPoint[1]};

		if (satisfiesPenalty(prevInt, midInt, pixels, difference) == 0) {
			//{increaseX, decreaseX, increaseY, decreaseY, takePrev}
			boolean [] answer = {false, false, !yUp, !yDn, true};
			return answer;
		}
		else {
			//{increaseX, decreaseX, increaseY, decreaseY, takePrev}
			boolean [] answer = {false, false, yUp, yDn, false};
			return answer;
		}
		 */
	}



	/**
	 * Return true if the given point is contained inside the boundary.
	 * Return false otherwise.
	 * See: http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
	 * @param point The point to check
	 * @return true if the point is inside the boundary, false otherwise
	 *
	 */
	private boolean polygonContains(int[] point, ArrayList<Integer[]> points) {
		/* int i;
      int j;
      boolean result = false;
      for (i = 0, j = points.size() - 1; i < points.size(); j = i++) { 
        if ((points.get(i)[1] > point[1]) != (points.get(j)[1] > point[1]) &&
            (point[0] < (points.get(j)[0] - points.get(i)[0]) * (point[1] - points.get(i)[1]) / (points.get(j)[1]-points.get(i)[1]) + points.get(i)[0])) {
          result = !result;
         }
      }
      return result;
		 */
		int[] xs = new int [points.size()];
		int[] ys = new int [points.size()];
		for (int i = 0 ; i<points.size(); i++) {
			xs[i]= points.get(i)[0];
			ys[i]= points.get(i)[1];
		}
		Polygon pol = new Polygon(xs, ys, points.size()); 
		if (pol.contains(point[0], point[1]))
			return true;
		else return false;
	}

}
