import java.util.ArrayList;


public class ShapeExpander {
	/**
	 * Performs the expansion of the shape using circles.
	 */
	
	public ShapeExpander() {}

	//changed places of xStart and yStart
	public int[][] runK(double dt, int endRadius, int numPts, int xStart, int yStart,
			double[] difference, int[][] pixels, double stepK, double constant, 
			String fs) {
		
		/*
		//uncomment to see where the center point is
		System.out.println("print statement from ShapeExpander.java\n" +
		"Printing only the starting point");
		int[][] shape = new int[2][2];
		shape [0][0] = xStart;
		shape [0][1] = yStart;
		shape [1][0] = xStart+1;
		shape [1][1] = yStart+1;
		//end of uncomment to see where the center point is
		*/
		
		
		//uncomment to see shape we found
		// image points we have found
		int[][] prevResult = new int[numPts][2]; // points of previous circle
		int[][] curResult = new int[numPts][2]; // points of current circle
		int[][] shape = new int[numPts][2]; // points of the shape we have found so far
		ArrayList<Integer> filledPosInShape = new ArrayList<Integer>();

		
		
		// get circles.
		// changes radius of circle
		for (double curRadius = 0; curRadius <= endRadius; curRadius = curRadius + dt) {
			prevResult = curResult;
			curResult = calculate(dt, curRadius, numPts, xStart, yStart, constant);
			if (curRadius != 0)
				findPtsToStop(prevResult, curResult, difference, pixels, shape,
						filledPosInShape);

			if (filledPosInShape.size() == numPts)
				break;
		}
		
		//printPts(shape);
		double fullShape = Double.parseDouble(fs);
		if (fullShape > 0) {
			RecursiveShapeExpander cont = new RecursiveShapeExpander();
			shape = cont.runK1(shape, difference, pixels, fullShape,
					xStart, yStart);
		}
		//end of uncomment to see shape we found
		
		return shape;
	}
	
	/**
	 * Helper method. Prints all points in this shape.
	 * @param shape - the shape whose points we are printing.
	 */
	private void printPts(int[][] shape) {
		System.out.println("Points of the full shape: " + shape.length);
		for (int i =0; i<shape.length; i++) {
			System.out.println("i = " + (i+1) + " pt: " + shape[i][0] + ", " + 
					shape[i][1]);	
		}
		
	}


	/**
	 * Calculate the x and y values of the equation when i starts at 0 and ends
	 * at endI put each x,y pair in the 2d array arrToReturn.
	 */
	public int[][] calculate(double dt, double curRadius, int numPts, int xstart,
			int ystart, double constant) {
		// initialize variables which are defined later in the loop
		double f = 0.0;
		double g = 0.0;
		double xCoord = 0.0;
		double yCoord = 0.0;
		double power = 0.0;
		double a = 0.1;//curRadius*Math.PI/numPts;
		double s = 2 * Math.PI/numPts;
		int arrToReturn[][] = new int[numPts][2];

		for (int i = 0; i < numPts; i++) {
			g = constant * i;//1000 * a * (s *i);//(2* Math.PI/numPts) *i ; //constant * i;
			// f = R - curK*dt;
			// i * (curRadius*Pi/numPts) * Pi/numPts
			//power = (i * a * s) - (4 * Math.pow(a, 2) * dt);
			f = curRadius*dt;
			//f = Math.pow(Math.E, power); //curRadius * dt; //instead of dt do: e ^{0.1*((2pi/n)*i) - 4*0.01*dt}
			xCoord = (f) * Math.cos(g);
			yCoord = (f) * Math.sin(g);
			//System.out.println("the power = " + power + " f = " + f + 
			//		" cos = " + Math.cos(g) + " sin  = " + Math.sin(g) );
			//System.out.println("g = " + g);
			xCoord += xstart;
			yCoord += ystart;
			int[] p = { (int) Math.rint(xCoord), (int) Math.rint(yCoord) };

			arrToReturn[i] = p;
		}
		// contains all pts of this circle
		return arrToReturn;
	}
	
	/*
	 * Given a penalty value, an image and two pts of the two previous images,
	 * finds if they differ by the set penalty value or more prev[][] - endIX2
	 * array containing the pts of the previous circle [][]cur - endIX2 array
	 * containing the pts of the current circle difference - if two points in
	 * cur and prev differ by a #>difference, we stop that point image - the
	 * original image we are working with. insertPoints[][] - endI X 2 array
	 * containing all the pts that are stopped allready
	 */
	private int[][] findPtsToStop(int[][] prev, int[][] cur,
			double[] difference, int[][] pixels, int[][] shape,
			ArrayList<Integer> filledPosInShape) {
		int[] p1 = null;
		int[] p2 = null;
		for (int i = 0; i < prev.length; i++) {// go through the # of rows
	
			p1 = prev[i];
			p2 = cur[i];
			int val = satisfiesPenalty(p1, p2, pixels, difference);
	
			if (val == -1) {
				if (p1[0] < 0)
					p1[0] = 0;
				if (p1[1] < 0)
					p1[1] = 0;
			}
			// if this pt is not stopped yet, and penalty is true add the pt. to
			// the array
			if (!filledPosInShape.contains(i) && (val == 1 || val == -1)) {// stop
																			// shrinking
																			// this
																			// pt.
				// keep track of the pixels we have stopped from movement
				filledPosInShape.add(i);
				shape[i][0] = p1[0];
				shape[i][1] = p1[1];
			}
		}
		return shape;
	}
	
	/*
	 * If the difference in the vals of the 2 pts is larger then the threshold,
	 * return true p1 - a 2 cell array such that p1[0] = xCoord of p1 and p1 [1]
	 * = y coord of p1 p2 - a 2 cell array such that p2[0] = xCoord of p1 and p2
	 * [1] = y coord of p2 image - the original image difference - a 3 cell
	 * array such that difference[0] = the difference in the reds difference[1]
	 * = the difference in the blues difference[2] = the difference in the
	 * greens
	 * 
	 * returns 1 if the point must be stopped
	 */
	 int satisfiesPenalty(int[] p1, int[] p2, int[][] pixels,
			double[] difference) {
//		 System.out.println("pixels[0].length" + pixels[0].length); //x //c
//		 System.out.println("pixels.length" + pixels.length); //y //r

		 int c = pixels[0].length-1;
		 int r = pixels.length-1;
		// if any coordinate is negative, stop the point
		if (p1[0] < 0 || p1[1] < 0 || p2[0] < 0 || p2[1] < 0
				|| p2[0] >= pixels[0].length || p2[1] >= pixels.length
				|| p1[0] >= pixels[0].length || p1[1] >= pixels.length) {
			return -1;
		}
                                           //r     //c
		int p1Red = ImageData.getRed(pixels[r-p1[1]][c-p1[0]]);
		int p1Blue = ImageData.getBlue(pixels[r-p1[1]][c-p1[0]]);
		int p1Green = ImageData.getGreen(pixels[r-p1[1]][c-p1[0]]);

		// System.out.println(p2[0] + ", " + p2[1]);

		int p2Red = ImageData.getRed(pixels[r-p2[1]][c-p2[0]]);
		int p2Blue = ImageData.getBlue(pixels[r-p2[1]][c-p2[0]]);
		int p2Green = ImageData.getGreen(pixels[r-p2[1]][c-p2[0]]);

		boolean reds = (Math.abs(p1Red - p2Red) >= difference[0]);
		boolean blues = (Math.abs(p1Blue - p2Blue) >= difference[1]);
		boolean greens = (Math.abs(p1Green - p2Green) >= difference[2]);

		//System.out.println("colors 1: " + p1Red + " " + p1Blue + " " + p1Green);
		//System.out.println("colors 2: " + p2Red + " " + p2Blue + " " + p2Green );
		if (reds && blues && greens) {
			return 1;
		}
		return 0;
	}

}
