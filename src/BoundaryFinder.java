import java.awt.*;
import java.awt.image.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.lang.Math;
import java.io.*;
//import ObjectParams.java;


import com.sun.image.codec.jpeg.*;


/**
 * Finds the boundary of an object inside an image.
 * 
 * @author Nona Sirakova
 */
public class BoundaryFinder {
	public static double R = 0.0;
	static int endK;

	static int endI;
	static private double step;
	static double dt;
	ImageData data = new ImageData("", "");
	ConcaveSetFinder csFinder = new ConcaveSetFinder();

	/**
	 * Applys the effect to image data.
	 * 
	 * @param pixels
	 *            the source two-dimensional array of image data.
	 * 
	 * @return a two-dimensional array of image data that is the result of
	 *         applying the effect to the source image.
	 */

	/**
	 * @return a String describing the effect.
	 */
	public String getDescription() {
		String ds= "Draw Shape";
		String concaveSet = "Draw Concave Set";
		return ds; 
	}
	


	public BufferedImage apply(BufferedImage image) {
		paths();
		listParamFileOptions();
		String paramFile = "Param_Files/";
		Scanner sc = new Scanner(System.in);
		paramFile += sc.nextLine();

		ArrayList<ObjectParams> params = initializeParamsFromFile(paramFile);
		printCOM(params);//prints center of mass of each object in the frame

		//set options according to the parameter file.
		String unfilteredImage = ""; // ask user for this information
		String orName = "";  // ask user for this information
		endK = 0; // number of circles to draw in image
		dt = 1;
		String unfilteredPath = params.get(0).unfilteredPath;
		int xcoordinate = params.get(0).xcoordinate;
		int ycoordinate = params.get(0).ycoordinate;
		ImageData.red = params.get(0).red;
		ImageData.blue = params.get(0).blue;
		ImageData.green = params.get(0).green;
		endI = params.get(0).getEndI();
		dt = params.get(0).dt;
		String circles = params.get(0).circles;
		String cs = params.get(0).concaveSet; //concave set
		String nums = params.get(0).numbers;
		String fullShape = params.get(0).fullShape;
		step = dt / 100; // configure step to be dt/100
		double[] diff = { ImageData.red, ImageData.blue, ImageData.green };
		double constant = 1.155 * Math.PI * 1.73205 / endI;

		// find all files in folder
		ArrayList<String> filteredFilesInFolder = findAllFilesInFolder(data.path);
		ArrayList<String> unfilteredFilesInFolder = findAllFilesInFolder(unfilteredPath);

		for (int j = 0; j < filteredFilesInFolder.size(); j++) {

			File f = new File(data.path + "\\" + filteredFilesInFolder.get(j));
			
			// update the image
			image = data.createBufferedImage(f, null,false);
			unfilteredImage = unfilteredPath + "\\"
					+ unfilteredFilesInFolder.get(j);
			orName = unfilteredFilesInFolder.get(j);
			
			// cleave off the extension on all 3 files and attach the new ending
			// and extension
			int dot = orName.lastIndexOf(".");
			orName = orName.substring(0, dot);
			
			//int dot2 = filteredFilesInFolder.get(j).lastIndexOf(".");
			//String prefix = filteredFilesInFolder.get(j).substring(0, dot2);

			if (xcoordinate <= 0 || ycoordinate <= 0 || orName.length() == 0
					|| unfilteredImage.length() == 0) {
				illegalParams(unfilteredImage, orName, xcoordinate, ycoordinate);
			}

			// end asking for info
			Stopwatch clock = new Stopwatch();
			clock.start();

			R = calculateR(image);
			//set params
			ArrayList<Integer> csXs;
			ArrayList<Integer> csYs;
			ArrayList<ArrayList<Integer>> everybdCSxs = new ArrayList<ArrayList<Integer>>();
			ArrayList<ArrayList<Integer>> everybdCSys = new ArrayList<ArrayList<Integer>>();
			//ArrayList<int[]> everybdCSxs = new ArrayList<int[]>();
			//ArrayList<int[]> everybdCSys = new ArrayList<int[]>();
			int[][] pixels = ImageData.imageToPixels(image);
			ArrayList<ArrayList<Integer>> everybdxs = new ArrayList<ArrayList<Integer>>();
			ArrayList<ArrayList<Integer>> everybdys = new ArrayList<ArrayList<Integer>>();
			boolean drawNumbers = false;
			
			for (int k = 0; k < params.size(); k++) {
				unfilteredPath = params.get(k).unfilteredPath;
				xcoordinate = params.get(k).xcoordinate;
				ycoordinate = params.get(k).ycoordinate;
				ImageData.red = params.get(k).red;
				ImageData.blue = params.get(k).blue;
				ImageData.green = params.get(k).green;
				endI = params.get(k).getEndI();
				dt = params.get(k).dt;
				circles = params.get(k).circles;
				cs = params.get(k).concaveSet;
				nums = params.get(k).numbers;
				fullShape = params.get(k).fullShape;
				// # of circles we will make is no more then the length of the
				// diagonal
				// diagonal = (length of the diagonal)^2 and endK = diagonal
				double diagonal = pixels.length * pixels.length
						+ pixels[0].length * pixels[0].length;
				endK = (int) Math.ceil(Math.sqrt(diagonal));
				//done setting params
				
				// return a 2d array containing in col 1 the x and in col 2 the
				// y of the pt
				
				int c = pixels[0].length-1;
				int r = pixels.length-1;
				xcoordinate = c - xcoordinate;
				ycoordinate = r - ycoordinate;
				ShapeExpander se = new ShapeExpander();
				int[][] points = se.runK(dt, endK, endI, xcoordinate, ycoordinate,
						diff, pixels, step, constant, fullShape);

				// divide up the points into an array of xs and ys
				extractXsAndYs(params, k, points, c, r);
								
				if (circles.equals("y") || circles.equals("Y")) {
					data.drawCircles(this, params, xcoordinate, ycoordinate, constant, k);
				}
				
				
				if (nums.equals("y") || nums.equals("Y")) {
					drawNumbers = true; 
				}
				// put the shapes of all objects of the img in two arrayLists
				everybdxs.add(params.get(k).xs);
				everybdys.add(params.get(k).ys);

				// find the center of mass (COM). For the x coordinate of the
				// center of mass - add up all x coordinates of the image we 
				// have caught and divide the result by the number of coordinates taken.
				// for the y coordinate - add up all y coordinates and divide by 
				// the number of coordinates taken
				double[] centerOfMass = findCOM(params.get(k).xs, params.get(k).ys);

				params.get(k).xcoordinate = Math.round((float) centerOfMass[0]);
				params.get(k).ycoordinate = Math.round((float) centerOfMass[1]);
				
				//use center of mass and everybdxs. and everybdys.
				//if params contains a yes at location 9 find concave set
				//else don't find concave set
				if (cs.equals("y") || cs.equals("Y")) {
					int [] xArr = ImageData.arrListToArray(params.get(k).xs);
					int [] yArr = ImageData.arrListToArray(params.get(k).ys);
					int [][]convexSetPoints = csFinder.findConcaveSet(params.get(k).xcoordinate, 
							params.get(k).ycoordinate, xArr, yArr);
					csXs = new ArrayList<Integer>();
					csYs = new ArrayList<Integer>();
					
					for (int curIndex = 0 ; curIndex < convexSetPoints.length; curIndex++ ) {
						csXs.add( convexSetPoints[curIndex][0]);
						csYs.add(convexSetPoints[curIndex][1]);
						
					}
					everybdCSxs.add(csXs);
					everybdCSys.add(csYs);
				}
				

			}// end the for (int k = 0; k < params.size(); k++) loop
			
			// draw shape in original image
			askImg(pixels, everybdxs, everybdys, unfilteredImage, orName + "_Contour.jpg",
					"yellow", false, drawNumbers);
			askImg(pixels, everybdxs, everybdys, unfilteredImage, orName + "_Contour_White.jpg",
					"red", true, drawNumbers);
			putPointsInFile(everybdxs, everybdys, orName + "_Contour.txt");
			
			if (cs.equals("y") || cs.equals("Y")) {
				askImg(pixels, everybdCSxs, everybdCSys, unfilteredImage, orName + "_Convex_Core.jpg", 
						"yellow", false, drawNumbers);
				askImg(pixels, everybdCSxs, everybdCSys, unfilteredImage, orName + "_Convex_Core_White.jpg", 
						"blue", true, drawNumbers);
				//everybdCSxs and ys = 
				askImg(pixels, everybdCSxs, everybdCSys, unfilteredImage, orName + "_Dist_Graph.jpg", 
						"blue", true, drawNumbers);
				
				putPointsInFile(everybdCSxs, everybdCSys, orName + "_Convex_Core.txt");
			}

			clock.stop();
			clock.printImageInformation(params, filteredFilesInFolder, j, this);
		}////end for loop

		return image;
	}

	/**
	 * Prints center of mass coordinates of each object in the frame.
	 * @param params
	 */
	public void printCOM(ArrayList<ObjectParams> params) {
		System.out.println("Contenter of mass of each object in frame: ");
		for (int i = 0; i < params.size(); i++) {
			params.get(i).printRecord();
			System.out.println("\n");
		}
	}

	/**
	 * Lists what each optoion in the parameter file stands for.
	 * Finally asks for the name of the parameter file we want to use.
	 */
	public void listParamFileOptions() {
		System.out.println("Format of File with Parameters:"
						+ "\n"
						+ "Line 1: Path to unfiltered images folder string\n" 
						+ "Line 2: x coordinate integer\n" 
						+ "Line 3: y coordinate integer\n"
						+ "Line 4: red integer \n"
						+ "Line 5: blue integer\n" 
						+ "Line 6: green integern\n"
						+ "Line 7: endI integer (number of straigth segments per circle\n"
						+ "Line 8: dt double (distance between circles)\n"
						+ "Line 9: circles (y n)\n"
						+ "Line 10: Create file listing all points in image (y n)\n"
						+ "Line 11: Enumerate points on image (y n) \n"
						+ "params of different objects are separated by ...\n"
						+ "comments are everything after ... but on the same line as the ...\n");
		System.out.println("Name of file with parameters\n");
	}

	/**
	 * Calculate where in the image the current point should be drawn. 
	 * In the image,
	 * @param params object containing the final parameters to be drawn. 
	 * @param k - the k-th point will be extracted.
	 * @param points list of (x, y) coordinates
	 * @param c - number of columns
	 * @param r - number of rows
	 */
	public void extractXsAndYs(ArrayList<ObjectParams> params, int k,
			int[][] points, int c, int r) {
		for (int i = 0; i < points.length; i++) {
			
			
			if (i >= params.get(k).xs.size()) {
				params.get(k).xs.add(c - points[i][0]);
				params.get(k).ys.add(r - points[i][1]);
			} else {
				params.get(k).xs.set(i, c- points[i][0]);
				params.get(k).ys.set(i, r - points[i][1]);
			}
			/*if (i >= params.get(k).xs.size()) {
				params.get(k).xs.add(points[i][0]);
				params.get(k).ys.add(points[i][1]);
			} else {
				params.get(k).xs.set(i, points[i][0]);
				params.get(k).ys.set(i, points[i][1]);
			}*/
				
		}
	}

	/**
	 * Displays an error message.
	 */
	public void illegalParams(String unfilteredImage, String orName,
			int xcoordinate, int ycoordinate) {
		System.out.println("Some parameter is not assigned");
		System.out.println("x: " + xcoordinate + "\ny: " + ycoordinate
				+ "\n unfilteredImage: " + unfilteredImage
				+ "\norName: " + orName);
		System.err.print("Program will terminate");
		System.exit(1);
	}

	/**
	 * Get each bug's parameters and save them in a separate obj. Return an
	 * arrayList of these objects. Each one contains a bug's parameters
	 */
	private ArrayList<ObjectParams> initializeParamsFromFile(
			String inputFileName) {

		ArrayList<ObjectParams> things = new ArrayList<ObjectParams>();
		FileReader fileToReadFrom = null;
		String line = "";
		try // tries to read from the file, if the file does not exist goes to
			// the catch statements
		{
			fileToReadFrom = new FileReader(inputFileName);
			BufferedReader bf = new BufferedReader(fileToReadFrom);
			while (bf.ready()) {
				ObjectParams bug = new ObjectParams();

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.unfilteredPath = line;

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.xcoordinate = Integer.parseInt(line);

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.ycoordinate = Integer.parseInt(line);

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.red = Integer.parseInt(line);

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.blue = Integer.parseInt(line);

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.green = Integer.parseInt(line);

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				bug.changeEndI(Integer.parseInt(line));

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.dt = Double.parseDouble(line);

				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.circles = line;
				
				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.concaveSet = line;
				
				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.numbers = line;
				
				line = bf.readLine();
				while (line.startsWith("..."))
					line = bf.readLine();
				line.trim();
				bug.fullShape = line;
				
				things.add(bug);
			}

		}

		catch (FileNotFoundException e)// if the file does not exist
		{
			System.err.println("\tFile was not found readText Method");
		} catch (IOException e)// if the file does not exist
		{
			System.err.println("\tUnexpected IOException. From Main Method");
		}

		finally // reads file into string and closes the file
		{

			try {
				fileToReadFrom.close();
			} catch (FileNotFoundException e) {
				System.err
						.println("\tFile not found at closing time readText Method");
			} catch (IOException e) {
				System.err
						.println("\tUnexpected IO error at closing time readText Method");
			} finally {
				return things;
			}
		}// end of finally clause

	}

	private void paths() {
		JIP foldPath = new JIP();
		data.path = foldPath.folder;

		File fldr = new File(data.path);
		data.pathParent = fldr.getParent();
	}

	/** finds all files in a folder and puts returns a list of them in a folder.
	 * 
	 * @param pathstr absolute path to folder
	 * @return a list of all files in the folder
	 */
	public ArrayList<String> findAllFilesInFolder(String pathstr) {
		int lastSlash = pathstr.lastIndexOf("\\");
		String folder = pathstr.substring(lastSlash + 1);
		// arraylist will contain all files
		ArrayList<String> allfilesInFldr = new ArrayList<String>();

		String files;
		File fldr = new File(pathstr);

		File[] listOfFiles = fldr.listFiles();
		System.out.println("Files in " + folder + " folder: ");
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile()) {
				files = listOfFiles[i].getName();
				allfilesInFldr.add(files);
				System.out.println(files);
			}
		}
		System.out.println("\n");

		return allfilesInFldr;
	}

	/**
	 * Finds the center of mass.
	 * @param xs array of all x coordinates
	 * @param ys array of all y coordinates
	 * @return center of mass as an array with format {x, y}
	 */
	public double[] findCOM(ArrayList<Integer> xs, ArrayList<Integer> ys) {
		double xCenter = 0.0;
		double yCenter = 0.0;
		for (int i = 0; i < xs.size(); i++) {
			xCenter += xs.get(i);
			yCenter += ys.get(i);
		}
		xCenter = xCenter / xs.size();
		yCenter = yCenter / ys.size();
		double[] result = { xCenter, yCenter };
		return result;
	}

	/** Put the contours of the object into a blank image and into the original
	 * image.
	 */
	private void askImg(int[][] originalImg, ArrayList<ArrayList<Integer>> xPts,
			ArrayList<ArrayList<Integer>> yPts, String orImgPath, String orName, String color,
			boolean whiteBackground, boolean numberPoints) {

		// gets the file for the original image and draws the contours into that
		// image
		File orImg = new File(orImgPath);
		if (orImg == null)
			System.out.println("file is null");
		BufferedImage original = data.createBufferedImage(orImg, originalImg, whiteBackground);
		if (original == null)
			System.out.println("original is null");
		
		else {
			System.out.println("initial point that will be drawn: " + xPts.get(0).get(0) + " " +
					yPts.get(0).get(0));
			for (int i = 0; i < xPts.size(); i++) {
				//data.drawShape(original, xPts.get(i), yPts.get(i), color, numberPoints);
				data.drawShape(original,  xPts.get(i), yPts.get(i), color, numberPoints);
			}
		}

		// Encode unfiltered image as a JPEG
		try {
			String saveName = "Experiment_Results/" + orName;
			//putPointsInFile(xPts, yPts, saveName);
			FileOutputStream fos = new FileOutputStream(saveName);
			JPEGImageEncoder jpeg = JPEGCodec.createJPEGEncoder(fos);
			try {
				jpeg.encode(original);
				fos.close();
			} catch (IOException e) {
				System.out.println("IOException");
			}
		} catch (FileNotFoundException e) {
			System.out.println("file was not found");
		}
	}

	private void putPointsInFile(ArrayList<ArrayList<Integer>> xPts, ArrayList<ArrayList<Integer>> yPts,
			String saveName) {
		
		String newExtension = ".txt";
		String removeOldExtension = saveName.substring(0, saveName.lastIndexOf('.'));
		//String fileName = removeOldExtension + newExtension;
		String fileName = "Experiment_Results/" + saveName;
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(fileName));
			writer.write("Number of objects detected: " + xPts.size());
			for (int i = 0; i < xPts.size(); i++) { //number of objects
				writer.write("\n_____________________________\n");
				writer.write("Object #" + (i+1) + "\n");
				for (int j = 0; j < xPts.get(i).size(); j++) { // points per object
					writer.write("" + (j+1) + ": (" + xPts.get(i).get(j) + ", " + 
						yPts.get(i).get(j) + ")\n" );
				}
			}
			//Close writer
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}	
	}

	/* Calculates the val of r for the given image */
	private double calculateR(BufferedImage image) {
		int[][] pts = ImageData.imageToPixels(image);
		int length = pts[0].length;
		int width = pts.length;
		double result = (Math.pow((double) width, 2.0) + Math.pow(
				(double) length, 2.0));
		result = (Math.sqrt(result) / 2);
		return result;
	}

}
