import java.awt.Font;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import javax.imageio.ImageIO;


public class ImageData {
	/*This class takes colors from an image, reads and writes to 
	 * an image an draws lines in an image. 
	 */
	public String path;
	public String pathParent;
	static int green;
	static int blue;
	static int red;

	public ImageData(String path, String pathParent) {
		this.path = path;
		this.pathParent = pathParent;
	}

	/**
	 * Extracts the red component of a pixel.
	 * 
	 * @param pixel an integer pixel
	 * 
	 * @return the red component [0-255] of the pixel.
	 */
	public static int getRed(int pixel) {
		return pixel >> 16 & 0xff;
	}

	/**
	 * Extracts the green component of a pixel.
	 * 
	 * @param pixel an integer pixel
	 * 
	 * @return the green component [0-255] of the pixel.
	 */
	public static int getGreen(int pixel) {
		return pixel >> 8 & 0xff;
	}

	/**
	 * Extracts the blue component of a pixel.
	 * 
	 * @param pixel
	 *            an integer pixel
	 * 
	 * @return the blue component [0-255] of the pixel.
	 */
	public static int getBlue(int pixel) {
		return pixel & 0xff;
	}

	/**
	 * Constructs a pixel from RGB components.
	 * 
	 * @param red
	 *            the red component [0-255]
	 * @param green
	 *            the green component [0-255]
	 * @param blue
	 *            the blue component [0-255]
	 * 
	 * @return the packed integer pixel.
	 */
	public static int makePixel(int red, int green, int blue) {
		return (red & 0xff) << 16 | (green & 0xff) << 8 | (blue & 0xff);
	}

	/**
	 * Converts a two-dimensional array of image data to a BufferedImage.
	 * 
	 * @param pixels
	 *            the source two-dimensional array of image data
	 * 
	 * @return a BufferedImage representing the source image.
	 * 
	 * @throws IllegalArgumentException
	 *             if the source image is ill-defined.
	 */
	public static BufferedImage pixelsToImage(int[][] pixels)
			throws IllegalArgumentException {
		if (pixels == null) {
			throw new IllegalArgumentException();
		}
	
		int width = pixels[0].length;
		int height = pixels.length;
		BufferedImage image = new BufferedImage(width, height,
				BufferedImage.TYPE_INT_RGB);
		for (int row = 0; row < height; row++) {
			image.setRGB(0, row, width, 1, pixels[row], 0, width);
		}
		return image;
	}

	/**
	 * Converts a BufferedImage to a two-dimensional array of image data.
	 * 
	 * @param image
	 *            the source BufferedImage
	 * 
	 * @return a two-dimensional array of image data representing the source
	 *         image
	 * 
	 * @throws IllegalArgumentException
	 *             if the source image is ill-defined.
	 */
	public static int[][] imageToPixels(BufferedImage image)
			throws IllegalArgumentException {
		if (image == null) {
			throw new IllegalArgumentException();
		}
	
		int width = image.getWidth();
		int height = image.getHeight();
		int[][] pixels = new int[height][width];
		for (int row = 0; row < height; row++) {
			image.getRGB(0, row, width, 1, pixels[row], 0, width);
		}
		return pixels;
	}

	/**
	 * Makes a buffered image when given the pathname of the original image
	 */
	public BufferedImage createBufferedImage(File file, int [][] originalImg, 
			boolean whiteBackground) {
		BufferedImage img;
		try {
			img = ImageIO.read(file);
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		if (img == null)
			System.out.println("IMAGE IS NULL. THAT'S THE PROBLEM");
		if (img.getType() == BufferedImage.TYPE_INT_RGB) {
			if (whiteBackground) img = makeBackgroundWhite (img, originalImg);
			return img;
		}
	
		BufferedImage nImg = new BufferedImage(img.getWidth(null),
				img.getHeight(null), BufferedImage.TYPE_INT_RGB);
		
		if (whiteBackground) img = makeBackgroundWhite (img, originalImg);
		
		nImg.createGraphics().drawImage(img, 0, 0, null);
		return nImg;
	}

	public BufferedImage drawLine(BufferedImage image, int x1, int y1, int x2, int y2) {
		Graphics g = image.getGraphics();
		g.setColor(java.awt.Color.BLUE);
		g.drawLine(x1, y1, x2, y2);
		return image;
	}
	
	/**
	 * Makes the background white
	 * @param img the image to make the background white
	 * @param originalImg the dimensions of the original image
	 * @return
	 */
	public BufferedImage makeBackgroundWhite (BufferedImage img, int [][] originalImg) {
		Graphics g = img.getGraphics();
		g.setColor(java.awt.Color.WHITE);
		g.clearRect(0, 0, originalImg[0].length, originalImg.length);
		g.fillRect(0, 0, originalImg[0].length, originalImg.length);
		
		return img;
	}

	/** Draws a connected line segment between the given points 
	 * @param image - the pixels of the image. (0, 0) is in the upper left hand corner.
	 * @param xs TODO
	 * @param ys TODO
	 * @param color TODO*/
	public BufferedImage drawShape(BufferedImage image, ArrayList<Integer> xs, ArrayList<Integer> ys,
			String color, boolean numberPoints) {
		Graphics g = image.getGraphics();
		if (color.equals("blue"))
			g.setColor(java.awt.Color.BLUE);
		else if (color.equals("red"))
			g.setColor(java.awt.Color.RED);
		else if (color.equals("yellow"))
			g.setColor(java.awt.Color.YELLOW);
		else if (color.equals("green"))
			g.setColor(java.awt.Color.GREEN);
	
		//This seems to flip the picture along the x = y line. 
		int []xarr = arrListToArray(xs); 
		int []yarr = arrListToArray(ys);
		g.drawPolygon(xarr, yarr, xs.size());
		//g.drawPolygon(xs, ys, xs.size());
		
		if (numberPoints == true) {
		//number the points
		g.setColor(java.awt.Color.GREEN);
		Font f = new Font("serif", Font.PLAIN, 15);
		g.setFont(f);
		for (int i =0; i<xs.size(); i++) {
			int j = i+1; 
			char[] data = ("" + j).toCharArray();
			//char []data = {'c', 'd'};
			g.drawChars(data, 0, data.length, xs.get(i), ys.get(i));
		}
		}
		
		
	
		return image;
	}

	public static int[] arrListToArray(ArrayList<Integer> arrList) {
		int[] arr = new int [arrList.size()];
		for (int i =0; i<arrList.size(); i++) {
			arr[i] = arrList.get(i);
		}
		return arr;
	}

	/**
	 * if you answered y or Y to draw circles question, 
	 * this code will draw the circles which the program uses 
	 * to sample the colors at those positions.
	 * 
	 * @param imageEffect TODO
	 * @param params
	 * @param xcoordinate
	 * @param ycoordinate
	 * @param constant
	 * @param k
	 */
	public void drawCircles(BoundaryFinder imageEffect, ArrayList<ObjectParams> params,
			int xcoordinate, int ycoordinate, double constant, int k) {
		/*int[][] points;
		ArrayList<Integer> xlist = new ArrayList<Integer>();
		ArrayList<Integer> ylist = new ArrayList<Integer>();
		for (double curK = 0; curK <= BoundaryFinder.endK; curK = curK + BoundaryFinder.dt) {
	
			
			points = imageEffect.calculate(BoundaryFinder.dt, curK, BoundaryFinder.endI, xcoordinate,
					ycoordinate, constant);
			for (int i = 0; i < points.length; i++) {
				xlist.add(points[i][0]);
				ylist.add(points[i][1]);
			}
		}
	
		params.get(k).xs = new ArrayList <Integer>();
		params.get(k).ys = new ArrayList <Integer>();
	
		// divide up the points into an array of xs and ys
		// mid num = endK
		for (int i = 0; i < xlist.size(); i++) {
			if (i>= params.get(k).xs.size()) {
			params.get(k).xs.add(xlist.get(i));
			params.get(k).ys.add(ylist.get(i));
			} else {
				params.get(k).xs.set(i, xlist.get(i));
				params.get(k).ys.set(i, ylist.get(i));
			}
		}
	*/	
	}
}