import java.util.ArrayList;



public class ObjectParams {
	public String unfilteredPath;
	public int xcoordinate ;
	public int ycoordinate;
	public int red, blue, green; 
	private int endI; //number of points
	public double dt;
	//public int [] xs;
	public ArrayList <Integer> xs;
	//public int [] ys;
	public ArrayList <Integer> ys;
	public String circles;
	public String concaveSet;
	public String numbers;
	public String fullShape;
	
	public ObjectParams(){
		unfilteredPath = "";
		xcoordinate = 0 ;
		ycoordinate = 0;
		red = 0; blue =  0; green = 0; endI = 0;
		 dt = 0.0;
		circles = "";
		concaveSet = "";
		numbers = "";
		fullShape = "";
	}
	public void printRecord(){
		System.out.println("(x,y) (" + 
				xcoordinate+ ", " + ycoordinate +")\nred " +  red +
				" blue "+ blue + " green " +  green + "\nendI " + endI + "\ndt "  + dt + "\ncircles " + 
				circles + "\nconcave set: " + concaveSet + "\nUnfilteredPath: " + unfilteredPath + 
				"\nnumbers:" + numbers + "\nfull shape:" + fullShape);
	}
	public void changeEndI(int num){
		endI = num;
		//xs = new int [endI];
		//ys = new int [endI];
		xs = new ArrayList<Integer>();
		ys = new ArrayList<Integer>();
	}
	public int getEndI (){
		return endI;
	}
}
