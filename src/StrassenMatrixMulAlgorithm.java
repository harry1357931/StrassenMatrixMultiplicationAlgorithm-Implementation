/* Strassen's Matrix Multiplication Algorithm Implementation 
 * Description: 
 *   1) Multiply two Matrices using Strassen's algorithm
 *      and Computes the Sum of elements of MatC.
 *   2) Computes Product of Two Matrices faster than ordinary algorithm  
 *
 ************************
 * Standard Input Format: Two Matrices
 ************************
        [ 5 0 1 7 ]
		[ -1 11 0 9 ]
		[ 45 -3 8 7 ]
		[ 5 0 -4 7 ]
		
		[ 1 -2 0 7 ]
		[ 1 2 3 7 ]
		[ 5 1 0 4 ]
		[ 0 4 2 3 ]
 **************************     
 * @param MatA  First Matrix
 * @param MatB  Second Matrix
 * @param MatC  To Store Product of Matrix A and B
 * @param SumOfElements Stores Sum of Elements of Matrix produced by product of AB = C
 * @param NumberOfMatrixProducts Stores the number of times product is computed by equation 4.8             
 * 
 * @author Gurpreet Singh
 */
import java.util.Scanner;
import java.util.StringTokenizer;

public class StrassenMatrixMulAlgorithm {
    public static int[][] MatA, MatB, MatC;
	public static StringTokenizer myTokens;
	public static String FileName = "even.txt";
	public static int NumOfRows = NumberOfRows(FileName);       // N = NumOfRows
	public static long SumOfElements=0, NumberOfMatrixProducts=0;
	
	public static void main(String[] args) {
		//SecondaryInput();   // for Standard Input, comment the readFile line(next line) if using this
      	readFile(FileName, NumOfRows);     // Reading Matrices A and B
		MatC = Strassen_Matrix_Multiply(MatA, MatB);
		System.out.println("The Number of Matrix Products using equation 4.8 are "+NumberOfMatrixProducts+"\n");
		for(int i=0; i<MatC.length; i++){    // To find Sum of MatC
			for(int j=0; j<MatC.length; j++){
				SumOfElements += MatC[i][j];
			}
		}
		System.out.println("The Sum of Elements of AB is "+SumOfElements);
		System.out.println("The Product of Matrix A and B is: \n");
		for(int i=0; i< MatC.length; i++){
			for(int j=0; j< MatC.length; j++){
				System.out.print(MatC[i][j]+"  ");
			}
			System.out.println();
		}
	}  // main
	
	public static void SecondaryInput(){
		    MatA = new int[2][2];         
	        MatB = new int[2][2];    
		    Scanner scanner = new Scanner(System.in);
	        String NextInput = scanner.next();
		
	        if(NextInput.equalsIgnoreCase("[")){
		    	int j=0;
		    	while(NextInput.equalsIgnoreCase("[") == true){
		    	   NextInput = scanner.next();
		    	   int count=0;	
		    	   while(NextInput.equalsIgnoreCase("]")== false){
		    		   MatA[j][count++] = Integer.parseInt(NextInput);
		    		   NextInput = scanner.next();
		    	   }
		    	   j++;
		    	   if(j==2) 
 		    		   break;
		    	   NextInput = scanner.next();
		    	} // outer while loop
		    }
		    System.out.println();
		    for(int i=0; i< 2; i++){
		    	for(int j=0; j<2; j++){
		           System.out.print(MatA[i][j]+" ");
		    	}
	            System.out.println();
		    }
		    
		    // Mat B
		    NextInput = scanner.next();
            if(NextInput.equalsIgnoreCase("[")){
		    	
		    	int j=0;
		    	while(NextInput.equalsIgnoreCase("[") == true){
		    	   NextInput = scanner.next();
		    	   int count=0;	
		    	   while(NextInput.equalsIgnoreCase("]")== false){
		    		   MatB[j][count++] = Integer.parseInt(NextInput);
		    		   NextInput = scanner.next();
		    	   }
		    	   j++;
		    	   if(j==2)    
 		    		   break;
		    	   NextInput = scanner.next();
		    	} // outer while loop
		    }
		    
		    System.out.println();
		    for(int i=0; i< 2; i++){
		    	for(int j=0; j<2; j++){
		           System.out.print(MatB[i][j]+" ");
		    	}
	            System.out.println();
		    }
	}		    
	
	public static int[][] Strassen_Matrix_Multiply(int[][] A, int[][] B){
		 int n = A.length;
		 int[][] C = new int[n][n];
		 
		 if(n < Min(NumOfRows/2.0, 400) || n == 1){ 
			 for(int i=0; i<n; i++){
				 for(int j=0; j<n; j++){
					 C[i][j] = 0;
					 for(int k=0; k<n; k++){
						  C[i][j] += A[i][k]*B[k][j];       // Eqn.---> 4.8 Naive Matrix Multiplication	
					 }
			     }
				 NumberOfMatrixProducts++;
			}
		 }// if ends here...
		 else if( n >= Min(NumOfRows/2.0, 400) && n%2==0){    // if n is even and... 
			 
			 int[][] S,T,U,V,W,X,Y,Z;
			 S = new int[n/2][n/2]; T = new int[n/2][n/2]; U = new int[n/2][n/2]; V = new int[n/2][n/2];
			 W = new int[n/2][n/2]; X = new int[n/2][n/2]; Y = new int[n/2][n/2]; Z = new int[n/2][n/2];
			 
			 for(int i=0; i<n; i++){                 // Compute S,T,U,V,W,X,Y,Z
				 for(int j=0; j<n; j++){
					 if(i<n/2 && j<n/2){
						 S[i][j] = A[i][j];
						 W[i][j] = B[i][j];
					 }
					 else if(i<n/2 && j>=n/2){
						 T[i][j-n/2] = A[i][j];
						 X[i][j-n/2] = B[i][j];
					 }
					 if(i>=n/2 && j<n/2){
						 U[i-n/2][j] = A[i][j];
						 Y[i-n/2][j] = B[i][j];
					 }
					 if(i>=n/2 && j>= n/2){
						 V[i-n/2][j-n/2] = A[i][j];
						 Z[i-n/2][j-n/2] = B[i][j];
					 }
			         		 
				 }
			 } // for loop i ends here...
			  
			 
			 int[][] P1 = Strassen_Matrix_Multiply(S, MatAddorSub(X,Z,true)), 
			   P2 = Strassen_Matrix_Multiply(MatAddorSub(S,T,false), Z), 
			   P3 = Strassen_Matrix_Multiply(MatAddorSub(U,V,false), W), 
			   P4 = Strassen_Matrix_Multiply(V, MatAddorSub(Y,W,true)), 
			   P5 = Strassen_Matrix_Multiply(MatAddorSub(S,V,false), MatAddorSub(W,Z,false)), 
			   P6 = Strassen_Matrix_Multiply(MatAddorSub(T,V,true), MatAddorSub(Y,Z,false)), 
			   P7 = Strassen_Matrix_Multiply(MatAddorSub(S,U,true), MatAddorSub(W,X,false));
			  
			 C = MatrixFromSevenProducts(n, P1, P2, P3, P4, P5, P6, P7);		  
		 }
		 else if(n >= Min(NumOfRows/2.0, 400) && n%2!=0 && n > 1){   // if n is odd
			 
			 int[][] PaddedMatA = new int[n+1][n+1], PaddedMatB = new int[n+1][n+1];
			 
			 for(int i=0; i< PaddedMatA.length; i++){
				 for(int j=0; j< PaddedMatA.length; j++){
					if(i< A.length && j< A.length){ 
					   PaddedMatA[i][j] = A[i][j];
					   PaddedMatB[i][j] = B[i][j];
					}
					
					else{
						PaddedMatA[i][j] = 0;
						PaddedMatB[i][j] = 0;
					}
				 }
			 }// for
			 
			 int[][] PaddedProduct = Strassen_Matrix_Multiply(PaddedMatA, PaddedMatB);
			 // UN pad it...and return the unpadded result
             for(int i=0; i< C.length; i++){
        	    for(int j=0; j< C.length; j++){
        		    C[i][j] = PaddedProduct[i][j];    
        	    }  
             } // for			
		 }
		return C;
	}   // Recursive Function
	
	public static int[][] MatrixFromSevenProducts(int n, int[][] P1, int[][] P2,int[][] P3,int[][] P4,int[][] P5,int[][] P6, int[][] P7 )
	{   int[][] C = new int[n][n];
	    
	    for(int i=0; i<n; i++){                       // Compute C
			 for(int j=0; j<n; j++){
				 if(i<n/2 && j<n/2){                  // Quad1
					 C[i][j] = P4[i][j] + P5[i][j] + P6[i][j] - P2[i][j] ;				 
				 }
				 else if(i<n/2 && j>=n/2){            // Quad2
					 C[i][j] = P1[i][j-n/2] + P2[i][j-n/2];
				 }
				 if(i>=n/2 && j<n/2){                 // Quad3
					 C[i][j] = P4[i-n/2][j] + P3[i-n/2][j];
				 }
				 if(i>=n/2 && j>= n/2){               // Quad4
					 C[i][j] = P1[i-n/2][j-n/2] + P5 [i-n/2][j-n/2] - P7[i-n/2][j-n/2] - P3[i-n/2][j-n/2];
				 }
			 }
		 } // for loop i ends here...
		return C;
	} 
	
	public static int[][] MatAddorSub(int[][] Mat1, int[][] Mat2, boolean Subtract){
		int[][] Mat3 = new int[Mat1.length][Mat1.length];
		
		for(int i=0; i<Mat1.length; i++){
			 for(int j=0; j<Mat1.length; j++){
				if(Subtract == true) 
				   Mat3[i][j] = Mat1[i][j] - Mat2[i][j];
				else
				   Mat3[i][j] = Mat1[i][j] + Mat2[i][j];
			 }
		 }
		
		return Mat3;
	}

	public static double Min(double a, int b)
	{
	   if(a>=b)	
		  return b;
	   
	   return a;
	}
	
    public static void readFile(String fileName, int MatrixDimension)        // length of array to be filled  
	{   
	    MatA = new int[MatrixDimension][MatrixDimension];         
	    MatB = new int[MatrixDimension][MatrixDimension];
	    TextFileInput tfi = new TextFileInput(fileName);    
		String line = tfi.readLine(); 
	    
		int count=0;
	    while(line!=null)
		{   
		    int j=0;
		    if(line.equalsIgnoreCase("")){ 
		    	line=tfi.readLine();
	    	    continue;
		    }else{
		    	myTokens = new StringTokenizer(line,"[");
	         	myTokens = new StringTokenizer(myTokens.nextToken(), "]");
	        	myTokens = new StringTokenizer(myTokens.nextToken(), " ");
		    	
		    	if(count < NumOfRows){
				    while(myTokens.hasMoreTokens()){
					   MatA[count][j] =  Integer.parseInt(myTokens.nextToken());
					   j++; 
					}
				    count++;
		    	}
		    	else{
		    		while(myTokens.hasMoreTokens()){ 
					   MatB[count-NumOfRows][j] =  Integer.parseInt(myTokens.nextToken());
					   j++; 
					}
					count++;
		    	 }
		     }
			      
			line=tfi.readLine();
		}  // while loop ends here
    } // Read File method ends here
	
	public static int NumberOfRows(String fileName)        // length of array to be filled  
	{   
	    int CountLine=0;
	    TextFileInput tfi = new TextFileInput(fileName); 
		String line = tfi.readLine(); 
	    while(line!=null)
		{   if(line.equalsIgnoreCase(""))
	    	    break;
	       
	    	line=tfi.readLine();
			CountLine++;
		}
	    return CountLine;  
		
    }  // Read File method ends here	    
}     // Class StrassenMatrixMultiplication
