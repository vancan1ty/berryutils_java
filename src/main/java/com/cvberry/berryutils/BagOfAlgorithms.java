package com.cvberry.berryutils;

import java.util.function.DoublePredicate;
import java.util.function.BiPredicate;
import java.util.function.DoubleBinaryOperator;

import java.util.Arrays;
import java.util.Random;

/** Currell Berry mathematical utility methods
 *  written for CS 4635 project 1, because we aren't 
 *  allowed to use a math library!
 *  */
public class BagOfAlgorithms {

    /** first comes a variety of algorithms for manipulating 2d arrays (matrices) */
    /*public static double[][] b_grayscaleToWBMatrix(double[][] gMatrix) {
        double mean = mat_mean(gMatrix);
        boolean[][] results = b_apply_to_binary(gMatrix, (x)-> (x < mean));
    }*/

    public static int[] mat_shape(double[][] matrix) {
       int height = matrix.length; 
       int width = matrix[0].length;
       return new int[]{height,width};
    }

    public static int[] mat_shape(boolean[][] matrix) {
       int height = matrix.length; 
       int width = matrix[0].length;
       return new int[]{height,width};
    }

    public static double mat_sum(double[][] matrix) {
        double sum = 0;
        int[] shapeArr = mat_shape(matrix);
        for (int i = 0; i < shapeArr[0]; i++) {
            for (int j = 0; j < shapeArr[1]; j++) {
                sum += matrix[i][j]; 
            }
        }
        return sum;
    }

    public static int mat_sum(boolean[][] matrix) {
        int sum = 0;
        int[] shapeArr = mat_shape(matrix);
        for (int i = 0; i < shapeArr[0]; i++) {
            for (int j = 0; j < shapeArr[1]; j++) {
                if (matrix[i][j]) {
                    sum += 1;
                }
            }
        }
        return sum;
    }

    public static double mat_mean(double[][] matrix) {
        double sum = mat_sum(matrix);
        int[] shapeArr = mat_shape(matrix);
        int totEntries = shapeArr[0]*shapeArr[1];
        return sum/totEntries;
    }

    public static boolean[][] mat_AND(boolean[][] m1, boolean[][] m2) {
        return mat_binary_op(m1,m2,(Boolean a, Boolean b) -> a & b);
    }

    public static boolean[][] mat_OR(boolean[][] m1, boolean[][] m2) {
        return mat_binary_op(m1,m2,(Boolean a, Boolean b) -> a | b);
    }

    public static boolean[][] mat_binary_op(boolean[][] m1, boolean[][] m2, BiPredicate<Boolean,Boolean> whatToDo) {
        assert(Arrays.equals(mat_shape(m1),mat_shape(m2)));
        int[] shapeArr = mat_shape(m1);
        boolean[][] out = new boolean[shapeArr[0]][shapeArr[1]];
        for (int i = 0; i < shapeArr[0]; i++) {
            for (int j = 0; j < shapeArr[1]; j++) {
                out[i][j] = whatToDo.test(m1[i][j],m2[i][j]);
            }
        }
        return out;
    }

    public static double[][] mat_add(double[][] m1, double[][] m2) {
        return mat_binary_op(m1, m2, (double a, double b) -> a+b);
    }

    public static double[][] mat_sub(double[][] m1, double[][] m2) {
        return mat_binary_op(m1, m2, (double a, double b) -> a-b);
    }
    
    public static double[][] mat_elementwiseMult(double[][] m1, double[][] m2) {
        return mat_binary_op(m1, m2, (double a, double b) -> a*b);
    }

    
    /**
     * use naive matrix multipication algorithm, n^3 time
     */
    public static double[][] mat_mult(double[][] m1, double[][] m2) {
        int[] shp1 = mat_shape(m1);
        int[] shp2 = mat_shape(m2);
        assert(shp1[1]==shp2[0]);

        double[][] out = new double[shp1[0]][shp2[1]];
        for (int m2col = 0; m2col < shp2[1]; m2col++) {
            for (int m1row = 0; m1row < shp1[0]; m1row++) {
                double value = 0;
                for (int index = 0; index < shp1[1]; index++) {
                    value += m1[m1row][index]*m2[index][m2col];
                }
                out[m1row][m2col] = value;
            }
        }
        return out;
    }

    public static double[] mat_vec_mult(double[][] mat, double[] vec) {
        int[] shp = mat_shape(mat);
        assert(shp[1] == vec.length);
        double[] out = new double[shp[0]];
        for (int m1row = 0; m1row < shp[0]; m1row++) {
            double value = 0;
            for (int index = 0; index < shp[1]; index++) {
                value += mat[m1row][index]*vec[index];
            }
            out[m1row] = value;
        }
        return out;
    }

    public static boolean mat_equals(double[][] m1, double[][] m2) {
        return Arrays.deepEquals(m1,m2); //CB I think this should work.
    }

    public static boolean mat_equals(boolean[][] m1, boolean[][] m2) {
        return Arrays.deepEquals(m1,m2); //CB I think this should work.
    }

    public static boolean mat_equals(double[][] m1, double[][] m2, double tolerance) {
        int[] shp1 = mat_shape(m1);
        if (!Arrays.equals(shp1,mat_shape(m2))) {
            return false;
        }
        for (int m =0; m < shp1[0]; m++) {
            for (int n = 0; n < shp1[1]; n++) {
                if (Math.abs(m2[m][n]-m1[m][n])>tolerance) {
                    //System.out.format("b: %f, a: %f, tolerance: %f\n",m2[m][n],m1[m][n],tolerance);
                    return false;
                }
                if(Double.isNaN(m2[m][n]) || Double.isNaN(m1[m][n])) {
                    return false;
                }
            }
        }

        return true;
    }

    public static boolean vec_equals(double[] v1, double[] v2, double tolerance) {
        if (v1.length != v2.length) {
            return false;
        }

        for (int i =0; i < v1.length; i++) {
            if (Math.abs(v1[i]-v2[i]) > tolerance) {
                return false;
            }
        }

        return true;
    }

    public static double[][] mat_binary_op(double[][] m1, double[][] m2, DoubleBinaryOperator whatToDo) {
        assert(Arrays.equals(mat_shape(m1),mat_shape(m2)));
        int[] shapeArr = mat_shape(m1);
        double[][] out = new double[shapeArr[0]][shapeArr[1]];
        for (int i = 0; i < shapeArr[0]; i++) {
            for (int j = 0; j < shapeArr[1]; j++) {
                out[i][j] = whatToDo.applyAsDouble(m1[i][j],m2[i][j]);
            }
        }
        return out;
    }

    public static boolean[][] mat_apply_predicate(double[][] mat, DoublePredicate predicate) {
        int[] shapeArr = mat_shape(mat);
        boolean[][] out = new boolean[shapeArr[0]][shapeArr[1]];
        for (int i = 0; i < shapeArr[0]; i++) {
            for (int j = 0; j < shapeArr[1]; j++) {
                out[i][j] = predicate.test(mat[i][j]);
            }
        }
        return out;
    }

    public static double[][] mat_fliplr(double[][] mat) {
        int[] shapeArr = mat_shape(mat);
        double[][] out = new double[shapeArr[0]][shapeArr[1]];
        int widthindex = mat[0].length-1;
        for (int i = 0; i < mat.length; i++) {
            int lindex=0;
            int rindex=widthindex;
            while ((rindex >=0) && (lindex <= widthindex)) {
                out[i][lindex] = mat[i][rindex];
                rindex--;
                lindex++;
            }
        }
        return out;
    }

    public static double[][] mat_flipud(double[][] mat) {
        int[] shapeArr = mat_shape(mat);
        double[][] out = new double[shapeArr[0]][shapeArr[1]];
        int heightindex = mat.length-1;
        int topindex=0;
        int bottomindex=heightindex;
        while ((bottomindex >=0) && (topindex <= heightindex)) {
           for (int j = 0; j < shapeArr[1]; j++) {
            out[topindex][j] = mat[bottomindex][j];
           }
           topindex++;
           bottomindex--;
        }
        return out;
    }

    public static double[][] mat_transpose(double[][] mat) {
        int[] shapeArr = mat_shape(mat);
        double[][] transpose = new double[shapeArr[1]][shapeArr[0]];
        for (int c = 0; c < shapeArr[0]; c++) {
            for(int d = 0 ; d < shapeArr[1] ; d++ ) {
                transpose[d][c] = mat[c][d];
            }
        }
        return transpose;
    }

    public static double[][] mat_print(double[][] mat) {
        mat_print(mat,"");
        return mat;
    }

    public static double[][] mat_print(double[][] mat, String message) {
        System.out.println(message);
        int[] shp = mat_shape(mat);
        for (int m = 0; m < shp[0]; m++) {
            for (int n = 0; n < shp[1]; n++) {
                System.out.format("%f ", mat[m][n]);
            }
            System.out.println();
        }
        return mat;
    }

    public static double[] vec_print(double[] vec) {
        vec_print(vec,"");
        return vec;
    }


    public static double[] vec_print(double[] vec, String message) {
        System.out.println(message);
        for (int m = 0; m < vec.length; m++) {
            System.out.format("%f ", vec[m]);
            System.out.println();
        }
        return vec;
    }


    public static double vec_dot_product(double[] a, double[] e) {
        assert(a.length==e.length);
        double out = 0;
        for(int i = 0; i < a.length; i++) {
            out += a[i]*e[i]; 
        }
        return out;
    }

    public static double[] vec_proj_a_onto_e(double[] a, double[] e) {
        assert(a.length==e.length);
        double ea = vec_dot_product(e,a);
        double ee = vec_dot_product(e,e);
        return vec_scalar_mult(ea/ee,e);
    }

    public static double[] vec_scalar_vec_op(
                                               double a,
                                               double[] b,
                                               DoubleBinaryOperator whatToDo) {
        double[] out = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            out[i] = whatToDo.applyAsDouble(a,b[i]);
        }
        return out;
    }

    public static double[] vec_scalar_mult(double a, double[] b) {
        return vec_scalar_vec_op(a,b, (double x, double y) -> x*y);
    }

    public static double[] vec_scalar_sub(double a, double[] b) {
        return vec_scalar_vec_op(a,b, (double x, double y) -> x-y);
    }

    public static double[] vec_scalar_add(double a, double[] b) {
        return vec_scalar_vec_op(a,b, (double x, double y) -> x+y);
    }

    public static double[] vec_scalar_divide(double a, double[] b) {
        return vec_scalar_vec_op(a,b, (double x, double y) -> x/y);
    }

    public static double[] vec_binary_elementwise_op(
                                               double[] a,
                                               double[] b,
                                               DoubleBinaryOperator whatToDo) {
        assert(a.length==b.length);
        double[] out = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            out[i] = whatToDo.applyAsDouble(a[i],b[i]);
        }
        return out;
    }

    public static double[] vec_mult(double[] a, double[] b) {
        return vec_binary_elementwise_op(a,b, (double x, double y) -> x*y);
    }

    public static double[] vec_sub(double[] a, double[] b) {
        return vec_binary_elementwise_op(a,b, (double x, double y) -> x-y);
    }

    public static double[] vec_add(double[] a, double[] b) {
        return vec_binary_elementwise_op(a,b, (double x, double y) -> x+y);
    }

    public static double[] vec_divide(double[] a, double[] b) {
        return vec_binary_elementwise_op(a,b, (double x, double y) -> x/y);
    }


    public static double vec_norm(double[] vec) {
        double accum = 0;
        for(int i = 0; i < vec.length; i++) {
            accum += vec[i]*vec[i];
        }
        return Math.sqrt(accum);
    }

    public static double[] mat_grab_column(double[][] mat, int col_index) {
        int[] shp = mat_shape(mat);
        double[] out = new double[shp[0]];
        for (int m = 0; m < shp[0]; m++) {
                out[m] = mat[m][col_index];
        }
        return out;
    }

    public static double[][] mat_splice_out_column(double[][] mat, int col) {
        int[] shp = mat_shape(mat);
        double[][] out = new double[shp[0]][shp[1]-1];
        //copy over cols before col
        for (int n = 0; n < col; n++) {
            for (int m = 0; m < shp[0]; m++) {
                out[m][n]=mat[m][n];
            }
        }
        //copy over cols after col
        for (int n = col+1; n < shp[1]; n++) {
            for (int m = 0; m < shp[0]; m++) {
                out[m][n-1]=mat[m][n];
            }
        }
        return out;
    }

    public static double[][] mat_splice_out_row(double[][] mat, int row) {
        int[] shp = mat_shape(mat);
        double[][] out = new double[shp[0]-1][shp[1]];
        //copy over rows before row
        for (int m = 0; m < row; m++) {
            for (int n = 0; n < shp[1]; n++) {
                out[m][n]=mat[m][n];
            }
        }
        //copy over rows after row
        for (int m = row+1; m < shp[0]; m++) {
            for (int n = 0; n < shp[1]; n++) {
                out[m-1][n]=mat[m][n];
            }
        }
        return out;
    }

    public static double[][][] qr_gschmidt(double[][] mat) {
        return qr_gschmidt(mat, 0.000001);
    }

    /** based on https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Numerical_stability
        @param mat the matrix to decompose
        @param singularThreshold the threshold for the norm of a vector under which to treat it as 0
     */
    public static double[][][] qr_gschmidt(double[][] mat, double singularThreshold) {
        int[] shpOrig = mat_shape(mat);

        double[][] Q = new double[shpOrig[0]][shpOrig[0]];
        int[] shpQ = mat_shape(Q);
        double[][] R = new double[shpOrig[0]][shpOrig[1]];

        int rank = 0;
        int qoffset = 0; //tracks how many columns back to track Q versus mat.  When a dependent column is found in mat,
        //this number is incremented
        for (int n = 0; n < shpOrig[1]; n++) {
            int nq = n-qoffset;
            //first, get original column
            double[] ocol = mat_grab_column(mat,n);
            //next, calculate evec: ocol less the projections of ocol onto each of the previous Q vectors
            double[] uk_est = Arrays.copyOf(ocol, ocol.length);

            if (nq < shpQ[1]) {
                for (int n1 = 0; n1 < nq; n1++) {
                    double[] currentQCol = mat_grab_column(Q,n1);
                    double Qcol_ocol_dp = vec_dot_product(ocol,currentQCol);
                    R[n1][n] = Qcol_ocol_dp;

                    //below is standard gram-schmidt
                    //uk_est = vec_sub(uk_est, vec_scalar_mult(Qcol_ocol_dp,currentQCol));
                    //below is modified gram-schmidt.  better stability supposedly.
                    uk_est = vec_sub(uk_est, vec_proj_a_onto_e(uk_est,currentQCol));
                }

                double uk_est_norm = vec_norm(uk_est);
                if (uk_est_norm < singularThreshold) {
                    //then we have a 0 vector, splice out from Q and R
                    //Q = mat_splice_out_column(Q,nq);
                    //R = mat_splice_out_row(R,nq);
                    //shpQ[1] = shpQ[1]-1;
                    qoffset+=1;
                    continue;
                }
                R[nq][n] = uk_est_norm;

                double[] evec = vec_scalar_mult(1.0/uk_est_norm, uk_est);
                rank +=1;
                for (int m = 0; m < shpOrig[0]; m++) {
                    Q[m][nq] = evec[m];
                }
            } else {
                for (int n1 = 0; n1 < shpQ[1]; n1++) {
                    double[] currentQCol = mat_grab_column(Q,n1);
                    double Qcol_ocol_dp = vec_dot_product(ocol,currentQCol);
                    R[n1][n] = Qcol_ocol_dp;
                }
            }
        }
        //now remove the entries which didn't get populated because of non-full rank
        for (int n = rank; n < shpQ[0]; n++) {
           Q = mat_splice_out_column(Q, rank);
           R = mat_splice_out_row(R, rank);
        }

        double[][][] out = new double[2][shpOrig[0]][shpOrig[1]];
        out[0] = Q;
        out[1] = R;
        return out;
    }

    /**
     * solves for x in R*x=y, where R is an upper triangular matrix.
     * can handle under-determined values of R (ie mxn, n > m)
     *
     * R: m <= n
     */
    public static double[] mat_backSubstUpperTriangularMatrix(double[][] R, double[] y) {
        int[] shp = mat_shape(R);
        double[][] tR = new double[shp[0]][shp[0]];
        for (int m = 0; m < shp[0]; m++) {
            for (int n = 0; n < shp[0]; n++) {
                tR[m][n] = R[m][n];
            }
        }

        int[] tShp = mat_shape(tR);

        double[] out = new double[shp[1]];

        //now, do back substitution
        for (int level = y.length-1; level >= 0; level--) {
            double lans = y[level];
            double subtractor = 0;
            for (int n = level+1; n<tShp[1]; n++) {
                subtractor += tR[level][n]*out[n];
            }
            out[level]=(lans-subtractor)/tR[level][level];
        }

        return out;
    }

    //CB TODO fix: get the system solver working.
    /** by default throw an error if the system is under or over determined */
    public static double[] mat_solveSystem(double[][] A, double[] x, double[] b) {
        return mat_solveSystem(A,x,b,false);
    }

    public static double[] mat_solveSystem(double[][] A, double[] x, double[] b, boolean lstSquaresSolution) {
        return x;
    }

    //next up -- add least squares solver support.
    //solves QRx=b, where QR is composed of independent vectors.
    //uses QR, because that is what I have coded.
    //will return a solution to underdetermined systems.
    public static double[] _mat_solveSystemNoLstSq(double[][] Q, double[][] R, double[] b) {
        int rank = mat_shape(Q)[1];
        int aheight = mat_shape(Q)[0]; //number of independent equations
        int awidth = mat_shape(R)[1];
        int num_unknowns = awidth;
        double[] out = new double[awidth];//but this gets replaced below

        if(rank < num_unknowns) {
            //fewer independent equations than unknowns
            //underdetermined, we use QR to give us one (of infinite) valid answers.
            double[] Qtb = mat_vec_mult(mat_transpose(Q),b);
            out = mat_backSubstUpperTriangularMatrix(R,Qtb);
        } else if (rank == num_unknowns) {
            if(aheight == num_unknowns) {
                //then we have a full rank system with unique solution.
                double[] Qtb = mat_vec_mult(mat_transpose(Q),b);
                out = mat_backSubstUpperTriangularMatrix(R,Qtb);
            } else if (aheight > num_unknowns) {
                //then we have an OVERdetermined system, more equations than unknowns.
            } else {
                //not possible. aheight cannot be greater than rank.
                throw new RuntimeException("it is impossible to get here...");
            }
        } else if (rank > num_unknowns) {
            //I believe this is impossible.
            throw new RuntimeException("it is impossible to get here...");
        }
        return out;
    }

    public static double[][] mat_rand(int height, int width) {
        Random random = new Random();
        long seed = random.nextLong();
        return mat_rand(height,width,seed);
    }

    public static double[][] mat_rand(int height, int width, long seed) {
        Random random = new Random(seed);
        double[][] out = new double[height][width];
        for (int m = 0; m < height; m++) {
            for (int n = 0; n < width; n++) {
                out[m][n] = random.nextDouble();
            }
        }
        return out;
    }

    public static double[] vec_rand(int length) {
        Random random = new Random();
        long seed = random.nextLong();
        return vec_rand(length,seed);
    }

    public static double[] vec_rand(int length, long seed) {
        Random random = new Random(seed);
        double[] out = new double[length];
        for (int m = 0; m < length; m++) {
            out[m] = random.nextDouble();
        }
        return out;
    }

    public static double[][] mat_get_sub_matrix(double[][] mat, int m_start, int m_end, int n_start, int n_end) {
       double[][] out = new double[m_end-m_start][n_end-n_start];
       for (int m = m_start; m < m_end; m++) {
           for (int n = n_start; n < n_end; n++) {
              out[m-m_start][n-n_start] = mat[m][n];
           }
       }
        return out;
    }


    //CB Graph stuff!
    //----------------------------------------------
    List<Integer> search_backwards(int[] prevarr, int index) {
        List<Integer> out = new LinkedList<>(); 
        int cursor = index;
        int prev = prevarr[cursor];
        while(cursor!=-2) {
            out.add(cursor);
            int temp = prevarr[cursor];
            prev = cursor;
            cursor = temp;
        }
        Collections.reverse(out);
        out.reverse();
        return out;
    }

    List<Integer> graph_search(List<List<Integer>> adjacents, int start, int end, Queue mqueue) {
        int[] prevarr = new int[adjacents.size()];
        for (int i = 0; i < adjacents.size(); i++) {
            prevarr[i]=-1;
        }
        prevarr[0]=-2;
        mqueue.add(start);
        while(!mqueue.isEmpty()) {
            int item = mqueue.remove();
            if (item==end) {
                return search_backwards(prevarr,end);
            } else {
                for (Integer a : adjacents.get(item)) {
                    if (prevarr[a]==-1) { //unvisited
                        prevarr[a] = item;
                        mqueue.add(a);
                    }
                }
            }
        }

        return null;
    }


//CB BigInteger example
import java.io.*;
import java.math.BigInteger;

public static class Base {
    public static void main(String [] args) {

        // This will reference one line at a time
        String line = null;

        try {
            BufferedReader bufferedReader = 
                new BufferedReader(new InputStreamReader(System.in));

            //discard first line
            bufferedReader.readLine();
            while((line = bufferedReader.readLine()) != null) {
                String validBaseStr = "";

                String[] segments = line.split(" ");
                String num1 = segments[0];
                String operator = segments[1];
                String num2 = segments[2];
                String resultnum = segments[4];

                for (int base = 1; base <= 36; base++) {
                    if (checkIfExpressionValidInBase(num1, num2, operator, resultnum, base))
                        {
                            if (base == 36) {
                                validBaseStr += "0";
                            } else {
                                validBaseStr += Integer.toString(base,36);
                            }
                        }
                }


                if(validBaseStr.length() == 0) {
                    System.out.println("invalid");
                } else {
                    System.out.println(validBaseStr);
                }

            }

            // Always close files.
            bufferedReader.close();         
        }
        catch(IOException ex) {
             ex.printStackTrace();
        }
    }

    //System.out.format("kval of 31: %d", computeKValueForNum(31));

    // public static void main(String[] args) {
    //     String tcase = "111";
    //     BigInteger num1 = new BigInteger(tcase, 2);
    //     System.out.println(java.lang.Character.MAX_RADIX);
    //     System.out.println(java.lang.Character.MIN_RADIX);

    //     System.out.format("11 valid %b\n", checkIfBase1Valid("11"));
    //     System.out.format("10 valid %b\n", checkIfBase1Valid("10"));
    //     System.out.println(num1);
    // }


    public static boolean checkIfExpressionValidInBase(String num1, String num2, String operand, String result, int base) {
        if (!checkIfStrValidInBase(num1, base) ||
            !checkIfStrValidInBase(num2, base) ||
            !checkIfStrValidInBase(result, base)) {
            return false;
        }

        if (base==1) {

            //check if computations add up
            int a = num1.length();
            int b = num2.length();
            int c = result.length();
            boolean valid;
            switch(operand) {
            case("+"): 
                valid = (a+b==c);
                break;
            case("-"): 
                valid = (a-b==c);
                break;
            case("*"): 
                valid = (a*b==c);
                break;
            case("/"): 
                valid = ((a/b==c) && (a%b==0));
                break;
            default:
                throw new RuntimeException("what!");
            }
            return valid;
        } else {
            //check if computations add up
            BigInteger a = new BigInteger(num1, base);
            BigInteger b = new BigInteger(num2, base);
            BigInteger c = new BigInteger(result, base);
            boolean valid;
            switch(operand) {
            case("+"): 
                valid = (a.add(b).equals(c));
                break;
            case("-"): 
                valid = (a.subtract(b).equals(c));
                break;
            case("*"): 
                valid = (a.multiply(b).equals(c));
                break;
            case("/"): 
                valid = (a.divide(b).equals(c) && a.remainder(b).equals(BigInteger.ZERO));
                break;
            default:
                throw new RuntimeException("what!");
            }

            return valid;
            }
        }
    

    public static boolean checkIfStrValidInBase(String str, int base) {
        boolean valid = true;
        if (base == 1) {
            valid = checkIfBase1Valid(str);
        } else {
            try {
                BigInteger i = new BigInteger(str, base);
            } catch (Exception e) {
                valid = false;
            }
        }

        return valid;
    }

    public static boolean checkIfBase1Valid(String str) {
        boolean valid = true;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) != '1') {
                valid = false;
                break;
            }
        }
        return valid;
    }

}

import java.util.*;
public static class Generator {

    public static void main(String[] args) {
        int numToPrint = Integer.parseInt(args[0]);
        Random random = new Random();
        System.out.println(numToPrint);
        for (int i = 0; i < numToPrint-1; i++) {
            int nextNum = Math.round(random.nextFloat()*2000);
            System.out.format("%d ", nextNum);
        }
        int nextNum = Math.round(random.nextFloat()*2000);
        System.out.format("%d\n",nextNum);
    }

    // public static void main(String[] args) {
    //     System.out.println('9');
    //     System.out.println('8');
    // }

}

}
