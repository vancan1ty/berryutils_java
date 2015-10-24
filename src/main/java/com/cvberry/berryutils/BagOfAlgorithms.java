package com.cvberry.berryutils;

import java.util.function.DoublePredicate;
import java.util.function.BiPredicate;
import java.util.function.DoubleBinaryOperator;

import java.util.Arrays;

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

    public static double[][][] qr_gschmidt(double[][] mat) {
        int[] shp = mat_shape(mat);

        double[][] Q = new double[shp[0]][shp[1]];
        double[][] R = new double[shp[0]][shp[1]];

        double[] col0 = mat_grab_column(mat,0);
        for (int n = 0; n < shp[1]; n++) {
            //first, get original column
            double[] ocol = mat_grab_column(mat,n);
            //next, calculate ncol: ocol less the projections of ocol onto each of the previous ncols
            double[] uk_est = Arrays.copyOf(ocol, ocol.length);
            for (int n1 = 0; n1 < n; n1++) {

                double[] currentQCol = mat_grab_column(Q,n1);
                double Qcol_ukest_dp = vec_dot_product(uk_est,currentQCol);

                double uk_est_magnitude = vec_dot_product(uk_est,uk_est);
                double Q_col_magnitude = vec_dot_product(currentQCol,currentQCol);

                double[] projection = vec_scalar_mult(Qcol_ukest_dp/uk_est_magnitude,uk_est);

                uk_est = vec_sub(uk_est, projection);
                double uk_est_norm = vec_norm(uk_est);

                R[n1][n] = Qcol_ukest_dp;

            }

            double uk_est_norm = vec_norm(uk_est);
            R[n][n] = uk_est_norm;

            double[] evec = vec_scalar_mult(1.0/uk_est_norm, uk_est);
            for (int m = 0; m < shp[0]; m++) {
                Q[m][n] = evec[m];
            }
        }
        mat_print(mat, "input");
        mat_print(Q, "Q");
        mat_print(R, "R");

        double[][][] out = new double[2][shp[0]][shp[1]];
        out[0] = Q;
        out[1] = R;
        return out;
    }
}
