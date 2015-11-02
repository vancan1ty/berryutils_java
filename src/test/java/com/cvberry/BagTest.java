package com.cvberry;

import org.junit.Test;
import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.Random;

import static com.cvberry.berryutils.BagOfAlgorithms.*;

public class BagTest 
{


    @Test
    public void testEquals()
    {
        double[][] mat1 = new double[][]{{1,2,3},{4,5,6},{7,8,9}};
        double[][] mat2 = new double[][]{{1,2,3},{4,5,6},{7,8,9}};
        double[][] mat3 = new double[][]{{9,8,7},{6,5,6},{7,8,9}};
        double[][] mat4 = new double[][]{{1.02,2,3},{4,5,6},{7,8,9}};

        assertTrue(mat_equals(mat1,mat2));
        assertTrue(mat_equals(mat1,mat2,0.00001));
        assertFalse(mat_equals(mat1,mat3));
        assertFalse(mat_equals(mat1,mat3,0.00001));
        assertFalse(mat_equals(mat1,mat4));
        assertFalse(mat_equals(mat1,mat4,0.00001));
        assertTrue(mat_equals(mat1,mat4,0.025));
                              
    }

    @Test
    public void testAdd()
    {
        double[][] mat1 = new double[][]{{1,2,3},{4,5,6},{7,8,9}};
        double[][] mat2 = new double[][]{{9,8,7},{6,5,4},{3,2,1}};
        double[][] properResult = new double[][]
            {{10,10,10},
             {10,10,10},
             {10,10,10}};
        assertTrue(mat_equals(mat_add(mat1,mat2), properResult));
    }

    @Test
    public void testSum()
    {
        double[][] mat1 = new double[][]{{1,2,3},{4,5,6},{7,8,9}};
        boolean[][] mat2 = new boolean[][]
            {{false,true,false},
             {false,false,true},
             {true,false,true}};

        assertTrue(mat_sum(mat1)==45.0);

        assertTrue(mat_sum(mat2)==4);
    }

    @Test
    public void testMean()
    {
        double[][] mat1 = new double[][]{{1,2,3},{4,5,6},{7,8,9}};

        assertTrue(mat_mean(mat1)==5.0);
    }



    @Test
    public void test_flipUD()
    {
        double[][] mat1 = new double[][]{{1,2,3},
                                         {4,5,6},
                                         {7,8,9}};
        double[][] tflipped = new double[][]{{7,8,9},
                                            {4,5,6},
                                            {1,2,3}};

        double[][] mat2 = new double[][]{{1,2,3},
                                         {1,1,1},
                                         {4,5,6},
                                         {7,8,9}};
        double[][] tflipped2 = new double[][]{{7,8,9},
                                            {4,5,6},
                                            {1,1,1},
                                            {1,2,3}};

        assertTrue(mat_equals(mat_flipud(mat1),tflipped));
        assertTrue(mat_equals(mat_flipud(mat2),tflipped2));
    }

    @Test
    public void test_flipLR()
    {
        double[][] mat1 = new double[][]{{1,2,3},
                                         {4,5,6},
                                         {7,8,9}};
        double[][] tflipped = new double[][]{{3,2,1},
                                             {6,5,4},
                                             {9,8,7}};
        assertTrue(mat_equals(mat_fliplr(mat1),tflipped));
    }

    @Test
    public void test_transpose()
    {
        double[][] mat1 = new double[][]{{1,2,3},
                                         {4,5,6},
                                         {7,8,9}};
        double[][] tflipped = new double[][]{{1,4,7},
                                             {2,5,8},
                                             {3,6,9}};
        assertTrue(mat_equals(mat_transpose(mat1),tflipped));
    }

    @Test
    public void test_OR()
    {
        boolean[][] mat1 = new boolean[][]
            {{false,true,false},
             {false,false,true},
             {true,false,false}};

        boolean[][] mat2 = new boolean[][]
            {{false,false,true},
             {true,true,true},
             {true,true,false}};

        boolean[][] tresult = new boolean[][]
            {{false,true,true},
             {true,true,true},
             {true,true,false}};
        assertTrue(mat_equals(tresult,
                              mat_OR(mat1,mat2)));
    }

    @Test
    public void test_AND()
    {
        boolean[][] mat1 = new boolean[][]
            {{false,true,false},
             {false,false,true},
             {true,false,false}};

        boolean[][] mat2 = new boolean[][]
            {{false,false,true},
             {true,true,true},
             {true,true,false}};

        boolean[][] tresult = new boolean[][]
            {{false,false,false},
             {false,false,true},
             {true,false,false}};
        assertTrue(mat_equals(tresult,
                              mat_AND(mat1,mat2)));
    }

    @Test
    public void test_dp() {
        double[] vec1 = new double[]{1, 0};
        double[] vec2 = new double[]{0, 1};
        double[] vec3 = new double[]{1/Math.sqrt(2), 1/Math.sqrt(2)};

        assertTrue(vec_dot_product(vec1,vec2)==0.0);
        assertTrue(vec_dot_product(vec1,vec1)==1.0);
        assertTrue(vec_dot_product(vec1,vec3)==1/Math.sqrt(2));
        assertEquals(1.0,vec_dot_product(vec3,vec3),0.00001);
    }

    @Test
    public void test_vec_proj() {
        double[] vec1 = new double[]{1, 0};
        double[] vec2 = new double[]{0, 1};
        double[] vec3 = new double[]{Math.sqrt(2), Math.sqrt(2)};
        double[] vec4 = new double[]{Math.sqrt(2), 0};
        double[] zeros = new double[]{0, 0};

        assertTrue(Arrays.equals(vec_proj_a_onto_e(vec1,vec2), zeros));
        assertTrue(Arrays.equals(vec_proj_a_onto_e(vec1,vec1), vec1));
        assertTrue(Arrays.equals(vec_proj_a_onto_e(vec3,vec1), vec4));

    }

    @Test
    public void test_qr1() {
        double[][] mat = new double[][]{{12, -51, 4},
                                        {6, 167, -68},
                                        {-4, 24, -41}};

        double[][][] results = qr_gschmidt(mat);

        double[][] Q = results[0];
        double[][] R = results[1];

        double[][] recomp = mat_mult(Q,R);
        //mat_print(recomp,"recomposed");
        assertTrue(mat_equals(recomp,mat,0.0001));
    }

    @Test
    public void test_qr2() {
        Random random = new Random();
        long seed = random.nextLong();
        System.out.format("qr2 -- using seed: %d\n", seed);
        for (int i = 0; i < 25; i++) {
            double[][] mat = mat_rand(8,8);
            double[][][] factorization = qr_gschmidt(mat);
            double[][] Q = factorization[0];
            double[][] R = factorization[1];
            double[][] recomp = mat_mult(Q,R);
            assertTrue(mat_equals(recomp,mat,0.0001));
        }
    }

    @Test
    public void test_qr_underdetermined() {
        double[][] mat1 = new double[][]{{1, 2, 0},
                                        {1, 3, 0},
                                        {0, 0, 0}};

        double[][] mat2 = new double[][]{{1, 2, 3},
                                         {1, 3, 5},
                                         {0, 0, 0}};

        double[][] mat3 = new double[][]{{1, 0, 2},
                                         {1, 0, 3},
                                         {0, 0, 0}};


        testQRFactorization(mat1);
        testQRFactorization(mat2);
        testQRFactorization(mat3);
    }

    @Test
    public void test_qr_skinny() {
        double[][] mat1 = new double[][]{{1, 2},
                                        {1, 3},
                                        {0, 0}};

        double[][] mat2 = new double[][]{{0, 1},
                                         {0, 0},
                                         {0, 1}};

        double[][] mat3 = new double[][]{{1, 1, 1},
                                         {1, 1, 0},
                                         {1, 0, 2},
                                         {1, 0, 0}
        };


        testQRFactorization(mat1);
        testQRFactorization(mat2);
        testQRFactorization(mat3);
    }

    @Test
    public void test_qr_wide() {
        double[][] mat1 = new double[][]{{1,2,1},
                                        {1,3,5}};

        double[][] mat2 = new double[][]{{0, 1,0},
                                         {0, 0,1}};

        double[][] mat3 = mat_transpose(new double[][]{{1, 1, 1},
                                         {1, 1, 0},
                                         {1, 0, 2},
                                         {1, 0, 0}
        });


        testQRFactorization(mat1);
        testQRFactorization(mat2);
        testQRFactorization(mat3);
    }


    public static void testQRFactorization(double[][] matrix) {
        double[][][] factorization = qr_gschmidt(matrix);
        double[][] Q = factorization[0];
        double[][] R = factorization[1];
        double[][] recomp = mat_mult(Q,R);
        // mat_print(matrix, "input");
        // mat_print(Q, "Q");
        // mat_print(R, "R");
        // mat_print(recomp, "recomp");
        // System.out.println();

        assertTrue(mat_equals(recomp,matrix,0.0001));
    }


    @Test
    public void test_mult() {
        double[][] mat = new double[][]{{2, 1, 12},
                                        {3, 1, -5},
                                        {-2, 5, -4}};

        double[][] omat = new double[][]{{3, 1, 12},
                                        {3, 1, -5},
                                        {-2, 5, -4}};


        double[][] matsquared = new double[][]{
                                        {-17, 63, -29},
                                        {19, -21, 51},
                                        {19, -17, -33}};

        assertTrue(mat_equals(mat_mult(mat,mat),matsquared,0.0001));
        assertFalse(mat_equals(mat_mult(omat,mat),matsquared,0.0001));
    }

    @Test
    public void test_mult2() {
        double[][] mat = new double[][]{{2, 1, 12},
                                        {3, 1, -5}};

        double[][] tresult = new double[][]
            {{149, -53},
             {-53, 35}};

        assertTrue(mat_equals(mat_mult(mat,mat_transpose(mat)),tresult,0.0001));
    }

    @Test
    public void test_qrTargeted() {
        double[][] mat = new double[][]{{2, 1, 2},
                                        {3, 1, 2},
                                        {-2, 5, 10}};

        testQRFactorization(mat);
    }

    @Test
    public void test_backSubstitution() {
        double[][] R1 = new double[][]
            {{1, 5, 3},
             {0, 2, 4},
             {0, 0, 9}};

        double[] y1 = new double[]{1,4,9};
        double[] ans1 = new double[]{-2,0,1};

        double[] out = mat_backSubstUpperTriangularMatrix(R1,y1);
        //vec_print(out,"result");

        assertTrue(vec_equals(out,ans1,0.0000001));

        double[][] R2 = new double[][]{
            {0.817453, 0.621355, 0.673689},
            {0.000000, 0.672652, -0.233282},
            {0.000000, 0.000000, 0.571856}};

    }

    @Test
    public void test_linearSolveSquare() {
        Random random = new Random();
        long seed = random.nextLong();
        System.out.format("linearSolveSquare -- using seed: %d\n", seed);
        for (int i = 0; i < 25; i++) {
            double[][] mat = mat_rand(8,8);
            double[] vec = vec_rand(mat_shape(mat)[0]);
            checkLinSolveSolution(mat,vec);
        }
    }

    @Test
    public void test_linearSolveWide() {
        Random random = new Random();
        long seed = random.nextLong();  //3580959054035785765l;
        System.out.format("linearSolveWide -- using seed: %d\n", seed);
        for (int i = 0; i < 25; i++) {
            int nWidth = random.nextInt(5)+5;
            double[][] mat = mat_rand(4,nWidth);
            double[] vec = vec_rand(mat_shape(mat)[0]);
            checkLinSolveSolution(mat,vec);
        }
    }


    //for underdetermined or rull-rank systems.
    public void checkLinSolveSolution(double[][] A, double[] b) {
        //mat_print(A,"A");
        double[][][] qr = qr_gschmidt(A);
        double[][] Q = qr[0];
        double[][] R = qr[1];
        double[] solution = _mat_solveSystemNoLstSq(Q,R,b);
        // mat_print(Q,"Q");
        // mat_print(R,"R");
        // vec_print(solution,"solution");
        // vec_print(b,"b");
        double[] reconstructed = mat_vec_mult(A,solution);
        // vec_print(reconstructed,"reconstructed");
        assertTrue(vec_equals(b,reconstructed,0.000001));
    }

    @Test
    public void test_mat_vec_mult() {
        double[] a = new double[]{1,2,3};
        double[][] mat = new double[][]{{2, 1, 2},
                                        {3, 1, 2},
                                        {-2, 5, 10}};
        double[] correct = new double[]{10,11,38};

        double[] out = mat_vec_mult(mat,a);
        assertTrue(vec_equals(correct,out,0.000001));
    }

}
