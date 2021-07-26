package com.company;

import java.util.Arrays;
import org.apache.commons.math3.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/*
This utility class contains a number of functions and algorithms needed for the effective handling
of matrices, which is quite cumbersome in original Java without usage of some external packages.
 */
public class Matrix {

    Mathematics mathUtil;
    // The dummy constructor
    Matrix() {
        this.mathUtil = new Mathematics();
    }


    // Print out a matrix of arbitrary dimension in an acceptable nice way

    void print(double[][] matrixA) {
        // first determine the dimensions of the involved matrix
        int rowA = matrixA.length;
        int colA = matrixA[0].length;

        System.out.println("The printed matrix has the dimension " + rowA + "x" + colA + " and the entries:");
        for (int i = 0; i < rowA; i++) {
            Arrays.stream(matrixA[i]).forEach(e -> System.out.printf( "%f   ",e ));
            System.out.print(" \n");
        }
    }

    // Print out a vector in an acceptable nice way

    void print (double[] vectorA) {
        // determine the dimension of the vector
        int dimA = vectorA.length;
        System.out.println("The printed vector has the dimension " + dimA + " and the elements:");
        for (int i = 0; i < dimA; i++) {
            System.out.println(vectorA[i]);
        }
    }

    // Calculate the transpose of a given matrix

    double[][] transpose(double[][] matrixA) {
        // first determine the dimensions of the involved matrix
        int rowA = matrixA.length;
        int colA = matrixA[0].length;

        double[][] matrixB = new double[colA][rowA];

        // columns of A will be the rows of B...
        for (int i = 0; i < colA; i++) {
            for (int j = 0; j < rowA; j++) {
                matrixB[i][j] = matrixA[j][i];
            }
        }
        return matrixB;

    }

    // Calculate the sum of two matrices

    double [][] sum(double[][] matrixA, double[][] matrixB) {
        // first determine the dimensions of the involved matrices
        int rowA = matrixA.length;
        int colA = matrixA[0].length;
        int rowB = matrixB.length;
        int colB = matrixB[0].length;

        double[][] matrixC;
        // only perform the summation if the dimensions are suitable (rowA=rowB, colA=colB)
        if (rowA == rowB && colA == colB) {
            matrixC = new double[rowA][colA];
            for (int i = 0; i < rowA; i++) {
                for (int j = 0; j < rowB; j++) {
                    matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
                }
            }
        } else {
            matrixC = new double[1][1];
            System.err.println("ERROR! The dimensions of the matrices are not equal!");
        }
        return matrixC;
    }

    // Calculate the product of two matrices (with two dimensions each)
    double[][] product(double[][] matrixA, double[][] matrixB) {
        // first determine the dimensions of the involved matrices
        int rowA = matrixA.length;
        int colA = matrixA[0].length;
        int rowB = matrixB.length;
        int colB = matrixB[0].length;

        double[][] matrixC;
        // only perform the product if the dimensions are suitable (colA = rowB)
        if (colA == rowB) {
            matrixC = new double[rowA][colB];
            for (int i = 0; i < rowA; i++) {
                for (int j = 0; j < colB; j++) {
                    for (int k = 0; k < colA; k++) {
                        matrixC[i][j] = matrixC[i][j] + matrixA[i][k] * matrixB[k][j];
                    }
                }
            }
        } else {
            matrixC = new double[1][1];
            System.err.println("ERROR! The number of columns of matrix1 must be equal to the number of " +
                    "rows of matrix2!");
        }
        return matrixC;

    }

    // calculate the product of a matrix and a vector, resulting in a vector
    double[] product(double[][] matrixA, double[] vectorB) {
        // first determine the dimensions of the matrix and the vector
        int rowA = matrixA.length;
        int colA = matrixA[0].length;
        int dimB = vectorB.length;

        double[] vectorC;
        // only perform the product if the dimensions are suitable (colA = dimB)
        if (colA == dimB) {
            vectorC = new double[dimB];
            for (int i = 0; i < rowA; i++) {
                for (int j = 0; j < dimB; j++) {
                    vectorC[i] = vectorC[i] + matrixA[i][j] * vectorB[j];
                }
            }
        } else {
            vectorC = new double[1];
            System.err.println("ERROR! The number of columns of the matrix must be equal to the number of " +
                    "rows of the vector!");
        }
        return vectorC;

    }

    // perform the QR decomposition of an arbitrary real quadratic matrix A
    // Two matrices are obtained as result: The orthogonal matrix Q and the upper triangular matrix R
    // Both will be stored in the returned 3D array: first dimension: Q, second dimension: R

    double[][][] decomposeQR(double[][] matrixA) {
        // first determine the dimensions of the given matrix
        int rowA = matrixA.length;
        int colA = matrixA[0].length;


        if (rowA != colA) {
            System.err.println("ERROR! The QR-decomposition was only implemented for square matrices!");
        }


        double[] wVec = new double[rowA];
        double[][] matrixR = new double[rowA][rowA];
        // Define a Q matrix to be the unit matrix at the beginning

        double[][] matrixQ = new double[rowA][rowA];
        for (int i = 0; i < rowA; i++) {
            matrixQ[i][i] = 1;
        }

        for (int i = 0; i < rowA; i++) {
            // first, determine the Householder vector for the current column
            // absV: the absolute value of the i'th column of the matrix (or rather of the current segment, see below)
            // its dimension is therefore reduced by one in each iteration...

            double[] vVec = new double[rowA - i];
            for (int j = 0; j < rowA-i; j++) {
                vVec[j] = matrixA[j+i][i];
            }
            double absV = this.mathUtil.vectorLength(vVec);


            // The  Householder vector is built as follows:
            // The unit vector e^1 (multiplied by the norm of the column times the signum of its first entry)
            // is substracted from the current column of the matrix.
            // In each outer iteration i, the matrix unter consideration is scaled down by the left column and
            // the upper row to consist of the quadratic segment in the lower right.
            // The same procedure as explained is then repeated for the segment; whereas the upper entries of the
            // Householder vector are filled with zeros

            for (int j = 0; j < rowA; j++) {
                if (j == i) {
                    wVec[j] = matrixA[i][i] + Math.signum(matrixA[i][i])*absV;
                } else if (j < i) {
                    wVec[j] = 0;
                } else {
                    wVec[j] = matrixA[j][i];
                }
            }

            // normalize the Householder vector
            wVec = this.mathUtil.normVector(wVec);



            // Apply the householder vector to each column of the matrix, separately!
            // This is no explicit matrix multiplication, such that the scaling is N^2 instead of N^3

            // In order to get the orthogonal Q matrix, apply the multiplications twice: once to the
            // intial matrix A and once to the Q matrix, which is simply the unit matrix at the beginning

            for (int j = 0; j < rowA; j++) {
                // calculate the scalar product of the actual column and the Householder vector
                double[] aj = new double[rowA];
                double[] qj = new double[rowA];

                for (int k = 0; k < rowA; k++) {
                    aj[k] = matrixA[k][j];
                    qj[k] = matrixQ[k][j];
                }
                double wDotA = 2 * this.mathUtil.dotProduct(wVec,aj);
                double wDotQ = 2 * this.mathUtil.dotProduct(wVec,qj);

                for (int k = 0; k < rowA; k++) {
                    matrixA[k][j] = matrixA[k][j] - wDotA * wVec[k];
                    matrixQ[k][j] = matrixQ[k][j] - wDotQ * wVec[k];
                }

            }

           ;
        }
        // The resulting matrix Q must be transposed since it is set together in the wrong direction
        // (should be Q = Q_N*Q_{N-1}*...*Q_2*Q_1 instead of Q_1*Q_2*...*Q_{N-1}*Q_N)

        matrixQ = this.transpose(matrixQ);

        // The return section: define the 3D array with the matrix R in the zeroth and the matrix Q in the first
        // outer dimension entry


        double[][][] result = new double[2][rowA][rowA];

        for (int i = 0; i < rowA; i++) {
            for (int j = 0; j < rowA; j++) {
                result[0][i][j] = matrixA[i][j];
                result[1][i][j] = matrixQ[i][j];
            }
        }

        return result;
    }

    /*
    This subroutine manages the calculation of eigenvalues and eigenvectors of a square matrix A.
    It uses the QR algorithm; both in the simplest form of repeated QR decompositions until convergence (for
    benchmark) and the more subtle form of initial conversation to Hessenberg form and subsequent QR
    transformations. The first methods requires M iterations of N^3 costs, the second requires one N^3 and M-1
    N^2 processes.
    Returned will be a matrix of dimension (N+1)xN, where the first row contains the eigenvalues and the
    respective columns contains the corresponding eigenvector.
    The correct eigenvectors will only be computed for symmetric matrices! However, since the Fock operator
    is hermitian, this requirement is unproblematic in the case of simple SCF cycles.
    For benchmark purposes, the external Apache Common Math library routines for Eigenvalue and -vector
    calculations are available as well.
     */
    double[][] eigenVV (double[][] matrixA, String method) {
        //String method = "external"; // simple or external

        double[][] matrixInput = matrixA;
        // Set the convergence criterion: the largest off diagonal element allowed when convergence is signaled.
        double crit = 1E-10;

        // Determine the dimensions of the input matrix A (only diagonal matrices allowed!)
        int rowA = matrixA.length;
        int colA = matrixA[0].length;

        if (rowA != colA) {
            System.err.println("ERROR! The QR-algorithm was only implemented for square matrices!");
        }



        // Initialize the matrix with all eigenvectors: it will be the product of all Q matrices
        // working on the unit matrix
        double[][] eigenvectors = new double[rowA][rowA];
        for (int i = 0; i < rowA; i++) {
            eigenvectors[i][i] = 1;
        }
        //this.print(eigenvectors);

        double[] eigenvalues = new double[rowA];

        if (method.equals("simple")) {


            // The largest off-diagonal value (see above)
            double residual = 0;
            int iterate = 0;
            int iterate2 = 0;
            do {
                iterate++;
                // First perform the QR decomposition of A (A_k = Q_k * R_k), then calculate the next iteration by
                // inverse multiplication: (A_{k+1} = R_k * Q_k).

                double[][][] resultQR = this.decomposeQR(matrixA);
                // Initialize the Q and R matrices resulting from the QR decomposition
                double[][] matrixQ = resultQR[1];
                double[][] matrixR = resultQR[0];


                // The matrix with the eigenvectors (product of all Q_i, where the new Q is multiplied from the right)
                //this.print(matrixQ);
                eigenvectors = this.product(eigenvectors, matrixQ);

                //this.print(eigenvectors);
                //System.exit(0);
                // The actualized matrix A
                matrixA = this.product(matrixR, matrixQ);

                //residual = 0;
                // Check the largest lower triangular off diagonal element of A_{k+1}
                //for (int i = 0; i < rowA; i++) {
                //    for (int j = 0; j < rowA; j++) {
                //        if (j < i && matrixA[i][j] > residual) {
                //            residual = matrixA[i][j];
                //        }
                //    }
                //}

                //System.out.println("Iteration: " + iterate + " residual " + residual + " crit " + crit);


                iterate2++;

                //} while (residual > crit);
            } while (iterate2 < 5000);

            // store the diagonal elements as eigenvalues
            for (int i = 0; i < rowA; i++) {
                eigenvalues[i] = matrixA[i][i];
            }

        // Benchmark purposes: Calculate the eigenvalues and eigenvectors of the matrix using the Apache
        // Commons Math Library

        } else if (method.equals("external")) {

            // Force the input matrix A to be symmetric, else, strange and useless errors are thrown -,-
            for (int i = 0; i < rowA; i++) {
                for (int j = i; j < rowA; j++) {
                    matrixA[i][j] = matrixA[j][i];
                }
            }

            RealMatrix matrixA_real = new org.apache.commons.math3.linear.Array2DRowRealMatrix(matrixA);
            double returnVal = 0;
            double[] eigenvaluesSort =
                    new org.apache.commons.math3.linear.EigenDecomposition(matrixA_real,returnVal).getRealEigenvalues();

            // Sort the eigenvalues: lowest one first
            for (int i =0; i < rowA; i++) {
                eigenvalues[rowA-i-1] = eigenvaluesSort[i];
                //eigenvalues[i] = eigenvaluesSort[i];
            }

            //System.out.println("test " + returnVal);
            double[][] eigenvectorsSort = new double[rowA][rowA];
            // The eigenvectors must be calculated once for each...
            for (int i = 0; i < rowA; i++) {


                RealVector eigenvectorI = new org.apache.commons.math3.linear.
                        EigenDecomposition(matrixA_real,returnVal).getEigenvector(i);
                for (int j = 0; j < rowA; j++) {
                    eigenvectorsSort[i][j] = eigenvectorI.getEntry(j);
                }
            }
            // Resort the entries of the eigenvector matrix...

            for (int i = 0; i < rowA; i++) {
                for (int j = 0; j < rowA; j++) {
                    //eigenvectors[rowA-j-1][i] = -eigenvectorsSort[j][i];
                    eigenvectors[i][rowA-j-1] = eigenvectorsSort[j][i];
                }
            }


        }



        double[][] result = new double[rowA+1][rowA];
        for (int i = 0; i < rowA; i++) {
            result[0][i] = eigenvalues[i];
            for (int j = 0; j < rowA; j++) {
                result[j+1][i] = eigenvectors[j][i];
            }
        }
        //System.out.println("The eigenvectors of the system are:");
        //this.print(eigenvectors);
        // Check the eigenvector matrix

        double[][] matrixTest = this.product(this.transpose(eigenvectors),this.product(matrixInput,eigenvectors));

        //System.out.println("The test eigenvalue matrix:");
        //this.print(matrixTest);

        return result;
    }


}
