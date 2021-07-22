package com.company;

/*
This is a dummy class for performing several mathematical operations that will be needed
in the other classes/methods
 */

public class Mathematics {

    private static final double PI = 3.1415926535897932384626;
    Mathematics() {

    }

    // Calculate the faculty of a given integer
    public int fact (int input) {
        int fact = 1; // the result
        if (input <= 1) {
            fact = 1;
        } else {
            for (int i = 1; i <= input; i++) {
                fact = fact *i;
            }
        }
        return fact;
    }

    // Calculate the double faculty of a given integer (n!!)
    public int factDouble (int input) {

        int factDouble = 0;
        // distinct between even and odd numbers! (even numbers first)
        if (input % 2 == 0) {
            int halfInput = input/2;
            factDouble = (int) (Math.pow(2, halfInput) * this.fact(halfInput));
        } else {
            int halfInput = (input+1)/2;
            factDouble = (int) (this.fact(2 * halfInput) / (Math.pow(2, halfInput) * this.fact(halfInput)));
        }
        return factDouble;
    }



    // Calculate the center of the product Gaussian resulting from two single Gaussians

    public double[] gaussProduct (double expA, double[] posA, double expB, double[] posB) {

        double [] productPos = new double[3];

        for (int i = 0; i < 3; i++) {
            productPos[i] = (expA * posA[i] + expB * posB[i]) / (expA + expB);
        }
        return productPos;

    }

    // Calculate the difference vector between two given points in space (must have the same dimension!)

    public double[] distVector (double[] vectorA, double[] vectorB) {
        int dim = vectorA.length;
        if (dim != vectorB.length) {
            System.err.println("Error in function distVector! The given vectors do not have the same dimension!");
        }
        double[] distVector = new double[dim];
        for (int i = 0; i < dim; i++) {
            distVector[i] = vectorA[i] - vectorB[i];
        }
        return distVector;
    }

    // Calculate the length of an arbitrary vector

    public double vectorLength(double[] vectorA) {
        double vectorLength = 0;
        for (double v : vectorA) {
            vectorLength = vectorLength + Math.pow(v, 2);
        }
        return Math.sqrt(vectorLength);
    }

    // Normalize a vector

    public double[] normVector(double[] vectorA) {
        int dim = vectorA.length;

        double lengthA = this.vectorLength(vectorA);

        for (int i = 0; i < dim; i++) {
            vectorA[i] = vectorA[i]/lengthA;
        }
        return vectorA;

    }

    // Calculate the dot product (inner product) of two vectors of the same length

    public double dotProduct(double[] vectorA, double[] vectorB) {
        int dimA = vectorA.length;
        int dimB = vectorB.length;

        double product = 0;
        if (dimA == dimB) {
            for (int i = 0; i < dimA; i++) {
                product = product + vectorA[i] * vectorB[i];
            }
        } else {
            System.err.println("ERROR! Both vectors must have the same length!");
        }
        return product;
    }

    // Calculate the dyadic tensor of two given vectors of the same dimension

    public double[][] dyadic(double[] vectorA, double[] vectorB) {
        int dimA = vectorA.length;
        int dimB = vectorB.length;

        double[][] product;

        // Check if the lengths of both vectors are the same
        if (dimA == dimB) {
            product = new double[dimA][dimA];
            for (int i = 0; i < dimA; i++) {
                for (int j = 0; j < dimA; j++) {
                    product[i][j] = vectorA[i]*vectorB[j];
                }
            }
        } else {
            product = new double[1][1];
            System.err.println("ERROR! Both vectors must have the same length!");
        }
        return product;
    }


    // Calculate the error function as a chain fractional expansion
    // The deviation from the exact function is always smaller than 1E-13!

    public double errorFunction(double x) {


        double erf = 0;
        // First case: x is small (smaller than 1.25): calculate a Taylor series:
        if (x < 0.5) {
            // For 13 series terms: precalculate the faculties!
            double[] pars = new double[]{1.0, 3.0, 10.0, 42.0, 216.0, 1320.0, 9360.0, 75600.0, 685440.0, 6894720.0,
                    76204800.0, 918086400.0, 11975040000.0, 168129561600.0, 2528170444800.0};
            for (int i = 1; i < 14; i++) {
                // distinguish between even and odd terms for the sign
                erf = erf - Math.pow(-1, i + 1) * Math.pow(x, 2 * i - 1)/pars[i - 1];
            }
            // Finally apply the prefactor of the whole Taylor series
            erf = erf * -2/Math.sqrt(PI);
        } else if (x < 3.5) {
            // values of the factorials
            double[] factorials = new double[]{1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0,
                    39916800.0, 479001600.0, 6227020800.0};

            // first, define the parameters that depend on the actual reference point of the series expansion
            double center;
            double erfValue;
            if (x < 1) {
                center = 0.75;
                erfValue = 0.71115563365351508;
            } else if (x < 1.5) {
                center = 1.25;
                erfValue = 0.92290012825645829;
            } else if (x < 2) {
                center = 1.75;
                erfValue = 0.98667167121918242;
            } else if (x < 2.5) {
                center = 2.25;
                erfValue = 0.99853728341331882;
            } else if (x < 3) {
                center = 2.75;
                erfValue = 0.99989937807788032;
            } else {
                center = 3.25;
                erfValue = 0.99999569722053627;
            }
            // Calculate the generalized Taylor series for the error function, depending on the point where
            // the series is located:
            // erf_a(x) = erf(a) - 2/sqrt(pi) exp(-a^2) * sum_{n_1}^10 (-1)^n H_{n-1}(a)/n! (x-a)^n

            double old1=0;
            double old2=0;
            double hermiteValue;
            for (int i = 1; i < 14; i++) {
                if (i < 3) {
                    hermiteValue = hermite(i - 1, center);
                } else {
                    hermiteValue = hermite(i-1,center,old1,old2);
                }
                erf = erf - Math.pow(-1, i) * hermiteValue / factorials[i-1] * Math.pow(x - center, i);
                old2=old1;
                old1=hermiteValue;
            }
            erf = erf * 2 / Math.sqrt(PI) * Math.exp(-center * center) + erfValue;

        } else {
            // Second case: x is large (larger than 1.25): calculate continued fraction expansion:
            int length = 15; // length of the expansion
            // Iterate from the bottom of the fraction expansion to the top: divide through the last term of the
            // series in each iteration...
            erf = x + length / 2.0;
            for (int i = 1; i < length; i++) {
                if (Math.floorMod(length-i,2) == 0) {
                    erf = 2 * x + (length - i) / erf;
                } else {
                    erf = x + (length - i) / erf;
                }
            }
            erf = 1 - 1 / Math.sqrt(PI) * Math.exp(-Math.pow(x, 2))  / erf;
        }

        return erf;

    }

    // Calculate the value of the n'th Hermite Polynomial for an arbitrary double number x
    // via the recursional definition of them (H_n = 2*x*H_{n-1}(x)-2(n-1)*H_{n-2}(x)

    public double hermite(int n, double x) {
        double value = 0;
        // Only calculate the first two functions explicitly, define all others by recursion!
        if (n == 0) {
            value = 1;
        } else if (n == 1) {
            value = 2 * x;
        } else {
            value = 2 * x * hermite(n-1,x) - 2 * (n - 1) * hermite(n-2,x);
        }
        return value;
    }
    // Especially for loop structures, where Hermite functions of ascending order are called:
    // instead of recursive calls: directly give the values of H_{n-1} and H_{n-2} as input parameters!
    // n : order of the polynomial to calculate
    // x : the x value to calculate
    // first: value of hermite (n-1,x)
    // second: value of hermite (n-2,x)

    public double hermite(int n, double x, double first, double second) {
        double value = 0;
        // Only calculate the first two functions explicitly, define all others by recursion!
        if (n < 2) {
            System.err.println("This version of hermite can only calculate functions for n >=2!");
        } else {
            value = 2 * x * first - 2 * (n - 1) * second;
        }
        return value;
    }





}
