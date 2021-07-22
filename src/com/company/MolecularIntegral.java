package com.company;

import java.text.DecimalFormat;
import java.util.Arrays;

/*
This class provides all tools necessary to calculate the different arrays containing
molecular integrals:
- overlap integral
- kinetic energy integral
- nuclear attraction integral
- two electron repulsion integral
All these integrals are attributes of this class and are calculated by invoking the
respective getter and setter methods
 */
public class MolecularIntegral {
    private static final double PI = 3.1415926535897932384626;
    private static final double ROOTPI = 1.7724538509055159;
    private static final double PI25 = 17.493419544621673;

    // Number of basis functions in the molecule, i.e. the dimension of all integral arrays
    int basisSize;
    // Number of atoms in the molecule, needed for the nuclear attraction integrals
    int nAtoms;
    // Coordinates of the atoms in the system
    double[][] atomCoordinates;
    // Atomic numbers of all atoms in the system
    int[] atomicNumbers;
    // Array with BasisFunction objects containing all set up basis functions in the system
    BasisFunction[] fullBasis;
    // The overlap integrals between pairs of basis functions
    double[][] overlapIntegral;
    // The kinetic energy integrals
    double[][] kineticIntegral;
    // The nuclear attraction integrals
    double[][] nuclearIntegral;
    // The number of individual repulsion integrals
    int repulsionSize;
    // The two electron repulsion integrals, in a condensed 1D storage form
    double[] repulsionIntegral;
    // Local object of the class Mathematics for some useful functions etc.
    Mathematics mathUtil;
    // Local object of the class Boys for evaluation of the Boys function
    Boys boysUnit;


    // The default constructor for initialization of all arrays
    MolecularIntegral(int size, int nAtoms) {
        // The order parameters: number of basis functions and atoms in the system
        this.basisSize = size;
        this.nAtoms = nAtoms;
        if (size < 1) {
            System.err.println("The basis set seems to be corrupted! No basis set dimension found...");
        } else if (nAtoms < 1) {
            System.err.println("The system seems to include to atoms! No atom number found...");
        }
        else {
            this.atomCoordinates = new double[nAtoms][3];
            this.atomicNumbers = new int[nAtoms];
            this.fullBasis = new BasisFunction[size];
            this.overlapIntegral = new double[size][size];
            this.kineticIntegral = new double[size][size];
            this.nuclearIntegral = new double[size][size];
            // calculate the highest possible index in the 1D storage system for nuclear repulsion integrals
            int index = (size-1)*((size-1)+1)/2+(size-1);
            this.repulsionSize = index*(index+1)/2+index+1;

            this.repulsionIntegral = new double[this.repulsionSize];

        }
        this.mathUtil = new Mathematics();
        this.boysUnit = new Boys();
    }
    // Setter method for the array of atomic coordinates

    public void setAtomCoordinates(double[][] atomCoordinates) {
        this.atomCoordinates = atomCoordinates;
    }

    // Setter method for the array of atomic numbers/core charges

    public void setAtomicNumbers(int[] atomicNumbers) {
        this.atomicNumbers = atomicNumbers;
    }


    // Setter method for the array of basis functions

    public void setFullBasis(BasisFunction[] fullBasis) {
        this.fullBasis = fullBasis;
    }

    // The calculation of the overlap matrix, designed to be a setter method


    public void setOverlapIntegral() {
        // The array is already initialized to be filled with zeros..

        // Outer loop: over all (contracted) basis functions
        for (int i = 0; i < this.basisSize; i++) {
            for (int j = i; j < this.basisSize; j++) {
                BasisFunction basisI = fullBasis[i];
                BasisFunction basisJ = fullBasis[j];

                // Inner loop: over primitive Gaussians of the current pair of basis functions
                for (int k = 0; k < basisI.contraction; k++) {
                    for (int l = 0; l < basisJ.contraction; l++) {

                        this.overlapIntegral[i][j] = this.overlapIntegral[i][j] + basisI.coeffs[k] * basisJ.coeffs[l] *
                                basisI.norms[k] * basisJ.norms[l] * overlap(basisI.exps[k], basisI.angularNum,
                                basisI.position, basisJ.exps[l], basisJ.angularNum, basisJ.position);
                    }
                }

                this.overlapIntegral[j][i] = this.overlapIntegral[i][j];
            }

        }

        // Debug: print out the overlap integral
        System.out.println("The overlap matrix has the dimension " + this.basisSize + " and the entries:");
        for (int i = 0; i < this.basisSize; i++) {
            //DecimalFormat df = new DecimalFormat("00.######E0");
            //Arrays.stream(overlapIntegral[i]).forEach(e -> System.out.print(df.format(e) + "   " ));
            Arrays.stream(overlapIntegral[i]).forEach(e -> System.out.printf( "%e   ",e ));
            System.out.print(" \n");
           // System.out.println(Arrays.toString(overlapIntegral[i]));
        }

    }


    // The calculation of the kinetic energy matrix, designed to be a setter method

    public void setKineticIntegral() {
        // The array is already initialized to be filled with zeros..

        // Outer loop: over all (contracted) basis functions
        for (int i = 0; i < this.basisSize; i++) {
            for (int j = 0; j <= i ; j++) {
                BasisFunction basisI = fullBasis[i];
                BasisFunction basisJ = fullBasis[j];

                // Inner loop: over primitive Gaussians of the current pair of basis functions
                for (int k = 0; k < basisI.contraction; k++) {
                    for (int l = 0; l < basisJ.contraction; l++) {

                        this.kineticIntegral[i][j] = this.kineticIntegral[i][j] + basisI.coeffs[k] * basisJ.coeffs[l] *
                                basisI.norms[k] * basisJ.norms[l] * kinetic(basisI.exps[k], basisI.angularNum,
                                basisI.position, basisJ.exps[l], basisJ.angularNum, basisJ.position);
                    }
                }

                // mirror the entries of the matrix along the diagonal

                this.kineticIntegral[j][i] = this.kineticIntegral[i][j];
            }

        }

        // Debug: print out the kinetic energy integral matrix
        System.out.println("The kinetic energy matrix has the dimension " + this.basisSize + " and the entries:");
        for (int i = 0; i < this.basisSize; i++) {
            Arrays.stream(kineticIntegral[i]).forEach(e -> System.out.printf( "%e   ",e ));
            System.out.print(" \n");
        }

    }


    // The calculation of the nuclear attraction energy matrix, designed to be a setter method

    public void setNuclearIntegral() {
        // The array is already initialized to be filled with zeros..

        // Outer loop: over all (contracted) basis functions
        for (int i = 0; i < this.basisSize; i++) {
            for (int j = i; j < this.basisSize ; j++) {

                BasisFunction basisI = fullBasis[i];
                BasisFunction basisJ = fullBasis[j];

                // Middle loop: over all nuclei with charges

                for (int n = 0; n < this.nAtoms; n++) {

                    double intAct = 0;

                    // Inner loop: over primitive Gaussians of the current pair of basis functions

                    for (int k = 0; k < basisI.contraction; k++) {
                        for (int l = 0; l < basisJ.contraction; l++) {

                            intAct = intAct + basisI.coeffs[k] * basisJ.coeffs[l] * basisI.norms[k] *
                                    basisJ.norms[l] * attract(basisI.exps[k], basisI.angularNum, basisI.position,
                                    basisJ.exps[l], basisJ.angularNum, basisJ.position, this.atomCoordinates[n]);

                         }
                    }



                    // Multiply actual Coulomb integral with charge of actual atom's core!

                    //System.out.println("int act " + intAct);
                    this.nuclearIntegral[i][j] = this.nuclearIntegral[i][j] - intAct * this.atomicNumbers[n];

                    //System.out.println("Value of the integral: " + this.nuclearIntegral[i][j]);

                }

                // mirror the entries of the matrix along the diagonal

                this.nuclearIntegral[j][i] = this.nuclearIntegral[i][j];
            }

        }

        // Debug: print out the kinetic energy integral matrix
        System.out.println("The nuclear attraction matrix has the dimension " + this.basisSize + " and the entries:");
        for (int i = 0; i < this.basisSize; i++) {
            Arrays.stream(nuclearIntegral[i]).forEach(e -> System.out.printf( "%e   ",e ));
            System.out.print(" \n");
        }

    }

    // The calculation of the two-electron repulsion integrals, designed as a setter method
    // In order to significantly reduce the calculation effort, the integrals will be stored in a
    // 1D array (instead of 4D), each combination of 8 identical integrals (due to symmetry) will only
    // be calculated and stored once
    public void setRepulsionIntegral() {
        // Set a dummy value for subsequent check if an entry has already been calculated

        double dummy=44444.44444;

        for (int i = 0; i < this.repulsionSize; i++) {
            this.repulsionIntegral[i] = dummy;
        }

        // Now start the big loop over all combinations of basis functions in the system

        int ij = 0;
        int kl = 0;
        int ijkl = 0;

        for (int i = 1; i <= this.basisSize; i++) {
            for (int j = 1; j <= this.basisSize; j++) {
                for (int k = 1; k <= this.basisSize; k++) {
                    for (int l = 1; l <= this.basisSize; l++) {
                        // Perform the (partial) mapping from 4D to 1D indices
                        if (i >= j) {
                            ij = (i - 1) * ((i - 1) + 1) / 2 + (j - 1);
                        } else {
                            ij = (j - 1) * ((j - 1) + 1) / 2 + (i - 1);
                        }

                        if (k >= l) {
                            kl = (k - 1) * ((k - 1) + 1) / 2 + (l - 1);
                        } else {
                            kl = (l - 1) * ((l - 1) + 1) / 2 + (k - 1);
                        }
                        // Now calculate the actual 1D array index

                        if (ij >= kl) {
                            ijkl = ij * (ij + 1) / 2 + kl;
                        } else {
                            ijkl = kl * (kl + 1) / 2 + ij;
                        }

                        // If the current matrix element has not been calculated so far (i.e., the default value
                        // is still existent), calculate it! Else, another identical permutation was already
                        // mapped to the current entry

                        if ((this.repulsionIntegral[ijkl] - 0.0001 < dummy) && (this.repulsionIntegral[ijkl] +
                                0.0001 > dummy)) {

                            // Define a local variable for incrementation
                            double eriAct = 0.0;
                            // Now take the actual basis function objects for the current combination

                            BasisFunction basisI = fullBasis[i - 1];
                            BasisFunction basisJ = fullBasis[j - 1];
                            BasisFunction basisK = fullBasis[k - 1];
                            BasisFunction basisL = fullBasis[l - 1];

                            // Start the inner loop over primitive Gaussians of all four current basis functions

                            for (int ia = 0; ia < basisI.contraction; ia++) {
                                for (int ib = 0; ib < basisJ.contraction; ib++) {
                                    for (int ic = 0; ic < basisK.contraction; ic++) {
                                        for (int id = 0; id < basisL.contraction; id++) {
                                            eriAct = eriAct + basisI.coeffs[ia] * basisJ.coeffs[ib] *
                                                    basisK.coeffs[ic] * basisL.coeffs[id] * basisI.norms[ia] *
                                                    basisJ.norms[ib] * basisK.norms[ic] * basisL.norms[id] *
                                                    repulsion(basisI.exps[ia], basisI.angularNum, basisI.position,
                                                            basisJ.exps[ib], basisJ.angularNum, basisJ.position,
                                                            basisK.exps[ic], basisK.angularNum, basisK.position,
                                                            basisL.exps[id], basisL.angularNum, basisL.position);
                                        }
                                    }
                                }
                            }
                            this.repulsionIntegral[ijkl] = eriAct;
                            //System.out.println("i " + i + " j " + j + " k " + k + " l " + l + "eri" + eriAct);
                        }
                    }
                }
            }
        }

    }


    /*
    This function calculates the overlap integral between two single (primitive) 3D Gaussian functions
    of arbitrary angular momentum.
    This integral is due to the exponential properties simply the product of three one-dimensional
    integrals between one-dimensional Gaussians (for x, y and z dimensions, each)
     */

    public double overlap(double expA, int[] angularA, double[] posA, double expB, int[] angularB, double[] posB) {

        // Determine shell angular momenta for both Gaussians
        double s1 = hermiteInt(angularA[0], angularB[0],0,posA[0]-posB[0], expA, expB);
        double s2 = hermiteInt(angularA[1], angularB[1],0,posA[1]-posB[1], expA, expB);
        double s3 = hermiteInt(angularA[2], angularB[2],0,posA[2]-posB[2], expA, expB);

        // Set all one dimensional overlaps together

        return Math.pow(PI / (expA + expB), 1.5) * s1 * s2 * s3;
    }

    /*
    This function calculates the kinetic energy integral between two single (primitive) 3D Gaussian functions
    of arbitrary angular momentum.
    The general expression for those integrals is well known as: T_ab = <G_a | -0.5*nabla^2 | G_b>
    Utilizing the properties of Hermite type Gaussian functions, the derivative of the left one can be rewritten
    in terms of (derivatives of) overlap integrals (now the 1D Gaussians can no longer be fully separated):
    T_ab = -0.5* (D_ij^2*S_kl*S_Mn + S_ij*D_kl^2*S_mn + S_ij*S_kl*D_mn^2)
    with D_ij^2=j(j-1)*S_(i,j-2)-2*beta*(2j+1)*S_ij + 4*beta^2*S_(i,j+2)
    with beta the orbital exponent of the second Gaussian G_b
     */

    public double kinetic(double expA, int[] angularA, double[] posA, double expB, int[] angularB, double[] posB) {

        // sum all shell angular momenta for the determination of the total sign

        double sign = Math.pow(-1, angularA[0] + angularA[1] + angularA[2] + angularB[0] + angularB[1] + angularB[2]);

        // the first term in the sum

        double term0 = expB * (2 * (angularB[0] + angularB[1] + angularB[2])+3) * overlap(expA, angularA, posA, expB,
                angularB, posB);
        // Apply the derivatives here by incrementing the angular momentum prefactors of the right Gaussian

        int[] angularIncL = new int[]{angularB[0] + 2, angularB[1], angularB[2]};
        int[] angularIncM = new int[]{angularB[0], angularB[1] + 2, angularB[2]};
        int[] angularIncN = new int[]{angularB[0], angularB[1], angularB[2] + 2};

        double term1 = -2 * Math.pow(expB, 2) * (overlap(expA, angularA, posA, expB, angularIncL, posB) +
                overlap(expA, angularA, posA, expB, angularIncM, posB) + overlap(expA, angularA, posA,
                expB, angularIncN, posB));

        // Apply the derivatives by decrementing the angular momentum prefactors of the right Gaussian

        int[] angularDecL = new int[]{angularB[0] - 2, angularB[1], angularB[2]};
        int[] angularDecM = new int[]{angularB[0], angularB[1] - 2, angularB[2]};
        int[] angularDecN = new int[]{angularB[0], angularB[1], angularB[2] - 2};

        double term2 = -0.5 * (angularB[0] * (angularB[0] - 1) * overlap(expA, angularA, posA, expB,
                angularDecL, posB) + angularB[1] * (angularB[1] - 1) * overlap(expA, angularA, posA, expB,
                angularDecM, posB) + angularB[2] * (angularB[2] - 1) * overlap(expA, angularA, posA, expB,
                angularDecN, posB));

        // Sum up all three parts

        return sign * (term0 + term1 + term2);
    }
    /*
    This function calculates the nuclear attraction integral between two single primitive Gaussians
    of arbitrary angular momentum and one of the nuclei in the system.
    Since the integral includes the expression 1/r_C, which prevents the whole integral to be factored
    out into the three spatial coordinates, the computation is somewhat more involved and cannot
    be done fully analytically. Instead, the Boys Function appears as an integral, that must be
    solved approximately (similar, but more general than the error function):
    F_j(r^2)=\int_0^1 t^{2j} exp(-r^2t^2) dt
    First, the center of both involved Gaussians must be calculated via the well known Gaussian
    product center rule. Then, an auxiliary Hermite Coulomb product center integral must be solved
    iteratively.
     */

    public double attract (double expA, int[] angularA, double[] posA, double expB, int[] angularB, double[] posB,
                           double[] posC) {
        // The sum of both exponential prefactors (prefactor for the product Gaussian)
        double expP = expA + expB;

        // Apply the Gaussian product rule for the position of the product Gaussian
        double[] posP = mathUtil.gaussProduct(expA,posA,expB,posB);

        // The distance vector between the actual nucleus and the product Gaussian
        double[] distPC = mathUtil.distVector(posP,posC);

        // The length of this distance vector
        double absPC = mathUtil.vectorLength(distPC);


        double attract = 0;
        for (int t = 0; t <= (angularA[0] + angularB[0] +1); t ++) {
            for (int u = 0; u <= (angularA[1] + angularB[1] + 1); u++) {
                for (int v = 0; v <= (angularA[2] + angularB[2] + 1); v++) {
                    attract = attract + hermiteInt(angularA[0], angularB[0], t,posA[0]-posB[0], expA, expB) *
                            hermiteInt(angularA[1], angularB[1], u,posA[1]-posB[1], expA, expB) *
                            hermiteInt(angularA[2], angularB[2], v,posA[2]-posB[2], expA, expB) *
                            hermiteCoulomb(t, u, v, 0, expP, distPC, absPC);

                }
            }
        }
        return attract * 2 * PI / expP;


    }

    /*
    This function calculates the electron repulsion integral between four single primitive Gaussian basis functions
    with the use of Hermite integration
     */

    public double repulsion (double expA, int[] angularA, double[] posA,double expB, int[] angularB, double[] posB,
                             double expC, int[] angularC, double[] posC,double expD, int[] angularD, double[] posD) {
        // Composite exponents for product Gaussians for a+b and c+d
        double expP = expA + expB;
        double expQ = expC + expD;
        // The composite-composite exponent...

        double alpha = expP * expQ / (expP + expQ);

        // Composite positions for product Gaussians of a+b and c+d
        double[] posP = mathUtil.gaussProduct(expA, posA, expB, posB);
        double[] posQ = mathUtil.gaussProduct(expC, posC, expD, posD);

        // Distance vector between both product Gaussian positions

        double[] distPQ = mathUtil.distVector(posP,posQ);
        // Absolute distance between the Gaussian positions

        double absPQ = mathUtil.vectorLength(distPQ);

        // Set repulsion to zero

        double rep = 0;

        for (int t = 0; t <= (angularA[0] + angularB[0] +1); t ++) {
            for (int u = 0; u <= (angularA[1] + angularB[1] + 1); u++) {
                for (int v = 0; v <= (angularA[2] + angularB[2] + 1); v++) {
                    for (int tau = 0; tau <= (angularC[0] + angularD[0] +1); tau ++) {
                        for (int nu = 0; nu <= (angularC[1] + angularD[1] +1); nu ++) {
                            for (int phi = 0; phi <= (angularC[2] + angularD[2] +1); phi ++) {
                                rep = rep +  hermiteInt(angularA[0], angularB[0], t,posA[0] - posB[0], expA, expB) *
                                        hermiteInt(angularA[1], angularB[1], u,posA[1] - posB[1], expA, expB) *
                                        hermiteInt(angularA[2], angularB[2], v,posA[2] - posB[2], expA, expB) *
                                        hermiteInt(angularC[0], angularD[0], tau,posC[0] - posD[0], expC, expD) *
                                        hermiteInt(angularC[1], angularD[1], nu,posC[1] - posD[1], expC, expD) *
                                        hermiteInt(angularC[2], angularD[2], phi,posC[2] - posD[2], expC, expD) *
                                        Math.pow(-1,tau + nu + phi) * hermiteCoulomb(t + tau, u + nu, v + phi,
                                        0, alpha, distPQ, absPQ);
                            }
                        }
                    }
                }
            }
        }
        // Normalize the resulting integral

        rep = rep * 2 * PI25 / (expP * expQ * Math.sqrt(expP + expQ));
        return rep;

    }





    /*
    The central function for evaluation of Hermite Gaussian coefficients (i.e., 1D-Gaussian integrals)
    This performs a recursive call depending on the angular momentum quantum numbers of the two
    involved 1D Gaussians (which it decrements) until both arrive at zero; then a simple exponent of the
    squared distance between the functions is taken
    The variables are defined as follows:
    - i : orbital angular momentum number of Gaussian 1
    - j : orbital angular momentum number of Gaussian 2
    - t : Number of nodes in Hermite integral (depends on the type of the integral, always zero for overlap integrals)
    - qx : Distance between the origins of Gaussian 1 und 2 (in bohr)
    - expA : orbital exponent of Gaussian 1
    - expB : orbital exponent of Gaussian 2
    */

    public double hermiteInt (int i, int j, int t, double qx, double expA, double expB) {

        double r; // The resulting value of the function

        // define two abbreviations
        double p = expA + expB;
        double q = expA * expB / p;

        // If the number of nodes is out of bounds: zero
        if ((t < 0) || (t > (i + j))) {
            r = 0;

            // If both angular momentum numbers are zero: base case
        } else if ((i == 0) && (j == 0) && (t == 0)) {
            r = Math.exp(-q * Math.pow(qx, 2));

            // If angular momentum j is already zero: decrement momentum index i
        } else if (j == 0) {
            r = 1.0 / ( 2.0* p) * hermiteInt(i - 1, j, t - 1, qx, expA, expB) - (q * qx / expA) *
                    hermiteInt(i - 1, j, t, qx, expA, expB) + (t + 1) *
                    hermiteInt(i - 1, j,t + 1, qx, expA, expB);

            // Last option: j is not zero, therefore decrement it
        } else {
            r = 1.0 / (2.0 * p) * hermiteInt(i, j - 1, t - 1, qx, expA, expB) - (q * qx / expB)*
                    hermiteInt(i, j - 1, t, qx, expA, expB) + (t + 1) *
                    hermiteInt(i, j - 1,t + 1, qx, expA, expB);
        }

        return r;
    }

    /*
    The function call for the Coulombg auxiliary Hermite integral recursive calculation
    as needed for the electron-nuclei attraction and electron-electron repulsion integrals.
    Depending on the three dimension-dependent quantum orbital angular momenta of the
    precalculated product Gaussian, the integral between it and the charge at distance distPC
    is calculated recursively, where the Boys function must be evaluated finally.
     */

    public double hermiteCoulomb (int t, int u, int v, int nPar, double expP, double[] distPC, double absPC) {

        // The resulting value of the function
        double r = 0;
        // first, calculate the second Boys function parameter t
        double xPar = expP * absPC * absPC;

        // If all orbital angular momentum quantum numbers are zero, evaluate the Boys function
        if ((t == 0) && (u == 0) && (v == 0)) {
            r = r + Math.pow(-2 * expP, nPar) * boysUnit.boysFunction(nPar, xPar);
        } else if ((t == 0) && (u == 0)) {
            if (v > 1) {
                r = r + (v - 1) * hermiteCoulomb(t, u, v - 2, nPar + 1, expP, distPC, absPC);
            }
            r = r + distPC[2] * hermiteCoulomb(t, u, v - 1, nPar + 1, expP, distPC, absPC);
        } else if (t == 0) {
            if (u > 1) {
                r = r + (u - 1) * hermiteCoulomb(t, u - 2, v, nPar + 1, expP, distPC, absPC);
            }
            r = r + distPC[1] * hermiteCoulomb(t, u - 1, v, nPar + 1, expP, distPC, absPC);
        } else {
            if (t > 1) {
                r = r + (t - 1) * hermiteCoulomb(t - 2, u, v, nPar + 1, expP, distPC, absPC);
            }
            r = r + distPC[0] * hermiteCoulomb(t - 1, u, v, nPar + 1, expP, distPC, absPC);
        }

        return r;
    }



}
