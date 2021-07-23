package com.company;



/*
This class contains everything that is needed to calculate the wave function and its energy for a given
system based on the Hartree-Fock approximation of no electron correlation.
So far, only closed-shell systems can be calculated based on the Roothaan-Hall equations. In the
future, the open-shell Pople-Nesbet equations will be implemented as well.
 */
public class HartreeFock {

    // The local instance of the molecular integral class, which contains all needed molecular
    // information and integrals for the SCF cycle.
    MolecularIntegral integrals;
    Matrix matrix;
    Mathematics mathUtil;


    // The obtained SCF energy
    double SCFEnergy;

    // The obtained MO coefficients
    double[][] MOCoeffs;


    // The default constructor: initialize the dummy objects for useful mathematical methods...
    HartreeFock() {

        this.matrix = new Matrix();
        this.mathUtil = new Mathematics();

    }

    // The setter method for setting all needed information...
    public void setIntegrals(MolecularIntegral integrals) {
        this.integrals = integrals;
    }

    // Separate method for initial calculation of the nuclear repulsion energy of the system.

    public double nuclearRepulsion() {
        double eNucRep = 0.0;
        double dist = 0.0;
        double[] rAB = new double[3];
        for (int i = 0; i < integrals.nAtoms; i ++) {
            for (int j = i+1; j < integrals.nAtoms; j++) {
                for (int k = 0; k < 3; k++) {
                    rAB[k] = integrals.atomCoordinates[i][k]-integrals.atomCoordinates[j][k];
                }
                dist = mathUtil.dotProduct(rAB,rAB);
                dist = Math.sqrt(dist);
                eNucRep = eNucRep + integrals.atomicNumbers[i]*integrals.atomicNumbers[j]/dist;
            }
        }
        return eNucRep;
    }

    // The method for applying the Hartree-Fock method to closed shell systems: The Roothaan-Hall equations will
    // solved iteratively during an SCF cycle. The HF energy and the MO coefficients of the molecular system
    // will then be obtained and returned indirectly by storing them to the respective attributes of the current
    // HF object.

    public double SCFClosedShell() {

        int dimension = this.integrals.basisSize;

        // Convergence criteria (will later be changeable by the user)
        double eCriterion = 1D-10;
        double densCriterion = 1D-10;
        int maxCycle = 200;


        // First, calculate the core Hamiltonian matrix (T + V(1e))
        double[][] matHCore = new double[dimension][dimension];

        matHCore = matrix.sum(this.integrals.kineticIntegral,this.integrals.nuclearIntegral);

        // For information, print out the obtained core Hamiltonian as initial guess to file.

        System.out.println("The core hamiltonian");
        matrix.print(matHCore);

        // Diagonalize the overlap matrix in order to get the orthogonalization matrix S

        double[][] dummyMatrix = this.matrix.eigenVV(this.integrals.overlapIntegral);


        int dimensionTest = 4;
        double[][] matrixIn = new double[][] {{1,8,2,4},{8,-5,-1,9},{2,-1,-11,7},{4,9,7,5}};
        double[][] matrixOut = this.matrix.eigenVV(matrixIn);

        System.out.println("The test eigenvectors and values");
        matrix.print(matrixOut);

        //System.exit(0);



        // Construct a matrix with the inverse square root eigenvalues of the core Hamiltonian on its diagonal

        double[][] lambdaMatrix = new double[dimension][dimension];
        for (int i = 0; i < dimension; i++) {
            lambdaMatrix[i][i] = 1/Math.sqrt(dummyMatrix[0][i]);
        }
        matrix.print(lambdaMatrix);

        // store the eigenvectors of the overlap matrix in the matrixLS array

        double[][] matrixLS = new double[dimension][dimension];
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrixLS[i][j] = dummyMatrix[i+1][j];
            }
        }

        System.out.println("The LS matrix:");
        matrix.print(matrixLS);
        // The transpose of the overlap eigenvectors
        double[][] transposeLS = matrix.transpose(matrixLS);

        // calculate the orthogonalization matrix S^{1/2}

        double[][] SHalf = matrix.product(matrixLS,matrix.product(lambdaMatrix,transposeLS));

        System.out.println("The S half matrix:");
        matrix.print(SHalf);

        // The transpose of the orthogonalization matrix S^{1/2}
        double[][] transposeSHalf = matrix.transpose(SHalf);

        System.out.println("The transpose S half matrix:");
        matrix.print(transposeSHalf);

        // Form the intial guess Fock matrix in the AO basis with the core Hamiltonian

        double[][] FZero = matrix.product(transposeSHalf,matrix.product(matHCore,SHalf));

        System.out.println("The Fock zero matrix:");
        matrix.print(FZero);

        //Diagonalize the Fock matrix
        dummyMatrix = this.matrix.eigenVV(FZero);



        // Extract the eigenvectors and eigenvalues

        double[][] FockVectorsOrtho = new double[dimension][dimension];
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                FockVectorsOrtho[i][j] = dummyMatrix[i+1][j];
            }
        }

        double[] orbitalEnergies = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            orbitalEnergies[i] = dummyMatrix[0][i];
        }




        System.out.println("The orbital energies");
        matrix.print(orbitalEnergies);
        // Get the eigenvectors in the original (non-orthogonal) AO basis

        double[][] FockVectorsNonOrtho = matrix.product(SHalf,FockVectorsOrtho);

        for (int i = 0; i < dimension; i++) {
            FockVectorsNonOrtho[i][0] = -FockVectorsNonOrtho[i][0];
            FockVectorsNonOrtho[i][2] = -FockVectorsNonOrtho[i][2];
        }
        System.out.println("The initial C matrix");
        matrix.print(FockVectorsNonOrtho);
        // Build the density matrix using the occupied MOs
        // First, calculate the number of electrons in the system from the atomic numbers

        int nElectrons = 0;
        for (int i = 0; i < integrals.nAtoms; i++) {
            nElectrons =  nElectrons + integrals.atomicNumbers[i];
        }

        System.out.println("The number of electrons is " + nElectrons);

        // The number of occupied orbitals (for closed shell systems, always half the electron number)
        int occupied = nElectrons/2;

        // Now build the density matrix

        double[][] densMat = new double[dimension][dimension];
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                for (int k = 0; k < occupied; k++) {
                    densMat[i][j] = densMat[i][j] + FockVectorsNonOrtho[i][k]*FockVectorsNonOrtho[j][k];
                }
            }
        }
        System.out.println("The density matrix");
        matrix.print(densMat);

        // Now calculate the intial HF energy

        double energyEl = 0.0;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                energyEl=energyEl+densMat[j][i]*(matHCore[j][i]+matHCore[j][i]);
            }
        }
        double eNucRep = this.nuclearRepulsion();
        double energyTot = energyEl + eNucRep;




        // NOW START THE SCF CYCLE!

        // first, store energy and density matrix of zeroth iteration for convergence test

        double[][] densOld = densMat;
        double energyTotOld = energyTot;
        int icycle = 0;
        double deltaE = 0.0;
        double rmsDens = 0.0;

        System.out.println("ITER.          E(elec)                  E(tot)            Delta(E)            RMS(D) ");
        System.out.println("  " + icycle + "  " + energyEl + "  " + energyTot + "  " + deltaE + "  " + rmsDens);

        icycle = 1;
        // The Fock matrix (nonorthogonal and orthogonal)
        double[][] FockMatrix = new double[dimension][dimension];
        double[][] FockMatrixPrime = new double[dimension][dimension];

        // Auxiliary loop indices for two-electron repulsion integrals
        int ij; int kl; int ik; int jl;
        int ijkl; int ikjl;
        double intCoulomb = 0.0; double intExchange = 0.0;

        do {
            // Compute the new Fock matrix
            for (int i = 0; i < integrals.basisSize; i++) {
                for (int j = 0; j < integrals.basisSize; j++) {
                    // First, add the one-electron integrals in the core Hamiltonian matrix
                    FockMatrix[i][j] = matHCore[i][j];
                    // Then, process the effectively stored two-electron integrals
                    for (int k = 0; k < integrals.basisSize; k++) {
                        for (int l = 0; l < integrals.basisSize; l++) {
                            // First, get the Coulomb integral from the 2 electron integral table:
                            // Calculate part compound indices
                            if (i >= j) {
                                ij = i * (i + 1) / 2 + j;
                            } else {
                                ij = j * (j + 1) / 2 + i;
                            }
                            if (k >= l) {
                                kl = k * (k + 1) / 2 + l;
                            } else {
                                kl = l * (l + 1) / 2 + k;
                            }
                            // Calculate actual 1D array index

                            if (ij == kl) {
                                ijkl = ij * (ij + 1) / 2 + kl ;
                            } else {
                                ijkl = kl * (kl + 1) / 2 + ij;
                            }
                            intCoulomb = integrals.repulsionIntegral[ijkl];

                            // Second, get the exchange integral from the 2 electron integral table:
                            // Now, lambda and nu (j and k) are exchanged!

                            if (i >= k) {
                                ik = i * (i + 1) / 2 + k;
                            } else {
                                ik = k * (k + 1) / 2 + i;
                            }

                            if (j >= l) {
                                jl = j * (j + 1) / 2 + l;
                            } else {
                                jl = l * (l + 1) / 2 + j;
                            }
                            // Calculate actual 1D array index

                            if (ik >= jl) {
                                ikjl = ik * (ik + 1) / 2 + jl;
                            } else {
                                ikjl = jl * (jl + 1) / 2 + ik;
                            }
                            intExchange = integrals.repulsionIntegral[ikjl];

                            // Add the two-electron part for this combination to the Fock matrix
                            FockMatrix[i][j] = FockMatrix[i][j] + densMat[k][l] * (2 * intCoulomb - intExchange);

                        }
                    }
                }
            }
            // Calculate the new density matrix by diagonalizing the new Fock matrix (always use the same
            // orthogonalization matrix from the overlap)

            FockMatrixPrime = matrix.product(transposeSHalf,matrix.product(FockMatrix,SHalf));

            // Diagonalize the new orthogonal Fock matrix in order to find the MO coefficients
            dummyMatrix = this.matrix.eigenVV(FockMatrixPrime);

            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    FockVectorsOrtho[i][j] = dummyMatrix[i+1][j];
                }
            }

            for (int i = 0; i < dimension; i++) {
                orbitalEnergies[i] = dummyMatrix[0][i];
            }

            // Transform the Fock eigenvectors into the original (non-orthogonal) AO basis
            // in order to get the MO coefficients of the eigenfunctions

            FockVectorsNonOrtho = matrix.product(SHalf,FockVectorsOrtho);

            // Calculate the density matrix from the MO coefficients


            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    densMat[i][j] = 0.0;
                    for (int k = 0; k < occupied; k++) {
                        densMat[i][j] = densMat[i][j] + FockVectorsNonOrtho[i][k]*FockVectorsNonOrtho[j][k];
                    }
                }
            }

            // Calculate the new SCF energy (use the nonorthogonal Fock matrix!)
            energyEl = 0.0;

            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < dimension; j++) {
                    energyEl=energyEl+densMat[j][i]*(matHCore[j][i]+FockMatrix[j][i]);
                }
            }
            energyTot = energyEl + eNucRep;

            // Check for convergence of the actual SCF cycle!
            // 1. calculate energy difference and compare to reference criterion
            // 2. calculate different of density matrices and compare to reference criterion

            deltaE = energyTot - energyTotOld;

            rmsDens = 0.0;
            for (int i = 0; i < integrals.basisSize; i++) {
                for (int j = 0; j < integrals.basisSize; j++) {
                    rmsDens = rmsDens + Math.pow((densMat[i][j] - densOld[i][j]),2);
                }
            }
            rmsDens = Math.sqrt(rmsDens);
            System.out.println("  " + icycle + "  " + energyEl + "  " + energyTot + "  " + deltaE + "  " + rmsDens);

            // Reset "old" energy and density
            energyTotOld = energyTot;
            densOld = densMat;

            icycle = icycle + 1;
            if (icycle > maxCycle) {
                System.err.println("ERROR! The SCF failed to converge!");
                System.exit(1);
            }
        } while (deltaE > eCriterion || rmsDens> densCriterion );

        System.out.println("CONGRATULATIONS! THE SCF HAS CONVERGED!");






        double result = 0;
        return result;
    }



}
