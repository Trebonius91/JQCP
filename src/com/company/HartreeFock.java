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

    // ONLY TEMPORARY, AS LONG AS THE EIGENVECTOR CALCULATION DOES NOT WORK
    Eigenvector_external external;


    // The obtained SCF energy
    double SCFEnergy;

    // The obtained MO coefficients
    double[][] MOCoeffs;


    // The default constructor: initialize the dummy objects for useful mathematical methods...
    HartreeFock() {

        this.matrix = new Matrix();
        this.mathUtil = new Mathematics();
        this.external = new Eigenvector_external();

    }

    // The setter method for setting all needed information...
    public void setIntegrals(MolecularIntegral integrals) {
        this.integrals = integrals;
    }

    // The method for applying the Hartree-Fock method to closed shell systems: The Roothaan-Hall equations will
    // solved iteratively during an SCF cycle. The HF energy and the MO coefficients of the molecular system
    // will then be obtained and returned indirectly by storing them to the respective attributes of the current
    // HF object.

    public double SCFClosedShell() {

        int dimension = this.integrals.basisSize;

        // First, calculate the core Hamiltonian matrix (T + V(1e))
        double[][] matHCore = new double[dimension][dimension];

        matHCore = matrix.sum(this.integrals.kineticIntegral,this.integrals.nuclearIntegral);

        // For information, print out the obtained core Hamiltonian as initial guess to file.

        System.out.println("The core hamiltonian");
        matrix.print(matHCore);

        // Diagonalize the overlap matrix in order to get the orthogonalization matrix S

        double[][] dummyMatrix = this.matrix.eigenVV(this.integrals.overlapIntegral);





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

        // Now calculate the intial HF energy

        double energy = 0.0;

        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                energy=energy+densMat[j][i]*(matHCore[j][i]+FZero[j][i]);
            }
        }
        System.out.println("Energy first" + energy);


        double result = 0;
        return result;


    }



}
