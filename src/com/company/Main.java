package com.company;
import javax.swing.*;
import java.io.FileWriter;
import java.io.IOException;
import org.apache.commons.math3.*;

public class Main {


    public static void main(String[] args) throws IOException {
	// write your code here
        System.out.println( "Hello World!");
    //    System.err.println( "Hello World!");

        //String s = JOptionPane.showInputDialog("Bitte Zahl eingeben");
        //int i = Integer.parseInt(s);
        //System.out.println(i*i);
       /*  int test;
        try {
            test =
        } */
        // the default input file for the QM calculations
        String filename="test/water/input.txt";
        System.out.println("Working Directory = " + System.getProperty("user.dir"));
        // Define a Math dummy object for several small function calls
        Mathematics mathUtil = new Mathematics();

        Molecule molecule = new Molecule(filename);

        molecule.checkContent();

        BasisSet basisSet = new BasisSet();

        basisSet.setNumAtoms(molecule.getNumAtoms());
        basisSet.setBasisName(molecule.basisName);
        basisSet.setAtomicNumbers(molecule.getAtomicNumbers());
        basisSet.setElementSymbols(molecule.elements);

        basisSet.setAtomCoordinates(molecule.getStructure());

        // Finally, set the whole basis set with all functions etc
        basisSet.setBasisList();

        // Calculate the molecular integrals
        // first, setup the object
        MolecularIntegral molecularIntegral = new MolecularIntegral(basisSet.basisNumber, molecule.numAtoms);

        molecularIntegral.setAtomCoordinates(molecule.getStructure());
        molecularIntegral.setAtomicNumbers(molecule.atomicNumbers);
        molecularIntegral.setFullBasis(basisSet.basisList);
        // then, calculate the overlap integral matrix

        molecularIntegral.setOverlapIntegral();
        molecularIntegral.setKineticIntegral();
        molecularIntegral.setNuclearIntegral();
        molecularIntegral.setRepulsionIntegral();

        // check Boys function

        /*FileWriter myWriter = new FileWriter("boys.dat");
        for (int i = 1; i <= 10; i++) {
            for (int k = 1; k <= 50; k++) {
                double boys = molecularIntegral.boysFunction(i-1,k*0.1);
                myWriter.write(i + " " + k + " " + boys + "\n");
            }
        }
         myWriter.close(); // please close me! */



        //double result = boys -2.5592187139634830E-002;
        //double boys = boysUnit.boysFunction(5,0.014);
        //System.out.println("boys " + boys);

        //Matrix matrix = new Matrix();


        //double[][] matrixD = new double[][] {{-261,209,-49},{-530,422,-98},{-800,631,-144}};

        //double[][] matrixD = new double[][] {{52,30,49,28},{30,50,8,44},{49,8,46,16},{28,44,16,22}};


        //matrix.print(vectorB);
        //double[][][] dummy = matrix.decomposeQR(matrixD);

        //double[][] eigenvalue = matrix.eigenVV(matrixD);
        //matrix.print(eigenvalue);


        // Call the Hartree Fock routine for the closed shell SCF
        HartreeFock hartreeFock = new HartreeFock();
        hartreeFock.setIntegrals(molecularIntegral);


        hartreeFock.SCFClosedShell();


    }


}
