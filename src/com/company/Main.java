package com.company;
import javax.swing.*;
import java.io.FileWriter;
import java.io.IOException;
import org.apache.commons.math3.*;

public class Main {


    public static void main(String[] args) throws IOException {

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

        System.out.println("Set the required basis set for the system... ");
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
        System.out.println(" ");
        System.out.println("Calculate the required molecular integrals...");
        MolecularIntegral molecularIntegral = new MolecularIntegral(basisSet.basisNumber, molecule.numAtoms);

        molecularIntegral.setAtomCoordinates(molecule.getStructure());
        molecularIntegral.setAtomicNumbers(molecule.atomicNumbers);
        molecularIntegral.setFullBasis(basisSet.basisList);
        // then, calculate the overlap integral matrix

        System.out.println("Calculate the overlap integrals...");
        molecularIntegral.setOverlapIntegral();
        System.out.println("Calculate the kinetic energy integrals...");
        molecularIntegral.setKineticIntegral();
        System.out.println("Calculate the electron-nuclear attraction integrals...");
        molecularIntegral.setNuclearIntegral();
        System.out.println("Calculate the electron-electron repulsion integrals...");
        molecularIntegral.setRepulsionIntegral();

        System.out.println("... done!");
        System.out.println(" ");

        // Call the Hartree Fock routine for the closed shell SCF

        System.out.println("Start the Hartree Fock SCF cycle...");
        HartreeFock hartreeFock = new HartreeFock();
        hartreeFock.setIntegrals(molecularIntegral);


        hartreeFock.SCFClosedShell();


    }


}
