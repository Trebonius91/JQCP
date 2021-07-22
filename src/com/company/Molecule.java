package com.company;

import java.io.FileNotFoundException;
import java.io.File;
import java.util.Scanner;

/*
This is the central class definition of the program package.
The chemical system to be computed is entirely translated into one single instance
of the Molecule class. Attributes are all properties needed to describe a chemical system
in a quantum mechanical sense.
Methods are all input/output as well as calculation routines for all the different
available quantum chemical methods: HF, MP2, CCSD(T) etc...
 */
public class Molecule {
    // The number of atoms in the system
    int numAtoms;
    // The total charge on the system, in multiples of the elementary charge
    int totCharge;
    // The spin multiplicity of the system (2S+1)
    int totSpin;
    // The name of the basis set to be calculated with
    String basisName;
    // The element symbols as Strings
    String elements[];
    // The element symbols as numbers (atomic numbers)
    int atomicNumbers[];
    // The xyz structure in bohr (to be converted from the Angstrom input/output)
    double structure[][];
    // The calculated HF energy of the system
    double energyHF;

    // The constructor which reads in the information from an input file
    Molecule(String filename) {

        try {
            File myFile = new File(filename);
            System.out.println(filename);
            Scanner myReader = new Scanner(myFile);


            while (myReader.hasNextLine()) {
                String currentLine = myReader.nextLine();
                String[] currentSplit = currentLine.split("\\s+");
                for (int i = 0; i < currentSplit.length; i++) {
                    // First, read in all general information (atom number etc.)
                    if (currentSplit[i].equals("%natoms")) {
                        this.numAtoms = Integer.valueOf(currentSplit[i+1]);
                        this.elements = new String[this.numAtoms];
                        this.atomicNumbers = new int[this.numAtoms];
                        this.structure = new double[this.numAtoms][3];
                    } else if (currentSplit[i].equals("%charge")) {
                        this.totCharge = Integer.valueOf(currentSplit[i+1]);
                    } else if (currentSplit[i].equals("%mult")) {
                        this.totSpin = Integer.valueOf(currentSplit[i+1]);
                    } else if (currentSplit[i].equals("%basis")) {
                        this.basisName = currentSplit[i+1];
                        // Second, read in the structure and elements of the molecule
                    } else if (currentSplit[i].equals("%geometry")) {
                        for (int j = 0; j < this.numAtoms; j++) {
                            String structureLine = myReader.nextLine();
                            String[] structureSplit = structureLine.split("\\s+");
                            elements[j]=structureSplit[0];
                            structure[j][0]=Double.valueOf(structureSplit[1]);
                            structure[j][1]=Double.valueOf(structureSplit[2]);
                            structure[j][2]=Double.valueOf(structureSplit[3]);
                        }
                    }
                }

            }
            myReader.close();

        } catch (FileNotFoundException e) {
            System.out.println("The input file " + filename  +" could not be handled!");
        }
        // Finally set the atomic numbers by converting them from the element symbols
        this.setAtomicNumbers();
        // Do a control output of the read-in molecule to the command line


    }

    public void checkContent() {
        System.out.println("The defined molecule has the following properties:");
        System.out.println("The number of atoms is: " + this.numAtoms);
        System.out.println("The total charge is: " + this.totCharge);
        System.out.println("The spin multiplicity is: " + this.totSpin);
        System.out.println("The basis set to be calculated with is: " + this.basisName);
        System.out.println("The structure of the molecule (Angstrom) is:");
        for (int i = 0; i < this.numAtoms; i++) {
            System.out.println(" " + elements[i] + " " + "(" + atomicNumbers[i] + ")" + " " + structure[i][0] +
                    " " + structure[i][1] + " " + structure[i][2]);
        }
    }

    //  Define several getter and setter methods...

    //    1. The number of atoms

    public void setNumAtoms (int num) {
        if (num > 0) {
            this.numAtoms = num;
        } else {
            System.err.println("The number of atoms must be at least one!");
        }
    }

    public int getNumAtoms() {
        return this.numAtoms;
    }

    //   2. The charge of the system

    public void setTotCharge (int num) {
        this.totCharge = num;
    }

    public int getTotCharge() {
        return this.totCharge;
    }

    //   3. The spin multiplicity of the system

    public void setTotSpin (int num) {
        if (num > 1) {
            this.totSpin = num;
        } else {
            System.err.println("The spin multiplicity can never be negative!");
        }
    }

    public int getTotSpin() {
        return this.totSpin;
    }

    //   4. The structure of the system

    public void setStructure (double xyz[][]) {
        boolean atomError = false;
        boolean dimensionError = false;
        // Check if the submitted structure has the correct dimensions!
        if (xyz.length != this.structure[0].length) {
            System.err.println("The number of atoms in the given array is not equal to that in the molecule!");
            atomError = true;
        }
        for (int i = 0; i < xyz[0].length; i++) {
            if (xyz[i].length != 3) {
                System.err.println("The dimension of the submitted array seems to be corrupted in line " + i + "!");
                dimensionError = true;
            }
        }
        if (!(atomError) && !(dimensionError)) {
            this.structure = xyz;
        } else {
            System.err.println("The structure of the molecule could not be modified!");
        }

    }

    public double[][] getStructure() {
        return this.structure;
    }

    // Method to obtain list with element numbers of all atoms in the molecule

    public void setAtomicNumbers() {
        // Define array for mapping of element symbols to atomic numbers
        String[] atomConverter = new String[] {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",
             "Mg", "Al", "Si", "P", "S", "Cl", "Ar"};

        // Loop through atomic numbers to check which element is present and set entry of array
        for (int i = 0; i < this.numAtoms; i++) {
            int k = 0;
            for (String symbol : atomConverter) {
                k++;
                if (this.elements[i].equals(symbol)) {
                    atomicNumbers[i] = k;
                }
            }
        }
    }

    public int[] getAtomicNumbers() {
        return this.atomicNumbers;
    }


}
