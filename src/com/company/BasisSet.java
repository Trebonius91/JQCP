package com.company;


import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/*
This class facilitates the efficient management of the basis set for the molecular system
in use.
It communicates with the Basis Set Exchange (BSE) site (https://www.basissetexchange.org/)
and retrieves the basis function parameters directly from there!
Consequently, each (Gaussian function) based basis set listed there can be taken
for the calculation, as long as all elements in the molecule have parametrized functions.

 */
public class BasisSet {
    // Number of atoms in the system
    int numAtoms;
    // Name of the basis set in usage
    String basisName;
    // Atomic numbers of all atoms in the system
    int[] atomicNumbers;
    // The element symbols of all atoms in the system
    String[] elementSymbols;
    // Nuclear coordinates of the atoms of the system
    double[][] atomCoordinates;
    // The main content of the class: a list with all basis functions!
    BasisFunction[] basisList;
    // The number of basis functions in the basis set
    int basisNumber;

    // The main constructor, mainly a default constructor

    BasisSet() {
        this.numAtoms = 0;
        this.basisName = "NULL";
    }

    // Getter/setter for the number of atoms

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

    // Getter/setter for the name of the basis set

    public void setBasisName (String name) {
        this.basisName = name;
    }

    public String getBasisName () {
        return this.basisName;
    }


    // Getter/setter for the array with atomic numbers

    public void setAtomicNumbers (int numbers[]) {
        if (this.numAtoms < 1) {
            System.err.println("The number of atoms must be greater than 1!");
            new IOException();
        }
        this.atomicNumbers = new int[this.numAtoms];
        for (int i = 0; i < this.numAtoms; i++) {
            if (numbers[i] < 1) {
                System.err.println("The atomic number at position " + i + "is not valid!");
            }
            this.atomicNumbers[i]=numbers[i];
        }
    }

    public int[] getAtomicNumbers() {
        return this.atomicNumbers;
    }


    // Getter/setter for atomic coordinates in the system

    public void setAtomCoordinates (double xyz[][]) {
        if (this.numAtoms < 1) {
            System.err.println("The number of atoms must be greater than 1!");
            new IOException();
        }
        this.atomCoordinates = new double[this.numAtoms][3];
        for (int i = 0; i < this.numAtoms; i++) {
            for (int j = 0; j < 3; j++) {
                this.atomCoordinates[i][j] = xyz[i][j];
            }
        }
    }

    public double[][] getAtomCoordinates() { return this.atomCoordinates; }


    // getter/setter for element symbols of the system


    public void setElementSymbols(String[] elementSymbols) {
        this.elementSymbols = elementSymbols;
    }

    public String[] getElementSymbols() {
        return elementSymbols;
    }

    // getter/setter for number of basis functions (settings is done during setBasisList anyway...)

    public void setBasisNumber(int number) {
        this.basisNumber = number;
    }

    public int getBasisNumber() {
        return basisNumber;
    }

    /*
    The main setter which reads the basis set information from the BSE site and allocates
    array sizes, contractions, parameters etc. from the files
     */

    public void setBasisList() {

        // Check if a reasonable number of atoms is set
        if (this.numAtoms < 1) {
            System.err.println("The number of atoms must be greater than 1!");
            new IOException();
        }

        // Loop through all atoms in the molecule and fetch the basis set information from the internet
        // Here, basis set information is only retrieved for each individual element in the system!
        // The total number of basis functions and other allocations will be done in the following..

        int[][] usedElements = new int[92][2]; // array for element indices already parametrized

        int indexLookup = 0;  // If the current element was already defined: position in the usedElements array
        int[] actElements = this.getAtomicNumbers();
        boolean defined = false; // If the element of the current atom was already defined
        String[] outputs = new String[this.numAtoms];

        int incrementUse = 0;  // how many elements are already defined
        // All basis set names must be in lower case for the URL
        String basisLower = this.basisName.toLowerCase();

        // Global precursor array for storage of all defined basis functions

        List<BasisFunction> basisList = new ArrayList<BasisFunction>();

        int basisSize = 0; // total number of (contracted) basis functions in the molecule

        for (int i = 0; i < this.numAtoms; i++) {
            for (int j = 0; j < incrementUse; j++) {
                if (actElements[i] == usedElements[j][0]) {
                    defined = true;
                    indexLookup = usedElements[j][1];
                    break;
                }

            }

            if (!(defined)) {


                // Call the respective BSE site for the actual element and basis set!
                HttpClient client = HttpClient.newHttpClient();
                HttpRequest request = HttpRequest.newBuilder()
                        .uri(URI.create("https://www.basissetexchange.org/basis/" + basisLower +
                                "/format/orca/?version=1&elements=" + actElements[i]))
                        .GET() // GET is default
                        .build();
                try {
                    HttpResponse<String> response = client.send(request,
                            HttpResponse.BodyHandlers.ofString());
                    //   System.out.println(response.body());
                    outputs[i] = response.body();
                } catch (IOException e) {
                    e.printStackTrace();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
                usedElements[incrementUse][0]= actElements[i];
                usedElements[incrementUse][1]= i;
                incrementUse++;
               // System.out.println("New element!");
            } else {
                outputs[i] = outputs[indexLookup];
            }
            defined = false;
            // Analyze the obtained sourcecode of the site and obtain the basis function parameters
            // if the respective element already occured, copy the output string...
            String[] basisSplit = outputs[i].split("\n");
            String[] basisFooter = new String[250]; // default size for parsing array
            boolean save = false;
            int k = 0;
            for (String basisLine : basisSplit) {
                if (basisLine.equals("$DATA")) {
                    save = true;
                } else if (basisLine.equals("$END")) {
                    save = false;
                }
                if (save) {
                    basisFooter[k] = basisLine;
                    k++;
                }
            }

            int footerLines = k;

            //for (int l = 0; l < footerLines; l++) {
            //   System.out.println(basisFooter[l]);
            //}

            // Now successively define all basis functions in the system by reading them from the website content

            boolean newFunction = false;
            int[][] angularMomentum = new int[10][3];
            int functionNumber = 0;  // the number of basis functions for the current shell

            String readMode = "";  // If only one or two (s/p) basis function types are read in


            for (int l = 0; l < footerLines; l++) {
                String[] currentSplit = basisFooter[l].split("\\s+");

                // consider different types of basis functions/angular momentum quantum numbers
                // S: one single s function
                if (currentSplit[0].equals("S")) {
                    newFunction = true;
                    angularMomentum[0] = new int[]{0, 0, 0};
                    functionNumber = 1;
                    readMode = "simple";
                // L: one s and three p functions (px, py, pz)
                } else if (currentSplit[0].equals("L")) {
                    newFunction = true;
                    angularMomentum[0] = new int[]{0, 0, 0};  // s
                    angularMomentum[1] = new int[]{1, 0, 0};  // px
                    angularMomentum[2] = new int[]{0, 1, 0};  // py
                    angularMomentum[3] = new int[]{0, 0, 1};  // pz
                    functionNumber = 4;
                    readMode = "double";
                // P: three p functions (px, py, pz)
                } else if (currentSplit[0].equals("P")) {
                    newFunction = true;
                    angularMomentum[0] = new int[]{1, 0, 0};  // px
                    angularMomentum[1] = new int[]{0, 1, 0};  // py
                    angularMomentum[2] = new int[]{0, 0, 1};  // pz
                    functionNumber = 3;
                    readMode = "simple";
                // D: five d functions (dz2, dxz, dyz, dx2y2, dxy)
                } else if (currentSplit[0].equals("D")) {
                    newFunction = true;
                    angularMomentum[0] = new int[]{0, 0, 2};  // dz2
                    angularMomentum[1] = new int[]{1, 0, 1};  // dxz
                    angularMomentum[2] = new int[]{0, 1, 1};  // dyz
                    angularMomentum[3] = new int[]{2, 2, 0};  // dx2y2
                    angularMomentum[4] = new int[]{1, 1, 0};  // dxy
                    functionNumber = 5;
                    readMode = "simple";
                }
                int actContraction = 0;
                // Read in the degree of contraction as central parameter!
                if (newFunction) {

                    actContraction = Integer.parseInt(currentSplit[1]);
                }
                // only add a new basis function if a line beginning with S/L etc was parsed
                if (newFunction) {
                    for (int n = 0; n < functionNumber; n++) {


                        // Generate the new basis function
                        BasisFunction actGaussian = new BasisFunction(actContraction);
                        actGaussian.setIndex(basisSize);
                        actGaussian.setPosition(this.atomCoordinates[i]);
                        actGaussian.setAngularNum(angularMomentum[n]);
                        actGaussian.setAtomIndex(i);
                        actGaussian.setAtomElement(this.elementSymbols[i]);

                        // Loop locally through the next Ncontraction lines to extract the exponents and contractions

                        double[] localCoeffs = new double[actContraction];
                        double[] localExps = new double[actContraction];

                        // For the usual formate: only one type of basis function per segment
                        if (readMode.equals("simple")) {
                            for (int m = 0; m < actContraction; m++) {
                                currentSplit = basisFooter[l + m + 1].split("\\s+");
                                localCoeffs[m] = Double.parseDouble(currentSplit[2]);
                                localExps[m] = Double.parseDouble(currentSplit[1]);
                            }
                            // For the special type for STO/6-31G... : s and p functions in three columns
                        } else if (readMode.equals("double")) {
                            if (actGaussian.getSymbol().equals("s")) {
                                for (int m = 0; m < actContraction; m++) {
                                    currentSplit = basisFooter[l + m + 1].split("\\s+");
                                    localCoeffs[m] = Double.parseDouble(currentSplit[2]);
                                    localExps[m] = Double.parseDouble(currentSplit[1]);
                                }
                            } else {
                                for (int m = 0; m < actContraction; m++) {
                                    currentSplit = basisFooter[l + m + 1].split("\\s+");
                                    localCoeffs[m] = Double.parseDouble(currentSplit[3]);
                                    localExps[m] = Double.parseDouble(currentSplit[1]);
                                }
                            }
                        }

                        // Parametrize contraction coefficients and primitive Gaussians of actual Basis function
                        actGaussian.setCoeffs(localCoeffs);
                        actGaussian.setExps(localExps);
                        actGaussian.setNorms();

                        actGaussian.checkBasisFunction();

                        // Add basis function to array
                        basisList.add(actGaussian);
                        // increment the basis function index
                        basisSize++;
                    }



                }
                newFunction = false;
            }



        }
        // Finally convert the intermediary linked list to the desired array attribute of the BasisSet class
        this.basisList = new BasisFunction[basisSize];
        basisList.toArray(this.basisList);
        // Also set the number of basis functions in the set
        this.basisNumber = basisSize;


    }
    // The getter function for the whole basis set of BasisFunction objects

    public BasisFunction[] getBasisList() {
        return basisList;
    }
}
