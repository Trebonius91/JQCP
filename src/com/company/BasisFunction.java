package com.company;

/*
This class contains the entire definition of a single (Gaussian Type) basis function, contained
of an arbitrary number of primitive Gaussians of arbitrary angular momentum
These objects are defined during the read in of the BSE data in BasisSet/setParameters
and subsequently used for all integral calculations
 */

import java.util.Arrays;
import java.lang.Math;

public class BasisFunction {

    // The value of Pi
    private static final double PI = 3.1415926535897932384626;
    // The number of the current basis function (unique identifier)
    int index;
    // The position of the basis function in xyz coordinates (bohr)
    double[] position;
    // The atom index on which the basis function is positioned (zero if anywhere else in space)
    int atomIndex;
    // The element of the atom at which the basis function is positioned
    String atomElement;
    // The three angular momentum quantum numbers (e.g. (0,0,0)=s etc.)
    int[] angularNum;
    // Number of contracted primitive Gaussian functions
    int contraction;
    // Contraction coefficients (prefactors in sum) of primitive Gaussians
    double[] coeffs;
    // Exponents of primitive Gaussians
    double[] exps;
    // The normalization constants of all primitive Gaussians in the current basis function
    double[] norms;
    // Letter to describe the actual basis function (s,px,py,pz etc)
    String symbol;
    // Local object of the class Mathematics for some useful functions etc.
    Mathematics mathUtil;

    // The default constructor: define the needed contractions for the array dimensions

    BasisFunction(int contraction) {

        if (contraction < 1) {
            System.err.println("No basis function can have a contraction of less than one Gaussian!");
        }
        this.index = 0;
        this.contraction = contraction;
        this.position = new double[3];
        this.angularNum = new int[3];
        this.coeffs = new double[contraction];
        this.exps = new double[contraction];
        this.norms = new double[contraction];
        this.atomIndex = 0;
        this.symbol = "";
        this.atomElement = "NULL";
        this.mathUtil = new Mathematics();

    }


    // getter/setter for basis function index

    public void setIndex (int index) {this.index = index;}

    public int getIndex() {
        return index;
    }

    // getter/setter for basis function position

    public void setPosition (double[] xyz) {
        this.position = xyz;
    }

    public double[] getTPosition() {
        return this.position;
    }

    // getter/setter for index of assiged atom

    public void setAtomIndex (int index) { this.atomIndex = index;}

    public int getAtomIndex() { return this.atomIndex;}

    // getter/setter for angular momentum quantum numbers

    public void setAngularNum (int[] numArray) {
        this.angularNum = numArray;
        if (this.angularNum[0] == 1) {
            this.setSymbol("px");
        } else if (this.angularNum[1] == 1) {
            this.setSymbol("py");
        } else if (this.angularNum[2] == 1) {
            this.setSymbol("pz");
        } else {
            this.setSymbol("s");
        }
    }

    public int[] getAngularNum() {return this.angularNum;}

    // getter/setter for number of contracted primitive Gaussians

    public void setContraction (int contract) {this.contraction = contract;}

    public int getContraction() {return this.contraction;}

    // getter/setter for contraction coefficients for primitive Gaussians
    // Set directly the symbol for the basis function!

    public void setCoeffs (double[] contractCoeffs) {
        this.coeffs = contractCoeffs;
    }


    public double[] getCoeffs() {
        return coeffs;
    }
    // getter/setter for exponents of primitive Gaussians

    public void setExps(double[] exps) {
        this.exps = exps;
    }

    public double[] getExps() {
        return exps;
    }

    // getter/setter for symbol of actual basis function

    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    public String getSymbol() {
        return symbol;
    }

    // getter/setter for atom element of actual basis function

    public void setAtomElement(String element) {
        this.atomElement = element;
    }

    public String getAtomElement() {
        return atomElement;
    }

    // getter/setter for the norms of primitive Gaussians in the basis function

    public void setNorms() {
        // check if exponents of primitive Gaussians were already set
        if (this.exps[0] == 0) {
            System.err.println("Please set the exponents of the primitive Gaussians first!");
        } else {
            for (int i = 0; i < this.contraction; i++) {
                int lmn_sum = this.angularNum[0] + this.angularNum[1] + this.angularNum[2];
                this.norms[i] = Math.sqrt(Math.pow(2, (2 * lmn_sum) + 1.5) * Math.pow(this.exps[i], (lmn_sum + 1.5))/
                        (mathUtil.fact(2 * this.angularNum[0] - 1) * mathUtil.fact(2 * 
                        this.angularNum[1] - 1) * mathUtil.fact(2 * this.angularNum[2] - 1) * Math.pow(PI, 1.5)));
            }
        }
    }

    public double[] getNorms() {
        return norms;
    }

    // Test method for checking the content of an arbitrary basis function

    public void checkBasisFunction() {
        System.out.println("----------------------------------------------------------");
        System.out.println("The basis function 'index "+ this.index + "' has the following properties:");
        System.out.println("The angular momentum numbers are: " + Arrays.toString(this.angularNum));
        System.out.println("Thus it is a " + this.symbol + "-function.");
        int printIndex = atomIndex + 1;
        System.out.println("It is localized at atom No. " + printIndex + " (" + this.atomElement + ")");
        System.out.println("The location is " + Arrays.toString(this.position));
        System.out.println("It contains " + this.contraction + " contracted primitive Gaussians:");
        System.out.println("No.      Exponent      Contraction Coefficient    Norm");
        for (int i=0; i < this.contraction; i++) {
            System.out.println(i+1 + "        " + this.exps[i] + "    " + this.coeffs[i] + "    " + this.norms[i]);
        }



    }

}
