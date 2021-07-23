package com.company;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

public class EigenvectorTest {

    static Matrix instance;

    @BeforeAll
    static void setup(){
        instance = new Matrix();
    }

    @Test
    public void testEigenvector(){
        int dimension = 4;
        double[][] matrixIn = new double[][] {{1,8,2,4},{8,-5,-1,9},{2,-1,-11,7},{4,9,7,5}};
        double[][] matrixOut = instance.eigenVV(matrixIn);

        //System.out.print

        Assertions.assertEquals(0.7468241328124314, val, 1E-12);
    }
}