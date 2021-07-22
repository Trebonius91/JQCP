package com.company;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

public class BoysTest {

    static Boys instance;

    @BeforeAll
    static void setup(){
        instance = new Boys();
    }

    @Test
    public void testBoys(){
        double val = instance.boysFunction(0,1.0);
        Assertions.assertEquals(0.7468241328124314, val, 1E-12);
    }
}
