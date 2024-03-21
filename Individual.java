package com.example.chessboard;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Individual {


    String[] parent ;
    float fitness;
   //static int randomChromosome ;



    public Individual(String[] parent , float fitness){
        this.parent = parent;
        this.fitness = fitness ;
    }
    /*public Individual(int randomChromosome){
        this.randomChromosome = randomChromosome ;
    }*/
    public Individual(String[] parent){
        this.parent = parent;
    }

    public float getFitness() {
        return fitness;
    }

    public void setFitness(float fitness) {
        this.fitness = fitness;
    }



    public String[] getParent() {
        return parent;
    }

    public void setParent(String[] parent) {
        this.parent = parent;
    }



}
