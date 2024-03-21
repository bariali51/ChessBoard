package com.example.chessboard;
import java.text.DecimalFormat;
import java.util.*;

public class GeneticAlgorithm {

    public static final int POPULATION_SIZE = 100;
    public static final int CHROMOSOME_SIZE = 8;
    public static final int CHROMOSOME_ARRAY_SIZE = 64;
    public static final float CROSSOVER_PROBABILITY = 0.8F;
    public static final float CHROMOSOME_PROBABILITY = 0.1F;
    // DecimalFormat df = new DecimalFormat("#.##");


    public void createPopulation() {


        List<String[]> maxAllPopulaion = new ArrayList<>();
        List<Float> maxAllPopulaionFintness = new ArrayList<>();
        List<Individual> orderMaxPopulation = new ArrayList<>();

for (int newPopulationIndice = 0 ; newPopulationIndice < POPULATION_SIZE ; newPopulationIndice++){
        // applique Genetic Algorithme
        System.out.println("--------- POPULATIONS_0"+newPopulationIndice+ "----------");

        List<Individual> population = createChromosome();
        int indice = 0;
        System.out.println("-------- Parent(Root)");
        while (indice < population.size()) {
            String[] chromosome = population.get(indice).parent;
            System.out.print(population.get(indice).fitness + "-->  |");
            for (int i = 0; i < chromosome.length; i++) {
                System.out.print(chromosome[i] + "|");
            }
            System.out.println();
            indice++;
        }
        System.out.println();

        List<Individual> order = new ArrayList<>();
        List<String[]> listMutation = new ArrayList<>();
        List<Float> listFitness = new ArrayList<>();
        List<Float> listFitnessMutation = new ArrayList<>();
        List<List<Individual>> listChild = new ArrayList<>();
        List<String[]> newPopulation = new ArrayList<>();


        for (int i = 0; i < CHROMOSOME_SIZE; i++) {

            listChild.add(singlePointCrossover(population));
        }
        for (int i = 0; i < CHROMOSOME_SIZE; i++) {
            newPopulation.add(listChild.get(i).get(0).parent);
            listFitness.add(listChild.get(i).get(0).fitness);
            newPopulation.add(listChild.get(i).get(1).parent);
            listFitness.add(listChild.get(i).get(1).fitness);
        }


        for (int i = 0; i < newPopulation.size(); i++) {
            System.out.print(listFitness.get(i) + "-->  |");
            for (int j = 0; j < 64; j++) {
                System.out.print(newPopulation.get(i)[j]);
            }
            System.out.println();
        }

        for (int i = 0; i < newPopulation.size(); i++) {

            //  System.out.println("--------- New_POPULATION "+(i+1)+"----------");

            listMutation.add(mutationOperator(newPopulation, i));

        }

        // Calculate Fitness
        for (int i = 0; i < newPopulation.size(); i++) {
            String[][] conver2D = covert1DTo2DChromosome(listMutation.get(i));
            listFitnessMutation.add(calculate_Fitness(conver2D));
        }
        Collections.sort(listFitness.reversed());

        // Order Chromosome By Fitness

        order = orderChromosome(listMutation, listFitnessMutation);

        //  System.out.println("------------- MaxPopulation_0"+newPopulationIndice+" --------");

        System.out.print(order.get(0).fitness + "--->  |");
        maxAllPopulaion.add(order.get(0).parent);
        maxAllPopulaionFintness.add(order.get(0).fitness);
        for (int j = 0; j < CHROMOSOME_ARRAY_SIZE; j++) {
            System.out.print(order.get(0).parent[j]);
        }

        System.out.println();

        order.clear();
        listFitnessMutation.clear();
        listFitness.clear();
        listMutation.clear();
        listChild.clear();
        population.clear();
        newPopulation.clear();


    }


    // order

       orderMaxPopulation = orderChromosome(maxAllPopulaion,maxAllPopulaionFintness);

        System.out.println("---------- AllMaxPopilation ------");

        for (int i = 0 ; i < orderMaxPopulation.size(); i++) {
          //  System.out.print(orderMaxPopulation.get(i).fitness + "--->  |");
            for (int j = 0; j < CHROMOSOME_ARRAY_SIZE; j++) {
                System.out.print(orderMaxPopulation.get(i).parent[j]);
            }
            System.out.println();
        }


  System.out.println();
   System.out.println("--------_Best Chromosome in Gentic Algorithm_-----");
        System.out.print(orderMaxPopulation.get(0).fitness + "--->  |");
        for (int j = 0; j < CHROMOSOME_ARRAY_SIZE; j++) {
            System.out.print(orderMaxPopulation.get(0).parent[j]);
        }
        System.out.println();

        System.out.println("--------_Chromosome 2D Final_-----");
        // Convert 1D to 2D
    String[][] chromsome2D = covert1DTo2DChromosome(orderMaxPopulation.get(0).parent);
        for (int row = 0; row < 8; row++) {
            for (int col = 0; col < 8; col++) {
                System.out.print(chromsome2D[row][col]);
            }
            System.out.println();
        }




}


    public String[][] covert1DTo2DChromosome(String[] chromosomes1D) {
        String[][] converteur = new String[8][8];

        int index = 0;
        for (int row = 0; row < 8; row++) {
            for (int col = 0; col < 8; col++) {
                converteur[row][col] = chromosomes1D[index++];
            }
        }



        return converteur;
    }

    // Selction

    public List<Individual> selctionChromosome(List<Individual> population) {
        List<Individual> selectedChromosomes = new ArrayList<>();

        // Calculate the total fitness of the population
        double totalFitness = 0.0;
        for (Individual chromosome : population) {
            totalFitness += chromosome.fitness;
        }
        //System.out.println("\nTotalFitness -------------------- >" + totalFitness);
        Random random = new Random();
        int numSelections = 2;
        // Perform roulette wheel selection
        for (int i = 0; i < numSelections; i++) {
            double spin = random.nextDouble() * totalFitness; // Spin

            double partialSum = 0.0;
            for (Individual chromosome : population) {
                partialSum += chromosome.fitness;
                if (partialSum >= spin) {
                    selectedChromosomes.add(chromosome);
                    break;
                }
            }
        }
        return selectedChromosomes;
    }


    // ----- crossOverMethode ---

    public List<Individual> singlePointCrossover(List<Individual> newGeneration) {
        List<Individual> selectedChromosomes = selctionChromosome(newGeneration);

        List<Individual> offspringChromosomes = new ArrayList<>();
        System.out.println("----- Cohose Tow Fitness -----");
        for (Individual chromosome : selectedChromosomes) {
            System.out.println("Fitness: " + chromosome.fitness);
        }

        // Perform crossover on selected chromosomes
        offspringChromosomes = crossover(selectedChromosomes.get(0), selectedChromosomes.get(1));
        System.out.println(offspringChromosomes.size());

        return offspringChromosomes;
    }

    public List<Individual> crossover(Individual parent1, Individual parent2) {
        String[] chromsomeX = parent1.parent;
        String[] chromsomeY = parent2.parent;
        List<Individual> offspringChromosomes = new ArrayList<>(); // ----- Childs ------

        int length = chromsomeX.length;

        Random random = new Random();
        float crossoverPoint = random.nextFloat(1); // Generate a random crossover point
        // System.out.println("--------------------_RandomCrossOver_------------- >: "+crossoverPoint);

        if (crossoverPoint < CROSSOVER_PROBABILITY) {

            List<String> chromsome1Left = new ArrayList<>();
            List<String> chromsome2Left = new ArrayList<>();
            List<String> chromsome1Right = new ArrayList<>();
            List<String> chromsome2Right = new ArrayList<>();

            String[] child1 = new String[length];
            String[] child2 = new String[length];

            // Perform single-point crossover
            for (int i = 0; i < (chromsomeX.length / 2); i++) {
                chromsome1Left.add(chromsomeX[i]);
                chromsome2Left.add(chromsomeY[i]);
            }
            for (int i = (chromsomeX.length / 2); i < (chromsomeX.length); i++) {
                chromsome1Right.add(chromsomeX[i]);
                chromsome2Right.add(chromsomeY[i]);
            }


            System.out.println("-------------------_ YES_CrossOver Single Point _-------------------------");

            //System.out.println("chromsome1Left----->"+chromsome1Left);
            //System.out.println("chromsome2Left----->" +chromsome2Left);
            //System.out.println("chromsome1Right---->"+chromsome1Right);
            // System.out.println("chromsome2Right---->"+chromsome2Right);

            chromsome1Left.addAll(chromsome2Right);
            chromsome2Left.addAll(chromsome1Right);
            for (int j = 0; j < chromsome1Left.size(); j++) {
                child1[j] = chromsome1Left.get(j); // Child1
                child2[j] = chromsome2Left.get(j); // Child2
            }


            // ---------- Covert newChromosome 1D to 2D --------
            String[][] new1Chromosome2D = covert1DTo2DChromosome(child1);
            String[][] new2Chromosome2D = covert1DTo2DChromosome(child2);


            // ---  Calculate Fitness newChromosome -----
            float springGenes1Fitness = calculate_Fitness(new1Chromosome2D);
            float springGenes2Fitness = calculate_Fitness(new2Chromosome2D);
            System.out.println("(" + parent1.fitness + "," + parent2.fitness + ")");

            offspringChromosomes.add(new Individual(child1, springGenes1Fitness)); // Set fitness to 0.0 for now
            offspringChromosomes.add(new Individual(child2, springGenes2Fitness)); // Set fitness to 0.0 for now


            System.out.print("----Parent1_--(" + parent1.fitness + ")--->|");
            for (int i = 0; i < chromsomeX.length; i++) {
                System.out.print(chromsomeX[i]);
            }
            System.out.print("|----Child1_----(" + (springGenes1Fitness) + ")---->|");
            for (int i = 0; i < chromsomeX.length; i++) {
                System.out.print(offspringChromosomes.get(0).parent[i]);
            }
            System.out.println();
            System.out.print("----Parent2_--(" + parent2.fitness + ")--->|");
            for (int i = 0; i < chromsomeY.length; i++) {
                System.out.print(chromsomeY[i]);
            }
            System.out.print("|----Child2_----(" + springGenes2Fitness + ")---->|");
            for (int i = 0; i < chromsomeY.length; i++) {
                System.out.print(offspringChromosomes.get(1).parent[i]);
            }
            System.out.println();


            chromsome1Left.clear();
            chromsome2Left.clear();
            chromsome1Right.clear();
            chromsome2Right.clear();
            return offspringChromosomes;

        } else {
            offspringChromosomes.add(new Individual(chromsomeX, parent1.fitness)); // Child 1
            offspringChromosomes.add(new Individual(chromsomeY, parent2.fitness)); // Child 2
            return offspringChromosomes;
        }


    }


    // mutation

    public String[] mutationOperator(List<String[]> newPopulation,int i) {
        String[] chromosomeMutation = new String[CHROMOSOME_ARRAY_SIZE];

        Random random = new Random();
        float p = (random.nextFloat(1));
            // System.out.println("\n------------------_Mutation_------------------");

            //  System.out.println("Probability--------->:"+p);
            if (p < CHROMOSOME_PROBABILITY) {
                System.out.println("\n--------------------- YES_Mutation_Chromosome: _ -------------");
                //System.out.println("ProbabilityChromosome _"+steps+"_----->:"+p);
                Random random1 = new Random();
                int x = random1.nextInt(CHROMOSOME_ARRAY_SIZE);
                int y = random1.nextInt(CHROMOSOME_ARRAY_SIZE);

                System.out.println("(" + (x) + "," + (y) + ")");
                String[] chromosome = newPopulation.get(i);
                System.out.println("(" + chromosome[x] + " , " + chromosome[y] + ")");
                System.out.println("----- Befor Mutation -----");
                for (int a = 0; a < chromosome.length; a++) {
                    System.out.print(chromosome[a]);
                }
                // Change Tow Genes on Chromosome
                String temp = chromosome[x];
                chromosome[x] = chromosome[y];
                chromosome[y] = temp;
                System.out.println();


                System.out.println("----- After Mutation -----");
                for (int k = 0; k < chromosome.length; k++) {
                    chromosomeMutation[k] = chromosome[k];
                }
                for (int k = 0; k < chromosome.length; k++) {
                    System.out.print(chromosomeMutation[k]);
                }
                System.out.println();
                // ----- newChromosome 2D ----
                String[][] new1Chromosome2D = covert1DTo2DChromosome(chromosome);
                // ---  Calculate Fitness newChromosome -----
                float chromosmeFitness = calculate_Fitness(new1Chromosome2D);
                return chromosomeMutation;
            } else {
                return newPopulation.get(i);
            }

   }







    public List<Individual> createChromosome(){


    List<Individual> listIndividual = new ArrayList<>();




    Random random = new Random();
    int numberChromosomes = CHROMOSOME_SIZE/*(1 + random.nextInt(5))*2*/;

  int indice = 0 ;
   while (indice!=numberChromosomes) {

     String[][] chromosome2D = createRandomlyChromosome2D();

       // Calculate Fitness
     float fitness = calculate_Fitness(chromosome2D);
     //System.out.println("\n-------------------------------------- >>Fitness:"+fitness);


     String chromosome[] = new String[CHROMOSOME_ARRAY_SIZE];

      int rows = chromosome2D.length;
      int cols = chromosome2D[0].length;
      int index = 0;
      for (int row = 0; row < rows; row++) {
          for (int col = 0; col < cols; col++) {
              chromosome[index++] = chromosome2D[row][col];
          }
      }




       listIndividual.add(new Individual(chromosome,fitness));


      indice++;
  }






   // Order Ensemble Chromosome1D on Stack by Fitness
    //System.out.println(listFitness);
   // List<String[]> listChromsomeOrder =  orderChromosome(listChromosome,listFitness);
   // List<Individual> listChromsomeOrder =  orderChromosome(listIndividual);
   // Collections.sort(listFitness.reversed());
   // System.out.println(listFitness);

    System.out.println("--------------------------------------------------------------------------------");






return listIndividual ;}


    public  List<Individual>orderChromosome(List<String[]> listStrings ,List<Float> listNumbers){

        List<Individual> pairs = new ArrayList<>();
        for (int i = 0; i < listStrings.size(); i++) {
            String str[] = listStrings.get(i);
            float number = listNumbers.get(i);
            pairs.add(new Individual(str, number));
        }

        // Sort the list of custom objects based on the numbers in descending order
        pairs.sort(Comparator.comparing((Individual pair) -> pair.fitness).reversed());

        // Extract the sorted strings from the list of custom objects
       // List<String[]> sortedStrings = new ArrayList<>();
        List<Individual> sortedStrings = new ArrayList<>();
       for (Individual pair : pairs) {
            sortedStrings.add(pair);
        }

  return sortedStrings; }



    public  String[][] createRandomlyChromosome2D() {

        String[][] chessBoard = new String[8][8];
        for (int row = 0; row < chessBoard.length; row++) { // --->
            for (int col = 0; col < chessBoard.length; col++) {  //  |
               chessBoard[row][col] = "E";
            }
            }






        // ------ Shuffle ----


 List<String> elements1 = new ArrayList<>();
 List<String> elements2 = new ArrayList<>();



        elements1.add("B1");
        elements1.add("R1");
        elements1.add("K1");
        elements1.add("Q1");
        elements1.add("Q1");

        elements2.add("B2");
        elements2.add("R2");
        elements2.add("K2");
        elements2.add("Q2");
        elements2.add("Q2");





     while ((elements1.size() + elements2.size()) < (chessBoard.length*chessBoard[0].length)){
           elements1.add("E");
           elements2.add("E");
     }




        // Shuffle the elements list
        Collections.shuffle(elements1);
        Collections.shuffle(elements2);
        elements1.addAll(elements2);

        int half = (chessBoard.length/2);
        int elementIndex = 0;
       for (int row = 0; row < chessBoard.length; row++) {
            for (int col = 0; col <chessBoard.length; col++) {
                    chessBoard[row][col] = elements1.get(elementIndex++);

            }


        }





        elements1.clear();
        elements2.clear();

        return chessBoard;}


  // Calculate Fitness
 public float  calculate_Fitness(String[][] chromosome2D) {


     float fitness = 0;
     float penalty  = 0 ;
     float nb_conflictsTotal;
     float nb_conflictsQueen1 = 0;
     float nb_conflictsQueen2 = 0;
     float nb_conflictsRook1 = 0;
     float nb_conflictsRook2 = 0;
     float nb_conflictsKinght1 = 0;
     float nb_conflictsKinght2 = 0;
     float nb_conflictsBishop1 = 0;
     float nb_conflictsBishop2 = 0;

     for (int row = 0; row < chromosome2D.length; row++) {
         for (int col = 0; col < chromosome2D.length; col++){

             // ----- Queen1 ------

             if (chromosome2D[row][col].equals("Q1")) {
               //  System.out.println("Q1(" + row + "," + col + ")");

                 //   -----Row ----
                 for (int x = 0; x < chromosome2D.length; x++) {

                    /**/
                      if ((chromosome2D[x][col].equals("Q2"))) {
                             //    System.out.println("YES(Q1)(Row:"+x+") -->  Q2");

                                 nb_conflictsQueen1++;


                     }
                     else if (chromosome2D[x][col].equals("R1")) {
                       //  System.out.println("YES(Q1)(Row) -->  R1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][col].equals("R2")) {
                        // System.out.println("YES(Q1)(Row) -->  R2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][col].equals("B1")) {
                        // System.out.println("YES(Q1)(Row) -->  B1");
                         nb_conflictsQueen1++;
                     } else if (chromosome2D[x][col].equals("B2")) {
                         //System.out.println("YES(Q1)(Row) -->  B2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][col].equals("K1")) {
                        // System.out.println("YES(Q1)(Row) -->  K1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][col].equals("K2")) {
                         //System.out.println("YES(Q1)(Row) -->  K2");
                         nb_conflictsQueen1++;
                     }
                 }

                 // ---- Colmun ----


                 for (int x = 0; x < chromosome2D.length; x++) {
                      if (chromosome2D[row][x].equals("Q1")) {

                             //System.out.println("YES(Q1)(Colmun) -->  Q1");
                             nb_conflictsQueen1++;
                         }
                      else if (chromosome2D[row][x].equals("Q2")) {
                        // System.out.println("YES(Q1)(Colmun) -->  Q2");

                         nb_conflictsQueen1++;
                     }
                      else if (chromosome2D[row][x].equals("R1")) {
                         // System.out.println("YES(Q1)(Colmun) -->  R1");
                          nb_conflictsQueen1++;
                      }
                      else if (chromosome2D[row][x].equals("R2")) {
                         // System.out.println("YES(Q1)(Colmun) -->  R2");
                          nb_conflictsQueen1++;
                      }
                      else if (chromosome2D[row][x].equals("B1")) {
                         // System.out.println("YES(Q1)(Colmun) -->  B1");
                          nb_conflictsQueen1++;
                      } else if (chromosome2D[row][x].equals("B2")) {
                         // System.out.println("YES(Q1)(Colmun) -->  B2");
                          nb_conflictsQueen1++;
                      }

                      else if (chromosome2D[row][x].equals("K1")) {
                        // System.out.println("YES(Q1)(Colmun) -->  K1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[row][x].equals("K2")) {
                        // System.out.println("YES(Q1)(Colmun) -->  K2");
                         nb_conflictsQueen1++;
                     }
                 }

                 //  ----- Diagonal ----

                 // upper Left Diagonal

                 for (int x = row, y = col; x >= 0 && y >= 0; x--, y--){

                      if (chromosome2D[x][y].equals("Q2")) {
                       //  System.out.println("YES(Q1)(upper_Left) -->  Q2");

                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("R1")) {
                        // System.out.println("YES(Q1)(upper_Left) -->  R1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("R2")) {
                         //System.out.println("YES(Q1)(upper_Left) -->  R2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("B1")) {
                         //System.out.println("YES(Q1)(upper_Left) -->  B1");
                         nb_conflictsQueen1++;
                     } else if (chromosome2D[x][y].equals("B2")) {
                       //  System.out.println("YES(Q1)(upper_Left) -->  B2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("K1")) {
                         //System.out.println("YES(Q1)(upper_Left) -->  K1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("K2")) {
                        // System.out.println("YES(Q1)(upper_Left) -->  K2");
                         nb_conflictsQueen1++;
                     }

                 }

                 // lower left Diagonal

                 for (int x = row, y = col; x < chromosome2D.length && y >= 0; x++, y--){
                    /* if (chromosome2D[x][y].equals("Q1")) {
                         queen11 = true ;
                            // System.out.println("YES(Q1)(lower_Left) -->  Q1");
                             nb_conflictsQueen1++;
                         }*/
                      if (chromosome2D[x][y].equals("Q2")) {

                        // System.out.println("YES(Q1)(lower_Left) -->  Q2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("R1")) {
                        // System.out.println("YES(Q1)(lower_Left) -->  R1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("R2")) {
                        // System.out.println("YES(Q1)(lower_Left) -->  R2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("B1")) {
                        // System.out.println("YES(Q1)(lower_Left) -->  B1");
                         nb_conflictsQueen1++;
                     } else if (chromosome2D[x][y].equals("B2")) {
                        // System.out.println("YES(Q1)(lower_Left) -->  B2");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("K1")) {
                        // System.out.println("YES(Q1)(lower_Left) -->  K1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("K2")) {
                        // System.out.println("YES(Q1)(lower_Left) -->  K2");
                         nb_conflictsQueen1++;
                     }

                 }






             }




             // --- Queen 2 ----

             else if (chromosome2D[row][col].equals("Q2")){
                // System.out.println("Q2("+row+","+col+")");

                 //   -----Row ----

                 for (int x = 0; x < chromosome2D.length; x++) {

                     if (chromosome2D[x][col].equals("Q1")) {
                             //System.out.println("YES(Q2)(Row) -->  Q1");
                             nb_conflictsQueen2++;

                 }
                       else if (chromosome2D[x][col].equals("R1")) {
                         //System.out.println("YES(Q2)(Row) -->  R1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][col].equals("R2")) {
                        // System.out.println("YES(Q2)(Row) -->  R2");
                         nb_conflictsQueen2++;
                     }
                        else if (chromosome2D[x][col].equals("B1")) {
                           // System.out.println("YES(Q2)(Row) -->  B1");
                            nb_conflictsQueen2++;
                        } else if (chromosome2D[x][col].equals("B2")) {
                            //  System.out.println("YES(Q2)(Row) -->  B2");
                            nb_conflictsQueen2++;
                        }

                       else if (chromosome2D[x][col].equals("K1")) {
                            // System.out.println("YES(Q2)(Row) -->  K1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][col].equals("K2")) {
                            //  System.out.println("YES(Q2)(Row) -->  K2");
                         nb_conflictsQueen2++;
                     }

                 }

                 // ---- Colmun ----


                 for (int x = 0; x < chromosome2D.length; x++) {

                       if (chromosome2D[x][col].equals("Q1")) {
                       //  System.out.println("YES(Q2)(Colmun) -->  Q1");
                         nb_conflictsQueen2++;
                 }


                      else if (chromosome2D[row][x].equals("R1")) {
                           // System.out.println("YES(Q2)(Colmun) -->  R1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[row][x].equals("R2")) {
                           // System.out.println("YES(Q2)(Colmun) -->  R2");
                         nb_conflictsQueen2++;
                     }
                       else if (chromosome2D[row][x].equals("B1")) {
                           //  System.out.println("YES(Q2)(Colmun) -->  B1");
                           nb_conflictsQueen2++;
                       } else if (chromosome2D[row][x].equals("B2")) {
                           //  System.out.println("YES(Q2)(Colmun) -->  B2");
                           nb_conflictsQueen2++;
                       }
                    else if (chromosome2D[row][x].equals("K1")) {
                           //  System.out.println("YES(Q2)(Colmun) -->  K1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[row][x].equals("K2")) {
                           //  System.out.println("YES(Q2)(Colmun) -->  K2");
                         nb_conflictsQueen2++;
                     }


                 }

                 //  ----- Diagonal ----

                 // upper Left Diagonal

                 for (int x = row, y = col; x >= 0 && y >= 0; x--, y--){
                     if (chromosome2D[row][x].equals("Q1")) {
                         //   System.out.println("YES(Q1)(upper_Left) -->  Q1");
                         nb_conflictsQueen1++;
                     }
                     else if (chromosome2D[x][y].equals("R1")) {
                          // System.out.println("YES(Q1)(upper_Left) -->  R1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("R2")) {
                          // System.out.println("YES(Q1)(upper_Left) -->  R2");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("B1")) {
                          //  System.out.println("YES(Q1)(upper_Left) -->  B1");
                         nb_conflictsQueen2++;
                     } else if (chromosome2D[x][y].equals("B2")) {
                          //  System.out.println("YES(Q1)(upper_Left) -->  B2");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("K1")) {
                          // System.out.println("YES(Q1)(upper_Left) -->  K1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("K2")) {
                          //  System.out.println("YES(Q1)(upper_Left) -->  K2");
                         nb_conflictsQueen2++;
                     }

                 }

                 // lower left Diagonal

                 for (int x = row, y = col; x < chromosome2D.length && y >= 0; x++, y--){
                     if (chromosome2D[row][x].equals("Q1")) {
                            // System.out.println("YES(Q1)(lower_Left) -->  Q1");
                             nb_conflictsQueen2++;
                         }
                        /* else if (chromosome2D[row][x].equals("Q2")) {
                            // System.out.println("YES(Q1)(lower_Left) -->  Q2");
                             nb_conflictsQueen2++;
                         }*/



                     else if (chromosome2D[x][y].equals("R1")) {
                          //  System.out.println("YES(Q1)(lower_Left) -->  R1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("R2")) {
                          //  System.out.println("YES(Q1)(lower_Left) -->  R2");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("B1")) {
                          //  System.out.println("YES(Q1)(lower_Left) -->  B1");
                         nb_conflictsQueen2++;
                     } else if (chromosome2D[x][y].equals("B2")) {
                          //  System.out.println("YES(Q1)(lower_Left) -->  B2");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("K1")) {
                          //  System.out.println("YES(Q1)(lower_Left) -->  K1");
                         nb_conflictsQueen2++;
                     }
                     else if (chromosome2D[x][y].equals("K2")) {
                          //  System.out.println("YES(Q1)(lower_Left) -->  K2");
                         nb_conflictsQueen2++;
                     }

                 }



             }


                else if (chromosome2D[row][col].equals("R1")){
                    // System.out.println("R1("+row+","+col+")");


                     // --- Rook 1 ----

                     // ---- colmun ----
                     for (int x = 0 ; x < chromosome2D.length ; x++){
                          if (chromosome2D[row][x].equals("R2")){
                              //  System.out.println("YES(R1)(Colmun) -->  R2");
                              nb_conflictsRook1++;
                          }
                          else if (chromosome2D[row][x].equals("B1")) {
                              //  System.out.println("YES(R1)(Colmun) -->  B1");
                              nb_conflictsRook1++;
                          }
                          else if (chromosome2D[row][x].equals("B2")) {
                              //  System.out.println("YES(R1)(Colmun) -->  B2");
                              nb_conflictsRook1++;
                          }
                          else if (chromosome2D[row][x].equals("K1")) {
                              //  System.out.println("YES(R1)(Colmun) -->  K1");
                              nb_conflictsRook1++;
                          }
                          else if (chromosome2D[row][x].equals("K2")) {
                              //  System.out.println("YES(R1)(Colmun) -->  K2");
                              nb_conflictsRook1++;
                          }

                     }

                     // ----- Row -----
                     for (int x = 0 ; x < chromosome2D.length ; x++) {
                         if (chromosome2D[x][col].equals("R2")) {
                             // System.out.println("YES(R1)(Row) -->  R2");
                             nb_conflictsRook1++;
                         } else if (chromosome2D[x][col].equals("B1")) {
                             //  System.out.println("YES(R1)(Row) -->  B1");
                             nb_conflictsRook1++;
                         } else if (chromosome2D[x][col].equals("B2")) {
                             //  System.out.println("YES(R1)(Row) -->  B2");
                             nb_conflictsRook1++;
                         } else if (chromosome2D[x][col].equals("K1")) {
                             //  System.out.println("YES(R1)(Row) -->  K1");
                             nb_conflictsRook1++;
                         } else if (chromosome2D[x][col].equals("K2")) {
                             //  System.out.println("YES(R1)(Row) -->  K2");
                             nb_conflictsRook1++;
                         }

                     }


                 }


                 //  ---- Rook 2 ----

               else if (chromosome2D[row][col].equals("R2")){
                            // System.out.println("R2("+row+","+col+")");

                     // ---- colmun ----
                     for (int x = 0 ; x < chromosome2D.length ; x++) {
                         if (chromosome2D[row][x].equals("R1")) {
                             // System.out.println("YES(R2)(Colmun) -->  R1");
                             nb_conflictsRook2++;
                         }

                      else if (chromosome2D[row][x].equals("B1")) {
                          // System.out.println("YES(R2)(Colmun) -->  B1");
                             nb_conflictsRook2++;
                     } else if (chromosome2D[row][x].equals("B2")) {
                          //  System.out.println("YES(R2)(Colmun) -->  B2");
                             nb_conflictsRook2++;
                     }
                       else  if (chromosome2D[row][x].equals("K1")) {
                          //  System.out.println("YES(R2)(Colmun) -->  K1");
                             nb_conflictsRook2++;
                         }
                         else if (chromosome2D[row][x].equals("K2")) {
                          //  System.out.println("YES(R2)(Colmun) -->  K2");
                             nb_conflictsRook2++;
                         }
                     }



                     // ----- Row -----
                     for (int x = 0 ; x < chromosome2D.length ; x++){
                         if (chromosome2D[row][x].equals("R1")) {
                             // System.out.println("YES(R2)(Colmun) -->  R1");
                             nb_conflictsRook2++;
                         }
                       else  if (chromosome2D[x][col].equals("B1")) {
                             // System.out.println("YES(R2)(Row) -->  B1");
                             nb_conflictsRook2++;
                         }
                         else if (chromosome2D[x][col].equals("B2")) {
                             // System.out.println("YES(R2)(Row) -->  B2");
                             nb_conflictsRook2++;
                         }
                        else if (chromosome2D[x][col].equals("K1")) {
                             //  System.out.println("YES(R2)(Row) -->  K1");
                             nb_conflictsRook2++;
                         }
                         else if (chromosome2D[x][col].equals("K2")) {
                             //  System.out.println("YES(R2)(Row) -->  K2");
                             nb_conflictsRook2++;
                         }

                     }

               }





                // --- Bishope ------



                 else if (chromosome2D[row][col].equals("B1")){
                     //System.out.println("B1("+row+","+col+")");

                    // Upper left Diagonal

                     for (int x = row, y = col; x >= 0 && y >= 0; x--, y--){
                          if (chromosome2D[x][y].equals("B2")){
                              // System.out.println("YES(B1)(upper_Left) -->  B2");
                              nb_conflictsBishop1++;
                          }
                         else if (chromosome2D[x][y].equals("K1")){
                              //  System.out.println("YES(B1)(upper_Left) -->  K1");
                             nb_conflictsBishop1++;
                         }
                          else if (chromosome2D[x][y].equals("K2")){
                              // System.out.println("YES(B1)(upper_Left) -->  K2");
                              nb_conflictsBishop1++;
                          }
                     }

                     // lower left Diagonal
                     for (int x = row, y = col; x < chromosome2D.length && y >= 0; x++, y--){
                         if (chromosome2D[x][y].equals("B2")){
                             // System.out.println("YES(B1)(lower_Left) -->  B2");
                             nb_conflictsBishop1++;
                         }
                         else if (chromosome2D[x][y].equals("K1")){
                             // System.out.println("YES(B1)(lower_Left) -->  K1");
                             nb_conflictsBishop1++;
                         }
                         else if (chromosome2D[x][y].equals("K2")){
                             //  System.out.println("YES(B1)(lower_Left) -->  K2");
                             nb_conflictsBishop1++;
                         }
                     }





                 }

                 else if (chromosome2D[row][col].equals("B2")){
                     //System.out.println("B1("+row+","+col+")");

                     // Upper left Diagonal

                     for (int x = row, y = col; x >= 0 && y >= 0; x--, y--){

                         if (chromosome2D[x][y].equals("K1")){
                             //  System.out.println("YES(B2)(upper_Left) -->  K1");
                             nb_conflictsBishop2++;
                         }
                         else if (chromosome2D[x][y].equals("K2")){
                             // System.out.println("YES(B2)(upper_Left) -->  K2");
                             nb_conflictsBishop2++;
                         }
                     }

                     // lower left Diagonal

                     for (int x = row, y = col; x < chromosome2D.length && y >= 0; x++, y--){

                         if (chromosome2D[x][y].equals("K1")){
                             //  System.out.println("YES(B2)(lower_Left) -->  K1");
                             nb_conflictsBishop1++;
                         }
                         else if (chromosome2D[x][y].equals("K2")){
                             //  System.out.println("YES(B2)(lower_Left) -->  K2");
                             nb_conflictsBishop1++;
                         }
                     }
                 }

                // ---- Kinght -----

                 else if (chromosome2D[row][col].equals("K1")){
                 // All possible moves of a knight
                 int X[] = { 2, 1, -1, -2, -2, -1, 1, 2 };
                 int Y[] = { 1, 2, 2, 1, -1, -2, -2, -1 };

                 // Position of knight after move
                 for (int i = 0; i < chromosome2D.length; i++) {
                     int x = row + X[i];
                     int y = col + Y[i];

                     if (x >= 0 && y >= 0 && x < chromosome2D.length && y < chromosome2D.length
                             && chromosome2D[x][y].equals("K2")){
                        // System.out.println("YES(K1) -->  K2");
                         nb_conflictsKinght1++;

                     }

                     else if (x >= 0 && y >= 0 && x < chromosome2D.length && y < chromosome2D.length
                             && chromosome2D[x][y].equals("Q1")){
                         //  System.out.println("YES(K1) -->  Q2");
                         nb_conflictsKinght1++;

                     }

                    else if (x >= 0 && y >= 0 && x < chromosome2D.length && y < chromosome2D.length
                             && chromosome2D[x][y].equals("Q2")){
                       //  System.out.println("YES(K1) -->  Q2");
                         nb_conflictsKinght1++;

                     }
                   else   if (x >= 0 && y >= 0 && x < chromosome2D.length && y < chromosome2D.length
                             && chromosome2D[x][y].equals("B1")){
                        // System.out.println("YES(K1) -->  B1");
                         nb_conflictsKinght1++;

                     }
                     else   if (x >= 0 && y >= 0 && x < chromosome2D.length && y < chromosome2D.length
                             && chromosome2D[x][y].equals("B2")){
                         // System.out.println("YES(K1) -->  B1");
                         nb_conflictsKinght1++;

                     }


                 }


                 }

         }
     }






     // Optimitation Use Penalty

     for (int row = 0; row < chromosome2D.length; row++) {
         int Counter = 0 ;
         for (int col = 0; col < chromosome2D.length; col++) {
            if (chromosome2D[row][col].equals("Q1")||chromosome2D.equals("Q2")){
                Counter++;
            }

         }

         if (Counter > 1){
             penalty++;
         }

     }





     nb_conflictsTotal = (float)(nb_conflictsQueen1+nb_conflictsQueen2+nb_conflictsRook1+nb_conflictsRook2+nb_conflictsBishop1
           +nb_conflictsBishop2+nb_conflictsKinght1+nb_conflictsKinght2
   );
     fitness = (float) ((1/(nb_conflictsTotal+1+penalty))*100) ;

 return fitness; }












































































}


