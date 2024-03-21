package com.example.chessboard;

import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.stage.Stage;

import java.io.IOException;

public class HelloApplication extends Application {
    @Override
    public void start(Stage stage) throws IOException {
        try {

            /*Parent root = FXMLLoader.load(getClass().getResource("/com/example/chessboard/chessboard.fxml"));
            Scene scene = new Scene(root);
            String css = this.getClass().getResource("chessboard.css").toExternalForm();


            //stage.initStyle(StageStyle.UNDECORATED);
            stage.setScene(scene);
            scene.getStylesheets().add(css);
            stage.show();*/

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm();
        //geneticAlgorithm.returnChessBoard();
         //Chromosomes chromosomes = geneticAlgorithm.initializingMatrice();
        geneticAlgorithm.createPopulation();



// launch();


    }
}