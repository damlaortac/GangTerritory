package com.company;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Main {

    public static void main(String[] args) {
        int[] gang_size = new int[3];
        gang_size[0] = 50000;
        gang_size[1] = 50000;
        gang_size[2] = 0;

        double[] graffiti_rates = new double[3];
        graffiti_rates[0] = 0.5;
        graffiti_rates[1] = 0.5;
        graffiti_rates[2] = 0.5;
        GangModel gm =
                new GangModel(100, 100, 3,  gang_size,
                        3*Math.pow(10,-5)  , 4, graffiti_rates, graffiti_rates, "/Users/damlaortac/Desktop/GangTerritory/GangTerritory/100x100/100x100_r.txt");
        //3*Math.pow(10,-5)
        //0.0000065
        int time = 100000;
        GangModel.saveToFile("t", time + "");
        GangModel.saveToFile("dt", "1");
        int number_of_selected_t = 10;
        String selected_t = "";



        for (int i = 0; i < time; i++) {
            boolean flag = false;
            if (i % (time / number_of_selected_t) == 0) {
                selected_t +=  (i * time / number_of_selected_t) + "\n";
                flag = true;
            }
            gm.calculateOrderParameter();
            gm.produceAndDecayGraffiti(flag);
            gm.randomWalk(flag);
            if (i%1000==0) System.out.println(i);
        }

        GangModel.saveToFile("number_of_selected_t", number_of_selected_t + "");
        GangModel.saveToFile("selected_dt", selected_t);
        GangModel.saveToFile("selected_t", selected_t);
    }
}
