package com.company;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Main {
    public static void main(String[] args) {

        int[] gang_size = new int[3];
        gang_size[0] = 100000;
        gang_size[1] = 100000;
        gang_size[2] = 0;

        double[] graffiti_rates = new double[3];
        graffiti_rates[0] = 0.5;
        graffiti_rates[1] = 0.5;
        graffiti_rates[2] = 0.5;
        GangModel gm =
                new GangModel(200, 200, 3,  gang_size,
                        3*Math.pow(10,-5)  , 0.55, graffiti_rates, graffiti_rates, "/Users/damlaortac/Desktop/GangTerritory/GangTerritory/applied/chicago.txt");
        //3*Math.pow(10,-5)
        //0.0000065



        int time = 100001;
        GangModel.saveToFile("t", time + "");
        GangModel.saveToFile("dt", "1");
        int number_of_selected_t = 5;
        String selected_t = "";

        double[] order_parameter_values = new double[time];

        for (int i = 0; i < time; i++) {
            boolean flag = false;
              if (i == 1 || i == 10000 || i == 20000 || i == 50000 || i == 100000) {
                  selected_t += i + "\n";
                  flag = true;
              }
//            if (i % (time / number_of_selected_t) == 0) {
//                selected_t +=  (i * time / number_of_selected_t) + "\n";
//                flag = true;
//            }
            gm.calculateOrderParameter();
            if (i == time - 1) {
                gm.calculateOrderParameterForLast();
            }
            gm.produceAndDecayGraffiti(flag);

            gm.randomWalk(flag);
            if (i%1000==0) System.out.println(i);
        }


    //            for (int i = 0; i < time; i++) {
    //                GangModel.saveToFile("order_parameter_averaged", order_parameter_values[i] + "");
    //            }

        GangModel.saveToFile("number_of_selected_t", number_of_selected_t + "");
        GangModel.saveToFile("selected_dt", selected_t);
        GangModel.saveToFile("selected_t", selected_t);

    }
}
