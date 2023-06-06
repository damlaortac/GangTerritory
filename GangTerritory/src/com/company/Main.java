package com.company;

import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
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

        ExecutorService srv = Executors.newFixedThreadPool(20);

        for (int j = 0; j < 50; j++) {
            GangModel gm =
                    new GangModel(100, 100, 3,  gang_size,
                            3*Math.pow(10,-5), 0, graffiti_rates, graffiti_rates, "/Users/damlaortac/Desktop/GangTerritory/GangTerritory/100x100/100x100_r.txt");
            srv.execute(gm);

            //0, 1, 2, 5, 15
        }
        srv.shutdown();
        try {
            srv.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }

    }
}
