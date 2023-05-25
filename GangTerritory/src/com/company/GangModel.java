package com.company;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;
import java.util.concurrent.ThreadLocalRandom;

public class GangModel {
    private int[][][] agent_lattice;
    private double[][][] graffiti_lattice;
    private Pixel[][] city_map;
    private double[][] hardness;
    private int area;

    private int lattice_row;
    private int lattice_column;
    private int N_gang;
    private int[] agent_sizes;
    private double beta;
    private double alpha;
    private double[] graffiti_rates;
    private double[] graffiti_decay_rates;


    public GangModel(int lattice_row, int lattice_column, int N_gang, int[] agent_sizes, double beta, double alpha, double[] graffiti_rates, double[] graffiti_decay_rates, String city_map_path) {
        this.lattice_row = lattice_row;
        this.lattice_column = lattice_column;
        this.N_gang = N_gang;
        this.agent_sizes = agent_sizes;
        this.beta = beta;
        this.alpha = alpha;
        this.graffiti_rates = graffiti_rates;
        this.graffiti_decay_rates = graffiti_decay_rates;

        //this.city_map = readCityMapFromCSV(city_map_path);
        this.hardness = calculateHardnessFromCityMap();

        this.area = lattice_row * lattice_column;
        for (int i = 0; i < lattice_row; i++) {
            for (int j = 0; j < lattice_column; j++) {
                if (hardness[i][j] == 1.0) area--;
            }
        }
        System.out.println(area);

        agent_lattice = new int[lattice_row][lattice_column][N_gang];
        graffiti_lattice = new double[lattice_row][lattice_column][N_gang];

        for (int i = 0; i < N_gang; i++) {
            int current_agent_size = agent_sizes[i];
            for (int j = 0; j < current_agent_size; j++) {
                int x = ThreadLocalRandom.current().nextInt(0, lattice_row);
                int y = ThreadLocalRandom.current().nextInt(0, lattice_column);
                while (hardness[x][y] == 1.0) {
                    x = ThreadLocalRandom.current().nextInt(0, lattice_row);
                    y = ThreadLocalRandom.current().nextInt(0, lattice_column);
                }
                agent_lattice[x][y][i]++;
            }
        }
        saveToFile("row", lattice_row + "");
        saveToFile("column", lattice_column + "");
        saveToFile("beta", beta + "");
        for (int k = 0; k < N_gang; k++) {
            saveToFile("N_" + k, agent_sizes[k] + "");
        }
    }

    public static void saveToFile(String file_name, String content) {
        File file = new File("Output/" + file_name + ".txt");
        try {
            file.createNewFile();
            FileWriter fos = new FileWriter(file);
            fos.append(content);
            fos.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Pixel[][] readCityMapFromCSV(String file_path) {
        File file = new File(file_path);
        Scanner sc = null;
        try {
            sc = new Scanner(file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        sc.useDelimiter("[\\s,]+");
        int row = Integer.parseInt(sc.next());
        int column = Integer.parseInt(sc.next());

        Pixel[][] pixels = new Pixel[row][column];

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                pixels[i][j] = new Pixel(Integer.parseInt(sc.next()), Integer.parseInt(sc.next()), Integer.parseInt(sc.next()));
            }
        }
        return pixels;
    }

    public double[][] calculateHardnessFromCityMap() {
        double[][] hardness_from_city_map = new double[lattice_row][lattice_column];
        return hardness_from_city_map;
//        for (int i = 0; i < lattice_row; i++) {
//            for (int j = 0; j < lattice_column; j++) {
//                Pixel current_pixel = city_map[i][j];
//                if (current_pixel.getR() == 255 && current_pixel.getG() == 255 && current_pixel.getB() == 255) {
//                    hardness_from_city_map[i][j] = 0.0;
//                }
//                else hardness_from_city_map[i][j] = 1.0;
//            }
//        }
//        return hardness_from_city_map;
    }

    public void produceAndDecayGraffiti(boolean save) {
        if (save) saveGraffitiToFile();
        for (int i = 0; i < lattice_row; i++) {
            for (int j = 0; j < lattice_column; j++) {
                for (int k = 0; k < N_gang; k++) {
                    graffiti_lattice[i][j][k] = graffiti_lattice[i][j][k] * (1 - graffiti_decay_rates[k]);
                    //double rnd = ThreadLocalRandom.current().nextDouble(0, 1);
                    //if (rnd <= graffiti_rates[k]) graffiti_lattice[i][j][k] += agent_lattice[i][j][k];

                    for (int count = 0; count < agent_lattice[i][j][k]; count++) {
                        double rnd = ThreadLocalRandom.current().nextDouble(0, 1);
                        if (rnd <= graffiti_rates[k]) graffiti_lattice[i][j][k]++;
                    }
                }

            }
        }
    }

    private void saveGraffitiToFile() {
        for (int k = 0; k < N_gang; k++) {
            File file = new File("Output/G_" + k + ".txt");
            try {
                file.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }

            try {
                FileWriter fos = new FileWriter(file, true);
                for (int i = 0; i < lattice_row; i++) {
                    for (int j = 0; j < lattice_column; j++) {
                        fos.append(graffiti_lattice[i][j][k] + " ");
                    }
                    fos.append("\n");
                }
                fos.append("\n");
                fos.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public void randomWalk(boolean save) {
        if (save) saveAgentToFile();
        int[][][] new_agent_lattice = new int[lattice_row][lattice_column][N_gang];

        for (int i = 0; i < lattice_row; i++) {
            for (int j = 0; j < lattice_column; j++) {
                double[] sum_of_opposing_gang_graffiti_densities_up =
                        isAllowed(i-1, j)? sumOfOpposingGangGraffitiDensities(i-1, j) : new double[N_gang];
                double[] sum_of_opposing_gang_graffiti_densities_down =
                        isAllowed(i+1, j)? sumOfOpposingGangGraffitiDensities(i+1, j) : new double[N_gang];
                double[] sum_of_opposing_gang_graffiti_densities_left =
                        isAllowed(i, j-1)? sumOfOpposingGangGraffitiDensities(i, j-1) : new double[N_gang];
                double[] sum_of_opposing_gang_graffiti_densities_right =
                        isAllowed(i, j+1)? sumOfOpposingGangGraffitiDensities(i, j+1) : new double[N_gang];

                for (int k = 0; k  < N_gang; k++) {
                    double[][][] temp = graffiti_lattice;
                    double factor_up = isAllowed(i-1, j)?
                            Math.exp(((-1) * beta * area * sum_of_opposing_gang_graffiti_densities_up[k])
                                    //- alpha * hardness[i-1][j]
                            ) : 0;

                    double factor_down = isAllowed(i+1, j)?
                            Math.exp(((-1) * beta * area * sum_of_opposing_gang_graffiti_densities_down[k])
                                    //- alpha * hardness[i+1][j]
                            ) : 0;

                    double factor_left = isAllowed(i, j-1)?
                            Math.exp(((-1) * beta * area * sum_of_opposing_gang_graffiti_densities_left[k])
                                    //- alpha * hardness[i][j-1]
                            ) : 0;

                    double factor_right = isAllowed(i, j+1)?
                            Math.exp(((-1) * beta * area * sum_of_opposing_gang_graffiti_densities_right[k])
                                    //- alpha * hardness[i][j+1]
                            ) : 0;

                    double sum_of_factors = factor_up + factor_down + factor_left + factor_right;

                    if (sum_of_factors == 0) {
                        System.out.println("AGENT CAN'T MOVE ANYWHERE???");
                        continue;
                    }
                    double probability_up = factor_up / sum_of_factors;
                    double probability_down = factor_down / sum_of_factors;
                    double probability_left = factor_left / sum_of_factors;
                    double probability_right = factor_right / sum_of_factors;


                    for (int agent = 0; agent < agent_lattice[i][j][k]; agent++) {
                        double rnd = ThreadLocalRandom.current().nextDouble(0, 1);
                        if (rnd < probability_up) {
                            new_agent_lattice[i-1][j][k]++;
                        }
                        else if (rnd <= probability_up + probability_down) {
                            new_agent_lattice[i+1][j][k]++;
                        }
                        else if (rnd <= probability_up + probability_down + probability_left) {
                            new_agent_lattice[i][j-1][k]++;
                        }
                        else if (rnd <= probability_up + probability_down + probability_left + probability_right) {
                            new_agent_lattice[i][j+1][k]++;
                        }
                        else {
                            System.out.println("ERROR WHERE AGENT CAN'T MOVE");
                        }
                    }
                }


            }
        }

        agent_lattice = new_agent_lattice;
    }

    private double[] sumOfOpposingGangGraffitiDensities(int x, int y) {
        double[] sum_of_opposing_gang_graffiti_densities = new double[N_gang];
        for (int k = 0; k  < N_gang; k++) {
            for (int opposing_gang = 0; opposing_gang  < N_gang; opposing_gang++) {
                if (k == opposing_gang) continue;
                sum_of_opposing_gang_graffiti_densities[k] += graffiti_lattice[x][y][opposing_gang]; /// area;
            }
        }
        return sum_of_opposing_gang_graffiti_densities;
    }
















    public void randomWalkOld(boolean save) {
        if (save) saveAgentToFile();
        for (int i = 0; i < lattice_row; i++) {
            for (int j = 0; j < lattice_column; j++) {
                double[] ePowerBetaDensityUp = new double[N_gang];
                double[] ePowerBetaDensityDown = new double[N_gang];
                double[] ePowerBetaDensityLeft = new double[N_gang];
                double[] ePowerBetaDensityRight = new double[N_gang];
                for (int k = 0; k < N_gang; k++) {
                    ePowerBetaDensityUp[k] = Math.exp( -1 * beta *
                            (isAllowed(i, j - 1) ? graffiti_lattice[i][j - 1][k] * area : 0));
                    ePowerBetaDensityDown[k] = Math.exp( -1 * beta *
                            (isAllowed(i, j + 1) ? graffiti_lattice[i][j + 1][k] * area : 0));
                    ePowerBetaDensityLeft[k] = Math.exp( -1 * beta *
                            (isAllowed(i - 1, j) ? graffiti_lattice[i - 1][j][k] * area : 0));
                    ePowerBetaDensityRight[k] = Math.exp( -1 * beta *
                            (isAllowed(i + 1, j) ? graffiti_lattice[i + 1][j][k] * area : 0));
                }
                for (int k = 0; k < N_gang; k++) {
                    double sum_of_opposing_gang_graffitti_up = 0;
                    double sum_of_opposing_gang_graffitti_down = 0;
                    double sum_of_opposing_gang_graffitti_left = 0;
                    double sum_of_opposing_gang_graffitti_right = 0;
                    for (int l = 0; l < N_gang; l++) {
                        if (k == l) continue;
                        sum_of_opposing_gang_graffitti_up += ePowerBetaDensityUp[l];
                        sum_of_opposing_gang_graffitti_down += ePowerBetaDensityDown[l];
                        sum_of_opposing_gang_graffitti_left += ePowerBetaDensityLeft[l];
                        sum_of_opposing_gang_graffitti_right += ePowerBetaDensityRight[l];
                    }
                    double sum_of_opposing_gang_graffitti_all =
                            sum_of_opposing_gang_graffitti_up +
                            sum_of_opposing_gang_graffitti_down +
                            sum_of_opposing_gang_graffitti_left +
                            sum_of_opposing_gang_graffitti_right;

                    for (int count = 0; count < agent_lattice[i][j][k]; count++) {
                        double probabilityUp = isAllowed(i, j-1)? sum_of_opposing_gang_graffitti_up / sum_of_opposing_gang_graffitti_all : 0;
                        double probabilityDown = isAllowed(i, j+1)? sum_of_opposing_gang_graffitti_down / sum_of_opposing_gang_graffitti_all : 0;
                        double probabilityLeft = isAllowed(i-1, j)? sum_of_opposing_gang_graffitti_left / sum_of_opposing_gang_graffitti_all : 0;
                        double probabilityRight = isAllowed(i+1, j)? sum_of_opposing_gang_graffitti_right / sum_of_opposing_gang_graffitti_all : 0;
                        double rnd = 0;
                        if (probabilityUp + probabilityDown + probabilityLeft + probabilityRight >= 0) {
                            rnd = ThreadLocalRandom.current().nextDouble(0, probabilityUp + probabilityDown + probabilityLeft + probabilityRight);
                        }

                        if (rnd < probabilityUp) {
                            agent_lattice[i][j][k]--;
                            agent_lattice[i][j - 1][k]++;
                        }
                        else if (rnd < probabilityUp + probabilityDown) {
                            agent_lattice[i][j][k]--;
                            agent_lattice[i][j + 1][k]++;
                        }
                        else if (rnd < probabilityUp + probabilityDown + probabilityLeft) {
                            agent_lattice[i][j][k]--;
                            agent_lattice[i - 1][j][k]++;
                        }
                        else if (rnd < probabilityUp + probabilityDown + probabilityLeft + probabilityRight) {
                            agent_lattice[i][j][k]--;
                            agent_lattice[i + 1][j][k]++;
                        }
                        else {
                            System.out.println("ERROR WHERE AGENT CAN'T MOVE");
                        }
                    }
                }

            }
        }

    }

    private void saveAgentToFile() {
        for (int k = 0; k < N_gang; k++) {
            File file = new File("Output/A_" + k + ".txt");
            try {
                file.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }

            try {
                FileWriter fos = new FileWriter(file, true);
                for (int i = 0; i < lattice_row; i++) {
                    for (int j = 0; j < lattice_column; j++) {
                        fos.append(agent_lattice[i][j][k] + " ");
                    }
                    fos.append("\n");
                }
                fos.append("\n");
                fos.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public boolean isAllowed(int x, int y) {
        if (x < 0 || y < 0 || x >= lattice_row || y >= lattice_column) return false;
        else if (hardness[x][y] == 1.0) return false;
        else return true;
    }

    public void calculateHamiltonian() {
        for (int i = 0; i < lattice_row; i++) {
            for (int j = 0; j < lattice_column; j++) {
                ThreadLocalRandom.current().nextInt(0, lattice_row);
            }
        }

    }
}
