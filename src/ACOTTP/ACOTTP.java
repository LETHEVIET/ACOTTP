/*
modifications:
- instance reader can read TTP files
- Ants.place_ant_TTP is new, so that all ants start in the same city (beneficial/detrimental?)
- local search for ACO: any useful?
*/


package ACOTTP;

import ttp.Optimisation.Optimisation;
import ttp.TTPSolution;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;

import static ttp.Optimisation.Optimisation.bitFlip;

/**
 * ACO algorithms for the TSP
 * <p>
 * This code is based on the ACOTSP project of Thomas Stuetzle.
 * It was initially ported from C to Java by Adrian Wilke.
 * <p>
 * Project website: http://adibaba.github.io/ACOTSPJava/
 * Source code: https://github.com/adibaba/ACOTSPJava/
 */
public class ACOTTP {

    public static boolean boosting = false;
//    public static ArrayList<Ants.ant_struct> logs = new ArrayList<Ants.ant_struct>();
    static ArrayList<Pair<TTPSolution, Long>> logs = new ArrayList<>();
    static long startTime;
    static long minutes_passed = 1;
    static long log_interval = 60 * 1000;
    static long currenttime;

    /*
     * ################################################
     * ########## ACO algorithms for the TSP ##########
     * ################################################
     *
     * Version: 1.0
     * File: main.c
     * Author: Thomas Stuetzle
     * Purpose: main routines and control for the ACO algorithms
     * Check: README and gpl.txt
     * Copyright (C) 2002 Thomas Stuetzle
     */

    /***************************************************************************
     * Program's name: acotsp
     *
     * Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the
     * symmetric TSP
     *
     * Copyright (C) 2004 Thomas Stuetzle
     *
     * This program is free software; you can redistribute it and/or modify
     * it under the terms of the GNU General Public License as published by
     * the Free Software Foundation; either version 2 of the License, or
     * (at your option) any later version.
     *
     * This program is distributed in the hope that it will be useful,
     * but WITHOUT ANY WARRANTY; without even the implied warranty of
     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
     * GNU General Public License for more details.
     *
     * You should have received a copy of the GNU General Public License
     * along with this program; if not, write to the Free Software
     * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
     *
     * email: stuetzle no@spam informatik.tu-darmstadt.de
     * mail address: Universitaet Darmstadt
     * Fachbereich Informatik
     * Hochschulstr. 10
     * D-64283 Darmstadt
     * Germany
     ***************************************************************************/

    static boolean termination_condition()
    /*
     * FUNCTION: checks whether termination condition is met
     * INPUT: none
     * OUTPUT: 0 if condition is not met, number neq 0 otherwise
     * (SIDE)EFFECTS: none
     */ {
        return (((InOut.n_tours >= InOut.max_tours) && (Timer.elapsed_time() >= InOut.max_time)) || (Ants.best_so_far_ant.fitness <= InOut.optimal));
    }

    static void construct_solutions()
        /*
         * FUNCTION: manage the solution construction phase
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: when finished, all ants of the colony have constructed a solution
         */ {
        int k; /* counter variable */
        int step; /* counter of the number of construction steps */

        for (k = 0; k < Ants.n_ants; k++) {
            Ants.ant_empty_memory(Ants.ant[k]);
        }

        step = 0;
        /* Place the ants on same initial city */
        for (k = 0; k < Ants.n_ants; k++)
            Ants.place_ant_TTP(Ants.ant[k], step);

        while (step < TTP.n - 1) {
            step++;
            for (k = 0; k < Ants.n_ants; k++) {
                Ants.neighbour_choose_and_move_to_next(Ants.ant[k], step);
            }
        }

        step = TTP.n;
        for (k = 0; k < Ants.n_ants; k++) {
            Ants.ant[k].tour[TTP.n] = Ants.ant[k].tour[0];
            TTP.compute_fitness(Ants.ant[k]);
        }
        InOut.n_tours += Ants.n_ants;
    }

    static void init_try(int ntry)
        /*
         * FUNCTION: initilialize variables appropriately when starting a trial
         * INPUT: trial number
         * OUTPUT: none
         * COMMENTS: none
         */ {
        Timer.start_timers();
        InOut.time_used = Timer.elapsed_time();
        InOut.time_passed = InOut.time_used;

        if (InOut.comp_report != null) {
            InOut.printToFile(InOut.comp_report, "Utilities.seed " + Utilities.seed);
        }

        /* Initialize variables concerning statistics etc. */

        InOut.n_tours = 1;
        InOut.iteration = 1;
        InOut.restart_iteration = 1;
        InOut.lambda = 0.05;
        Ants.best_so_far_ant.fitness = Integer.MAX_VALUE;
        InOut.found_best = 0;

        if (Ants.mmas_flag) {
            Ants.trail_max = 1. / ((Ants.rho) * Ants.nn_tour());
            Ants.trail_min = Ants.trail_max / (2. * TTP.n);
            Ants.init_pheromone_trails(Ants.trail_max);
        }

        /* Calculate combined information Ants.pheromone times heuristic information */
        Ants.compute_total_information();

        if (InOut.comp_report != null)
            InOut.printToFile(InOut.comp_report, "begin try " + ntry);
        if (InOut.stat_report != null)
            InOut.printToFile(InOut.stat_report, "begin try " + ntry);
    }

    static void local_search()
        /*
         * FUNCTION: manage the local search phase; apply local search to ALL ants; in
         * dependence of LocalSearch.ls_flag one of 2-opt, 2.5-opt, and 3-opt local search
         * is chosen.
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: all ants of the colony have locally optimal tours
         * COMMENTS: typically, best performance is obtained by applying local search
         * to all ants. It is known that some improvements (e.g. convergence
         * speed towards high quality solutions) may be obtained for some
         * ACO algorithms by applying local search to only some of the ants.
         * Overall best performance is typcially obtained by using 3-opt.
         */ {
        int k;

        // TRACE ( System.out.println("apply local search to all ants\n"); );

        for (k = 0; k < Ants.n_ants; k++) {
            switch (LocalSearch.ls_flag) {
                case 0: // TTP: don't do anything, needs to be considered in the switch statement
                    break;
                case 1:
                    LocalSearch.two_opt_first(Ants.ant[k].tour); /* 2-opt local search */
                    break;
                case 2:
                    LocalSearch.two_h_opt_first(Ants.ant[k].tour); /* 2.5-opt local search */
                    break;
                case 3:
                    LocalSearch.three_opt_first(Ants.ant[k].tour); /* 3-opt local search */
                    break;
                case 4:                                                             //TTP: randomly change between them
                    switch ((int) (Utilities.ran01(Utilities.seed) * 3 + 1)) {
                        case 1:
                            LocalSearch.two_opt_first(Ants.ant[k].tour); /* 2-opt local search */
                            break;
                        case 2:
                            LocalSearch.two_h_opt_first(Ants.ant[k].tour); /* 2.5-opt local search */
                            break;
                        case 3:
                            LocalSearch.three_opt_first(Ants.ant[k].tour); /* 3-opt local search */
                            break;
                        default:
                    }
                    break;

                default:
                    System.err.println("type of local search procedure not correctly specified");
                    System.exit(1);
            }
            if (termination_condition())
                return;
            TTP.compute_fitness(Ants.ant[k]);
        }
    }

    static void update_statistics()
        /*
         * FUNCTION: manage some statistical information about the trial, especially
         * if a new best solution (best-so-far or restart-best) is found and
         * adjust some parameters if a new best solution is found
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: restart-best and best-so-far ant may be updated; Ants.trail_min
         * and Ants.trail_max used by MMAS may be updated
         */ {

        int iteration_best_ant;
        double p_x; /* only used by MMAS */

        iteration_best_ant = Ants.find_best(); /* iteration_best_ant is a global variable */

        // try to boost the best:
        double obOld = Ants.ant[iteration_best_ant].fitness;

//        TTPSolution sol = null;

        //TTP: break added to account for slow local searches
        if ((Timer.elapsed_time() >= InOut.max_time - 5)) {

        } else if (obOld < 0.95 * Ants.best_so_far_ant.fitness && boosting) { // attempt boosting
//            System.out.print("boosting...");
            int[] tlocal = TTP.rotateTour(Ants.ant[iteration_best_ant].tour);

            Ants.ant[iteration_best_ant].sol = Optimisation.hillClimber(
                    TTP.currentTTPInstance,
                    Ants.ant[iteration_best_ant].sol.tspTour,
                    Ants.ant[iteration_best_ant].sol.packingPlan,
                    2,
                    10000,
                    60000); // 1+1EA
            Ants.ant[iteration_best_ant].sol = Optimisation.insertionReverse(
                    TTP.currentTTPInstance,
                    Ants.ant[iteration_best_ant].sol.tspTour,
                    Ants.ant[iteration_best_ant].sol.packingPlan,
                    60000,
                    1);
            Ants.ant[iteration_best_ant].sol = bitFlip(
                    TTP.currentTTPInstance,
                    Ants.ant[iteration_best_ant].sol,
                    60000,
                    1); // bitflip

            Ants.ant[iteration_best_ant].sol.computationTime = (long) (Timer.elapsed_time() * 1000);

            TTP.currentTTPInstance.evaluate(Ants.ant[iteration_best_ant].sol, false);
            double obNew = Ants.ant[iteration_best_ant].sol.getObjective() * -1;
            if (obNew < obOld) {
                Ants.ant[iteration_best_ant].fitness = obNew;
                // update hashmap?
                TTP.tours.put(Arrays.toString(tlocal), Ants.ant[iteration_best_ant].sol);
            }
        }

        if (Ants.ant[iteration_best_ant].fitness < Ants.best_so_far_ant.fitness) {

            InOut.time_used = Timer.elapsed_time(); /* best sol found after time_used */
            Ants.copy_from_to(Ants.ant[iteration_best_ant], Ants.best_so_far_ant);
            Ants.copy_from_to(Ants.ant[iteration_best_ant], Ants.restart_best_ant);

            InOut.found_best = InOut.iteration;
            InOut.restart_found_best = InOut.iteration;
            InOut.found_branching = InOut.node_branching(InOut.lambda);
            InOut.branching_factor = InOut.found_branching;
            if (Ants.mmas_flag) {
                if (LocalSearch.ls_flag == 0) {
                    p_x = Math.exp(Math.log(0.05) / TTP.n);
                    Ants.trail_min = (1. - p_x) / (p_x * (double) ((Ants.nn_ants + 1) / 2));
                    Ants.trail_max = 1. / ((Ants.rho) * Ants.best_so_far_ant.fitness);
                    Ants.trail_0 = Ants.trail_max;
                    Ants.trail_min = Ants.trail_max * Ants.trail_min;
                } else {
                    Ants.trail_max = 1. / ((Ants.rho) * Ants.best_so_far_ant.fitness);
                    Ants.trail_min = Ants.trail_max / (2. * TTP.n);
                    Ants.trail_0 = Ants.trail_max;
                }
            }
            InOut.write_report();
        }
        if (Ants.ant[iteration_best_ant].fitness < Ants.restart_best_ant.fitness) {
            Ants.copy_from_to(Ants.ant[iteration_best_ant], Ants.restart_best_ant);
            InOut.restart_found_best = InOut.iteration;
            System.out.println("restart best: " + Ants.restart_best_ant.fitness + " restart_found_best "
                    + InOut.restart_found_best + ", time " + Timer.elapsed_time());
        }
    }

    static void search_control_and_statistics()
        /*
         * FUNCTION: occasionally compute some statistics and check whether algorithm
         * is converged
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: restart-best and best-so-far ant may be updated; Ants.trail_min
         * and Ants.trail_max used by MMAS may be updated
         */ {
        // TRACE ( System.out.println("SEARCH CONTROL AND STATISTICS\n"); );

        if ((InOut.iteration % 100) == 0) {
            InOut.population_statistics();
            InOut.branching_factor = InOut.node_branching(InOut.lambda);
            System.out.println("best so far " + Ants.best_so_far_ant.fitness + ", iteration: " + InOut.iteration
                    + ", time " + Timer.elapsed_time() + ", b_fac " + InOut.branching_factor);

            if (Ants.mmas_flag && (InOut.branching_factor < InOut.branch_fac)
                    && (InOut.iteration - InOut.restart_found_best > 250)) {
                /*
                 * MAX-MIN Ant System was the first ACO algorithm to use
                 * Ants.pheromone trail re-initialisation as implemented
                 * here. Other ACO algorithms may also profit from this mechanism.
                 */
                System.out.println("INIT TRAILS!!!\n");
                Ants.restart_best_ant.fitness = Integer.MAX_VALUE;
                Ants.init_pheromone_trails(Ants.trail_max);
                Ants.compute_total_information();
                InOut.restart_iteration = InOut.iteration;
                InOut.restart_time = Timer.elapsed_time();
            }
            System.out.println("try " + InOut.n_try + " iteration " + InOut.iteration + ", b-fac "
                    + InOut.branching_factor);
        }
    }

    static void mmas_update()
        /*
         * FUNCTION: manage global Ants.pheromone deposit for MAX-MIN Ant System
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: either the iteration-best or the best-so-far ant deposit Ants.pheromone
         * on matrix "Ants.pheromone"
         */ {
        /*
         * we use default upper Ants.pheromone trail limit for MMAS and hence we
         * do not have to worry regarding keeping the upper limit
         */

        int iteration_best_ant;

        // TRACE ( System.out.println("MAX-MIN Ant System Ants.pheromone deposit\n"); );

        if (InOut.iteration % Ants.u_gb == 0) {
            iteration_best_ant = Ants.find_best();
            Ants.global_update_pheromone(Ants.ant[iteration_best_ant]);
        } else {
            if (Ants.u_gb == 1 && (InOut.iteration - InOut.restart_found_best > 50))
                Ants.global_update_pheromone(Ants.best_so_far_ant);
            else
                Ants.global_update_pheromone(Ants.restart_best_ant);
        }

        if (LocalSearch.ls_flag != 0) {
            /*
             * implement the schedule for Ants.u_gb as defined in the
             * Future Generation Computer Systems article or in Stuetzle's PhD thesis.
             * This schedule is only applied if local search is used.
             */
            if ((InOut.iteration - InOut.restart_iteration) < 25)
                Ants.u_gb = 25;
            else if ((InOut.iteration - InOut.restart_iteration) < 75)
                Ants.u_gb = 5;
            else if ((InOut.iteration - InOut.restart_iteration) < 125)
                Ants.u_gb = 3;
            else if ((InOut.iteration - InOut.restart_iteration) < 250)
                Ants.u_gb = 2;
            else
                Ants.u_gb = 1;
        } else
            Ants.u_gb = 25;

    }

    static void pheromone_trail_update()
        /*
         * FUNCTION: manage global Ants.pheromone trail update for the ACO algorithms
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: Ants.pheromone trails are evaporated and Ants.pheromones are deposited
         * according to the rules defined by the various ACO algorithms.
         */ {

        Ants.mmas_evaporation_nn_list();

        mmas_update();

        Ants.compute_nn_list_total_information();
    }

    /* --- main program ------------------------------------------------------ */

    public static void main(String[] args) {

        if (args.length == 0) {
            args = new String[]{
                    "--rho",
                    "0.5",
                    "--alpha",
                    "1",
                    "--beta",
                    "2",
                    "--ants",
                    "25",
                    "--time",
                    "45",
                    "--tours",
                    "100",
                    "--tries",
                    "1",
                    "--elitistants",
                    "100",
                    "--rasranks",
                    "6",
                    "--localsearch",
                    "4",                       // 0:no local search  1:2-opt  2:2.5-opt  3:3-opt   4:randomly switch between 1..3 for local search
//            "--boosting",
                    "-seed",
                    "321",
                    "-i",
                    "instances/a280_n279_bounded-strongly-corr_01.ttp"};
        }

        /*
         * FUNCTION: main control for running the ACO algorithms
         * INPUT: none
         * OUTPUT: none
         * (SIDE)EFFECTS: none
         * COMMENTS: this function controls the run of "max_tries" independent trials
         */
//	for (String argument : args) {
//	    System.out.println(argument);
//	}
        Timer.start_timers();

        InOut.init_program(args);

        Optimisation.setRandomNumberSeed(Utilities.seed); // needed to make the TTP computations deterministic

        TTP.instance.nn_list = TTP.compute_nn_lists();
        Ants.pheromone = Utilities.generate_double_matrix(TTP.n, TTP.n);
        Ants.total = Utilities.generate_double_matrix(TTP.n, TTP.n);


        /* TTP INIT START */
        boolean debugPrint = false;
        File ttpFile = new File(args[args.length - 1]);
        if (debugPrint) System.out.println(new File(".").getAbsolutePath());
        if (debugPrint)
            System.out.println("readInstanceFiles: " + ttpFile.getAbsolutePath() + " exists=" + ttpFile.exists());
        TTP.currentTTPInstance = new ttp.TTPInstance(ttpFile);
        /* TTP END */

        InOut.time_used = Timer.elapsed_time();
        System.out.println("Initialization took " + InOut.time_used + " seconds\n");

        init_try(InOut.n_try);

        startTime = System.currentTimeMillis();
        logs.add(new Pair<>(Ants.ant[0].sol, Long.parseLong(0+"")));

        while (!termination_condition()) {
            if (debugPrint) System.out.println(InOut.iteration + " start");

            construct_solutions();

            if (LocalSearch.ls_flag > 0)
                local_search();

            update_statistics();

            pheromone_trail_update();

            search_control_and_statistics();

            InOut.iteration++;
            if (debugPrint) System.out.println(InOut.iteration + " end");

            currenttime = System.currentTimeMillis();
            if ((currenttime - startTime) / log_interval >= minutes_passed){
                logs.add(new Pair<>(Ants.best_so_far_ant.sol, (currenttime - startTime)));
                minutes_passed += (currenttime - startTime) / log_interval - minutes_passed + 1;
            }
        }
        currenttime = System.currentTimeMillis();
        logs.add(new Pair<>(Ants.best_so_far_ant.sol, (currenttime - startTime)));

        InOut.exit_try(InOut.n_try);

        InOut.exit_program();

        // Added by AW
        double aw_best_tour_length = Utilities.best_of_vector(InOut.best_in_try, InOut.max_tries);
        String aw_best_tour = InOut.aw_best_tour_in_try[Utilities.aw_best_tour_index()];
        try {
            Writer w = new OutputStreamWriter(new FileOutputStream(InOut.solution_dir + "/tour." + TTP.instance.name), StandardCharsets.UTF_8);
            BufferedWriter out = new BufferedWriter(w);
            out.write(aw_best_tour_length + "\n");
            out.write(aw_best_tour);
            out.close();
        } catch (IOException e) {
            System.err.print("Could not write file tour." + TTP.instance.name + " " + e.getMessage());
            System.exit(1);
        }

        System.out.println();
        System.out.println("Best tour:");
        System.out.println(aw_best_tour_length);
        System.out.println(aw_best_tour);

        System.out.print(ttpFile.getName() + ": ");
        Ants.best_so_far_ant.sol.println();

        InOut.writelog(logs);
    }
}
