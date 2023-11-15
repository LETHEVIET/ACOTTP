/**
 * ACO algorithms for the TSP
 * Usage: ACO algorithms for the TSP [OPTION]... [ARGUMENT]...
 * Options:
 * -seed <arg>              # seed for the random number generator
 * -quiet                   reduce output to a minimum, no extra files
 * -a,--alpha <arg>         # alpha (influence of pheromone trails)
 * -b,--beta <arg>          # beta (influence of heuristic information)
 * -c,--elitistants <arg>   # number of elitist ants
 * -d,--dlb                 1 use don't look bits in local search
 * -e,--rho <arg>           # rho: pheromone trail evaporation
 * -f,--rasranks <arg>      # number of ranks in rank-based Ant System
 * -g,--nnants <arg>        # nearest neighbours in tour construction
 * -h,--help                display this help text and exit
 * -i,--tsplibfile <arg>    f inputfile (TSPLIB format necessary)
 * -k,--nnls <arg>          # No. of nearest neighbors for local search
 * -l,--localsearch <arg>   0:no local search  1:2-opt  2:2.5-opt  3:3-opt
 * -m,--ants <arg>          # number of ants
 * -o,--optimum <arg>       # stop if tour better or equal optimum is found
 * -q,--q0 <arg>            # q_0: prob. of best choice in tour construction
 * -r,--tries <arg>         # number of independent trials
 * -s,--tours <arg>         # number of steps in each trial
 * -t,--time <arg>          # maximum time for each trial
 * -u,--as                  apply basic Ant System
 * -v,--eas                 apply elitist Ant System
 * -w,--ras                 apply rank-based version of Ant System
 * -x,--mmas                apply MAX-MIN ant_colony system
 * -y,--bwas                apply best-worst ant_colony system
 * -z,--acs                 apply ant_colony colony system
 */

package de.adrianwilke.acotspjava.test;

import java.util.LinkedList;
import java.util.List;

import de.adrianwilke.acotspjava.AcoTsp;

public class Test {

    private static final String TEST_FILE =
//            "instances/d1291.tsp"
//            "instances/eil51.tsp"
            "instances/a280_n279_bounded-strongly-corr_01.ttp"
//            "instances/eil51_n50_bounded-strongly-corr_01.ttp" // 3938.79 in single run (aco better!)
//            "instances/eil51_n500_uncorr_10.ttp" // 64941.2 in single run (aco better!)
//              "instances/dsj1000_n9990_bounded-strongly-corr_07.ttp"
//            "instances/eil101_n100_bounded-strongly-corr_01.ttp" // 5000.1787 in a single run (aco better!)
//            "instances/eil101_n1000_uncorr_10.ttp" // 129291 in a single run (aco is better!)
//            "instances/rat783.tsp"
//            "instances/rat195_n1940_bounded-strongly-corr_03.ttp"                     // great average: 227965.3565 (ACO seems close)
//            "instances/rat195_n1940_bounded-strongly-corr_07.ttp"                     // great average: 359672.7401 (ACO seems close)
//            "instances/rat195_n582_bounded-strongly-corr_07.ttp"                      // great average: 110555.6726 (ACO seems close)
//            "instances/rat783_n2346_bounded-strongly-corr_03.ttp"                     // gread average: 263725.5363
//            "instances/rat783_n7820_uncorr-similar-weights_07.ttp"                    // great average: 699101.0504
//            "instances/pcb3038_n9111_uncorr-similar-weights_07.ttp"                   // great average: 873670.2233
            ;

    private static List<String> argList = new LinkedList<String>();

    public static void main(String[] args) {
        System.out.println(TEST_FILE);
        // System.out.println(new File("").getAbsolutePath());

//	argList.add("-h");

        // putTestParameters();
        // argList.add("-l");
        // argList.add("0");

        // ./AcoTsp --ants 20 --time 2 --tours 2 --tries 2 -i ../../d1291.tsp -l 0 --as

        // Comments: ACOTSP, ACOTSPJava

        // argList.add("--as");
        // -l 0: 62571,67053
        // -l 1: 52842,53925
        // -l 2: 52547,53464
        // -l 3: 51759,52609

        // argList.add("--eas");
        // -l 0: 84079,90539
        // -l 3: 51524,52335

        // argList.add("--ras");
        // -l 0: 68201,78482
        // -l 3: 51401,51890

        // argList.add("--mmas");
        // -l 0: 97732,98702
        // -l 3: 52179,52698

        // argList.add("--bwas");
        // -l 0:92378,94083
        // -l 3:51076,52068

        // argList.add("--acs");
        // -l 0:61881,63012
        // -l 3:50962,51401

        putTestParameters();

        AcoTsp.main(argList.toArray(new String[0]));

    }

    private static void putTestParameters() {
//        argList.add("--acs");

        argList.add("--rho");
        argList.add("0.5");

        argList.add("--alpha");
        argList.add("1");

        argList.add("--beta");
        argList.add("2");

        argList.add("--ants");
        argList.add("25");//25

        argList.add("--time");
        argList.add("600");

        argList.add("--tours");
        argList.add("100");

        argList.add("--tries");
        argList.add("1");

        argList.add("--elitistants");
        argList.add("100");

        argList.add("--rasranks");
        argList.add("6");

        argList.add("--localsearch");
        argList.add("3");                       // 0:no local search  1:2-opt  2:2.5-opt  3:3-opt

        argList.add("-i");
        argList.add(TEST_FILE);
    }

}
