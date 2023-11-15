package ACOTTP.test;

import ACOTTP.ACOTTP;

import java.util.LinkedList;
import java.util.List;

public class ACOTTP_CLI {

    private static final String TEST_FILE = "instances/a280_n279_bounded-strongly-corr_01.ttp";

    private static final List<String> argList = new LinkedList<String>();

    static long timelimit = 180;

    public static void main(String[] args) {
        System.out.println(TEST_FILE);

        putTestParameters();

        ACOTTP.main(argList.toArray(new String[0]));

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
        argList.add(timelimit + "");

        argList.add("--tours");
        argList.add("100");

        argList.add("--tries");
        argList.add("1");

        argList.add("--mmas");

        argList.add("--localsearch");
        argList.add("3");

        argList.add("--seed");
        argList.add("-1618462292");// 0:no local search  1:2-opt  2:2.5-opt  3:3-opt

        argList.add("-i");
        argList.add(TEST_FILE);

    }

}
