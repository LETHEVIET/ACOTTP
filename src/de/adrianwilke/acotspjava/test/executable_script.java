package de.adrianwilke.acotspjava.test;

import de.adrianwilke.acotspjava.AcoTsp;

import java.util.LinkedList;
import java.util.List;

public class executable_script {

    private static String input_path = "";
    private static String output_path = "";

    private static List<String> argList = new LinkedList<String>();

    public static void main(String[] args) {
        
        input_path = args[0];
        output_path = args[1];
        System.out.println(input_path);
        System.out.println(output_path);

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
        argList.add("10");

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
        argList.add(input_path);

        argList.add("-do");
        argList.add(output_path);
    }
}
