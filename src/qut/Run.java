package qut;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public class Run {

    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        Map<String, Sigma70Consensus> sequentialConsensus = Sequential.run("referenceGenes.list", "Ecoli");
        //Map<String, Sigma70Consensus> parallelConsensus = Parallel.run("referenceGenes.list", "Ecoli");
        compareByString(sequentialConsensus);
    }


    private static final Map<String, String> SEQUENTIAL_RESULTS = new HashMap<>();
    static {
        SEQUENTIAL_RESULTS.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        SEQUENTIAL_RESULTS.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        SEQUENTIAL_RESULTS.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        SEQUENTIAL_RESULTS.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        SEQUENTIAL_RESULTS.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        SEQUENTIAL_RESULTS.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        SEQUENTIAL_RESULTS.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        SEQUENTIAL_RESULTS.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        SEQUENTIAL_RESULTS.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");
    }


    private static void printer(Map<String, Sigma70Consensus> consensus) {
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

    private static void compareByString(Map<String, Sigma70Consensus> consensus) {
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet()) {
            String expected = SEQUENTIAL_RESULTS.get(entry.getKey());
            String received = entry.getValue().toString();

            if (Objects.equals(expected, received)) {
                System.out.println(entry.getKey() + " passed!");
            } else {
                System.out.println(entry.getKey() + " failed!");
                System.out.println("Expected: <" + expected + ">");
                System.out.println("Received: <" + received + ">");
            }
        }
    }
    private static void compareDirect(Map<String, Sigma70Consensus> sequentialConsensus, Map<String, Sigma70Consensus> parallelConsensus) {
        if (sequentialConsensus.equals(parallelConsensus)) {
            System.out.println("Direct compare passed!");
            System.out.println("Sequential:");
            printer(sequentialConsensus);
            System.out.println("Parallel:");
            printer(parallelConsensus);
        }

    }

}
