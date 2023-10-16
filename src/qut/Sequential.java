package qut;

import qut.*;
import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;

public class Sequential
{
    private static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private static Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
    private static final Matrix BLOSUM_62 = BLOSUM62.Load();
    private static byte[] complement = new byte['z'];

    static
    {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }

                    
    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true)
        {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    private static boolean Homologous(PeptideSequence A, PeptideSequence B)
    {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private static NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene)
    {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
           upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }

    private static Match PredictPromoter(NucleotideSequence upStreamRegion)
    {
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    private static void ProcessDir(List<String> list, File dir)
    {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private static List<String> ListGenbankFiles(String dir)
    {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private static GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    public static void run(String referenceFile, String dir) throws FileNotFoundException, IOException
    {             
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        for (String filename : ListGenbankFiles(dir))
        {
            System.out.println(filename);
            GenbankRecord record = Parse(filename);
            for (Gene referenceGene : referenceGenes)
            {
		System.out.println(referenceGene.name);
                for (Gene gene : record.genes)
                    if (Homologous(gene.sequence, referenceGene.sequence))
                    {
                        NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
	                Match prediction = PredictPromoter(upStreamRegion);
                        if (prediction != null)
                        {
                            consensus.get(referenceGene.name).addMatch(prediction);
                            consensus.get("all").addMatch(prediction);
                        }
                    }
            }
        }

        Map<String, String> sequentialResults = new HashMap<>();
        sequentialResults.put("all", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (5430 matches)");
        sequentialResults.put("fixB", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (965 matches)");
        sequentialResults.put("carA", " Consensus: -35: T T G A C A gap: 17.7 -10: T A T A A T  (1079 matches)");
        sequentialResults.put("fixA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (896 matches)");
        sequentialResults.put("caiF", " Consensus: -35: T T C A A A gap: 18.0 -10: T A T A A T  (11 matches)");
        sequentialResults.put("caiD", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (550 matches)");
        sequentialResults.put("yaaY", " Consensus: -35: T T G T C G gap: 18.0 -10: T A T A C T  (4 matches)");
        sequentialResults.put("nhaA", " Consensus: -35: T T G A C A gap: 17.6 -10: T A T A A T  (1879 matches)");
        sequentialResults.put("folA", " Consensus: -35: T T G A C A gap: 17.5 -10: T A T A A T  (46 matches)");

        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet()) {
            String expected = sequentialResults.get(entry.getKey());
            String received = entry.getValue().toString();

            if (expected == received) {
                System.out.println(entry.getKey() + " passed!");
            } else {
                System.out.println(entry.getKey() + " failed!");
                System.out.println("Expected: '" + expected + "'\nReceived: '" + received + "'");
            }
        }
    }

    public static void main(String[] args) throws FileNotFoundException, IOException
    {
        run("referenceGenes.list", "Ecoli");
    }
}
