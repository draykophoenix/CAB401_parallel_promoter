package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class ParallelByGenbankFile
{
    private static ConcurrentMap<String, Sigma70Consensus> consensus = new ConcurrentHashMap<String, Sigma70Consensus>();
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
        List<Gene> referenceGenes = new ArrayList<>();
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

    private static Match PredictPromoter(Series sigma70_pattern, NucleotideSequence upStreamRegion)
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
        List<String> list = new ArrayList<>();
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

    public ConcurrentMap<String, Sigma70Consensus> predict(String referenceFile, String dir, Integer numThreads) throws FileNotFoundException, IOException, InterruptedException {
        List<Gene> referenceGenes = Collections.unmodifiableList(ParseReferenceGenes(referenceFile));

        List<String> genbankFiles = ListGenbankFiles(dir);

        int numFiles = genbankFiles.size();
        List<Worker> workers = new ArrayList<>();

        for (int i = 0; i < numThreads; i++) {
            int work = (numFiles + numThreads - 1) / numThreads;
            // System.out.println(i + " - work: " + work);
            int start = work * i;
            // System.out.println(i + " - start: " + start);
            int end = Math.min((start + work), numFiles);
            // System.out.println(i + " - end: " + end);

            Worker worker = new Worker(genbankFiles.subList(start, end), referenceGenes);
            workers.add(worker);
            worker.start();
        }

        for (Worker worker : workers) {
            worker.join();
        }

        return consensus;
    }

    private class Worker extends Thread {
        private List<String> filenames;
        private List<Gene> referenceGenes;

        private Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
        public Worker(List<String> filenames, List<Gene> referenceGenes) {
            this.filenames = filenames;
            this.referenceGenes = referenceGenes;
        }
        public void run() {
            String threadName = Thread.currentThread().getName() + ": ";
            long innerStartTime = System.currentTimeMillis();
            for (String filename : filenames) {
                System.out.println(threadName + filename);
                GenbankRecord record = null;
                try {
                    record = Parse(filename);
                } catch (IOException e) {
                    e.printStackTrace();
                }

                for (Gene referenceGene : referenceGenes) {
                    System.out.println(threadName + referenceGene.name);
                    for (Gene gene : record.genes)
                        if (Homologous(gene.sequence, referenceGene.sequence)) {
                            NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                            Match prediction = PredictPromoter(sigma70_pattern, upStreamRegion);
                            if (prediction != null) {
                                consensus.get(referenceGene.name).addMatch(prediction);
                                consensus.get("all").addMatch(prediction);
                            }
                        }
                }
            }
            long innerEndTime = System.currentTimeMillis();
            long innerDuration = innerEndTime - innerStartTime;
            System.out.println(threadName + "Elapsed:" + innerDuration + "ms" );

        }
    }
}
