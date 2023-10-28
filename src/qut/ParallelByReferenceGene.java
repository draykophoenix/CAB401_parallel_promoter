package qut;

import edu.au.jacobi.pattern.Match;
import edu.au.jacobi.pattern.Series;
import jaligner.BLOSUM62;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

public class ParallelByReferenceGene
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
            return new NucleotideSequence(Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
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

    public ConcurrentMap<String, Sigma70Consensus> predict(String referenceFile, String dir, int numThreads) throws FileNotFoundException, IOException, InterruptedException {
        List<Gene> referenceGenes = Collections.unmodifiableList(ParseReferenceGenes(referenceFile));

        List<String> genbankFiles = ListGenbankFiles(dir);

        List<GenbankRecord> records = new ArrayList<>(); // UNSAFE?
        for (String filename : genbankFiles) {
            try {
                records.add(Parse(filename));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        List<Worker> workers = new ArrayList<>();

        int numFiles = genbankFiles.size();
        int numReferenceGenes = referenceGenes.size();
        int numJobs = numFiles * numReferenceGenes;

        Worker.setConstants(numThreads, numReferenceGenes, numJobs);
        for (int i = 0; i < numThreads; i++) {
            Worker worker = new Worker(referenceGenes, records, i);
            workers.add(worker);
            worker.start();
        }

        for (Worker worker : workers) {
            worker.join();
        }

        return consensus;
    }

    private class Worker extends Thread {
        private static int numJobs;
        private static int numThreads;
        private static int numReferenceGenes;
        private List<Gene> referenceGenes;
        private List<GenbankRecord> records;
        private int threadsIndex;

        private Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);
        public Worker(List<Gene> referenceGenes, List<GenbankRecord> records, int threadsIndex) {
            this.referenceGenes = referenceGenes;
            this.records = records;
            this.threadsIndex = threadsIndex;
        }
        public static void setConstants(int numThreads, int numReferenceGenes, int numJobs) {
            Worker.numThreads = numThreads;
            Worker.numReferenceGenes = numReferenceGenes;
            Worker.numJobs = numJobs;

        }
        public void run() {
            String threadName = Thread.currentThread().getName() + ": ";

            int division = numJobs / numThreads;
            int remainder = numJobs % numThreads;

            int work = division;
            // Extra work for first threads
            if (threadsIndex < remainder) work+= 1;
            // Shift created from extra work
            int shift = Math.min(threadsIndex, remainder);

//            System.out.println(threadName + " - work: " + work);
            int start = division * threadsIndex + shift;
//            System.out.println(threadName + " - start: " + start);
            int end = start + work;
//            System.out.println(threadName + " - end: " + end);


            for (int i = start; i < end; i++) {
                int posReferenceGene = i % numReferenceGenes;
                int posFile = i / numReferenceGenes;

                Gene referenceGene = referenceGenes.get(posReferenceGene);
                GenbankRecord record = records.get(posFile);

                System.out.println(threadName + "index: " + i + " file: " + posFile + " gene: " + posReferenceGene + " (" + referenceGene.name + ")");

                for (Gene gene : record.genes) {
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
        }
    }
}
