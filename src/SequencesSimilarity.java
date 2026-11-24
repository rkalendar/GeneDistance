import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.stream.IntStream;

public final class SequencesSimilarity {

    /**
     * The k-mer sets: 4/5/6/8/10 are the symmetric (palindromic) ones,
     * 41/61/81/101 the non-symmetric ones.
     */
    public static final int[] MODELS = {4, 5, 6, 8, 10, 41, 61, 81, 101};

    public static boolean isModel(int m) {
        for (int v : MODELS) {
            if (v == m) {
                return true;
            }
        }
        return false;
    }

    /**
     * The set itself, so that no count can drift away from the table it
     * describes.
     */
    private static String[] mers(int m) {
        switch (m) {
            case 5:
                return RE5MERS;
            case 6:
                return RE6MERS;
            case 8:
                return RE8MERS;
            case 10:
                return RE10MERS;
            case 41:
                return K4MERS;
            case 61:
                return R6MERS;
            case 81:
                return R8MERS;
            case 101:
                return R10MERS;
            default:
                return RE4MERS;   // 4
        }
    }

    public static String modelInfo(int m) {
        String s;
        switch (m) {
            case 5:
                s = "symmetric 5-mers";
                break;
            case 6:
                s = "symmetric 6-mers";
                break;
            case 8:
                s = "symmetric 8-mers";
                break;
            case 10:
                s = "symmetric 10-mers";
                break;
            case 41:
                s = "non-symmetric 4-mers";
                break;
            case 61:
                s = "non-symmetric 6-mers";
                break;
            case 81:
                s = "non-symmetric 8-mers";
                break;
            case 101:
                s = "non-symmetric 10-mers";
                break;
            default:
                s = "symmetric and mirror 4-mers";
                break;
        }
        return s + ", " + mers(m).length;
    }

    public SequencesSimilarity(int m) {
        if (!isModel(m)) {
            m = 4;
        }
        mod = m;
        String[] mers = mers(m);
        kmer = mers[0].length();

        pt = new HashMap<>(2 * mers.length);
        for (int i = 0; i < mers.length; i++) {
            pt.put(mers[i], i);
        }
        keys = new String[pt.size()];
        for (Map.Entry<String, Integer> entry : pt.entrySet()) {
            keys[entry.getValue()] = entry.getKey();
        }
        codeIndex = buildCodeIndex();
        buildSpacingFeatures();
    }

    /**
     * The features of the spacing measure: every word of the set pooled with
     * its reverse complement, so that the count of a feature is the same
     * whichever strand was sequenced. A palindromic word is its own reverse
     * complement, so the palindromic sets are left exactly as they are; in the
     * other sets a word and its reverse complement, counted separately before,
     * become one feature (the reverse complement is counted even when the set
     * does not list it).
     */
    private void buildSpacingFeatures() {
        if (codeIndex == null) {       // not a plain acgt set: nothing to pool
            spIndex = null;
            spKeys = keys;
            return;
        }
        final int[] tbl = new int[1 << (2 * kmer)];
        Arrays.fill(tbl, -1);
        final List<String> names = new ArrayList<>();
        final HashMap<String, Integer> seen = new HashMap<>();
        for (String w : keys) {
            final String rc = revcomp(w);
            final String rep = (w.compareTo(rc) <= 0) ? w : rc;
            if (seen.containsKey(rep)) {
                continue;
            }
            final String other = rep.equals(w) ? rc : w;   // the reverse complement of rep
            final int f = names.size();
            seen.put(rep, f);
            names.add(rep.equals(other) ? rep : rep + "/" + other);
            tbl[codeOf(rep)] = f;
            tbl[codeOf(other)] = f;
        }
        spIndex = tbl;
        spKeys = names.toArray(new String[0]);
    }

    private static int codeOf(String w) {
        int code = 0;
        for (int i = 0; i < w.length(); i++) {
            code = (code << 2) | BASE[w.charAt(i)];
        }
        return code;
    }

    /**
     * The number of features the spacing measure compares (words pooled with
     * their reverse complement).
     */
    public int getSpacingCount() {
        return spKeys.length;
    }

    public int getSequenceCount() {
        return seq == null ? 0 : seq.length;
    }

    public int getKmerCount() {
        return pt.size();
    }

    public int getModel() {
        return mod;
    }

    public String getReportFile() {
        return reportfile;
    }

    public String getMegaFile() {
        return megafile;
    }

    /**
     * Tags the report files of a mode of its own, so that it cannot overwrite
     * the report of another one. Call it before SetFolder.
     */
    public void SetTag(String t) {
        tag = t;
    }

    /**
     * Names the report files outright, for a mode that does not use a k-mer
     * set. It may be called after SetFolder, since a mode may only know its own
     * name once it has seen the sequences.
     */
    public void SetSuffix(String s) {
        suffix = s;
        if (basePath != null) {
            reportfile = basePath + suffix + ".xls";
            megafile = basePath + suffix + ".meg";
        }
    }

    /**
     * Names the report files for a streaming mode that never calls SetFolder.
     */
    public void SetPaths(String filePath, String sfx) {
        basePath = filePath;
        reportfile = filePath + sfx + ".xls";
        megafile = filePath + sfx + ".meg";
    }

    public int getLongestSequence() {
        int n = 0;
        for (String s : seq) {
            n = Math.max(n, s.length());
        }
        return n;
    }

    public void SetFolder(String baseName, String[] tfiles, String filePath) throws IOException {
        fPath = baseName;
        basePath = filePath;
        String sfx = (suffix != null) ? suffix : ("_k" + mod + tag);
        reportfile = filePath + sfx + ".xls";  //+ "result.xls"
        megafile = filePath + sfx + ".meg";    //+ "result.mega"

        final List<String> seqs = new ArrayList<>();
        final List<String> names = new ArrayList<>();

        for (String infile : tfiles) {
            if (infile != null) {
                System.out.println("Target file: " + infile);
            }
        }
        // Streamed record by record, so a FASTA larger than 2 GB - which
        // Files.readAllBytes cannot hold in one byte[] - is read all the same.
        // The records are identical, byte for byte, to the old whole-file reader.
        try {
            StreamingFasta.forEachRecord(tfiles, (name, s) -> {
                names.add(name);
                seqs.add(s);
            });
        } catch (IOException e) {
            System.err.println("Failed to read the input: " + e.getMessage());
        }
        if (seqs.isEmpty()) {
            System.out.println("There is no sequence(s).\nFile format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
        } else if (seqs.size() > 1) {
            System.out.println("Target FASTA sequences = " + seqs.size());
        }
        seq = seqs.toArray(new String[0]);
        sname = names.toArray(new String[0]);
    }

    /**
     * The whole analysis: the MEGA distance matrix (.meg) and the report (.xls)
     * with the k-mer counts, the k-mer distances and the similarity matrix.
     *
     * The two files used to be produced by two passes (RunTree and RunPattern)
     * that each rebuilt the same k-mer count matrix and each ran the same
     * O(nseq^2 * nkmers^2) comparison. The two metrics only differ in how the
     * terms of that comparison are summed, so they are now gathered together in
     * a single pass.
     */
    public void Run() throws IOException {
        final int nseq = seq.length;
        final int[] slen = new int[nseq];
        final int[][] m2 = new int[nseq][];
        IntStream.range(0, nseq).parallel().forEach(j -> {
            slen[j] = seq[j].length();
            m2[j] = countSpacingOne(seq[j]);
        });
        runFrom(sname, slen, m2);
    }

    /**
     * As {@link #Run()}, but the counts and lengths are built by streaming the
     * files one record at a time, so the sequences are never all in memory -
     * only the small per-record count vectors. This is sequential like -d2star,
     * so it needs no random access and streams directly (no SeqStore). Byte for
     * byte the same output on the same input.
     */
    public void RunStream(String[] files) throws IOException {
        final List<String> names = new ArrayList<>();
        final List<Integer> lens = new ArrayList<>();
        final List<int[]> counts = new ArrayList<>();
        StreamingFasta.forEachRecord(files, (name, s) -> {
            names.add(name);
            lens.add(s.length());
            counts.add(countSpacingOne(s));
        });
        final int n = names.size();
        final int[] slen = new int[n];
        for (int i = 0; i < n; i++) {
            slen[i] = lens.get(i);
        }
        runFrom(names.toArray(new String[0]), slen, counts.toArray(new int[0][]));
    }

    private void runFrom(String[] names, int[] slen, int[][] m2) throws IOException {
        final int nseq = names.length;
        final int nkmers = spKeys.length;
        final String[] keys = spKeys;   // the pooled features, not the raw words of the set

        // the names as they appear in the report (the .meg keeps the original ones)
        final String[] xname = new String[nseq];
        for (int j = 0; j < nseq; j++) {
            xname[j] = names[j].replace(' ', '_').replace(',', '_');
        }

        // lower triangles: similarity (%) for the report, distance for MEGA
        final double[][] similarity = new double[nseq][];
        final double[][] distance = new double[nseq][];
        for (int i = 1; i < nseq; i++) {
            similarity[i] = new double[i];
            distance[i] = new double[i];
        }

        try (BufferedWriter w = newWriter(reportfile)) {
            StringBuilder sr = new StringBuilder(1 << 17);

            // printing kmers number
            sr.append("ID/kmer");
            for (int n1 = 0; n1 < nkmers; n1++) {
                sr.append("\t").append(keys[n1]);
            }
            sr.append("\n");
            for (int j = 0; j < nseq; j++) {
                sr.append(xname[j]);
                final int[] row = m2[j];
                for (int n2 = 0; n2 < nkmers; n2++) {
                    sr.append("\t").append(row[n2]);
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n");

            // printing the average spacing L/q of every feature seen at least twice (0 otherwise).
            // It is shown for the reader; the comparison below works on the counts themselves,
            // since the sequence length cancels in every ratio of two spacings.
            sr.append("ID/kmer");
            for (int n1 = 0; n1 < nkmers; n1++) {
                sr.append("\t").append(keys[n1]);
            }
            sr.append("\n");
            for (int j = 0; j < nseq; j++) {
                sr.append(xname[j]);
                final int[] row = m2[j];
                final int len = slen[j];
                for (int i = 0; i < nkmers; i++) {
                    sr.append("\t").append((row[i] > minkmer) ? Math.round((double) len / row[i]) : 0);
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n");

            // comparative analysis of each vector
            pairwise(m2, nseq, nkmers, similarity, distance);

            sr.append("\n\n");
            for (int j = 0; j < nseq; j++) {
                sr.append("\t").append(xname[j]);
            }
            sr.append("\n");
            for (int i = 0; i < nseq; i++) {
                sr.append(xname[i]);
                for (int n = 0; n < nseq; n++) {
                    if (i == n) {
                        sr.append("\t").append("100");
                    } else {
                        sr.append("\t").append(Math.round(i > n ? similarity[i][n] : similarity[n][i]));
                    }
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n");
            w.write(sr.toString());
        }

        writeMega(names, distance, nseq);
    }

    private void writeMega(String[] names, double[][] distance, int nseq) throws IOException {
        try (BufferedWriter w = newWriter(megafile)) {
            StringBuilder mg = new StringBuilder(1 << 17);
            mg.append("#mega\n");
            mg.append("!Format DataType=Distance DataFormat=LowerLeft NTaxa=").append(nseq).append(";\n\n");
            for (int j = 0; j < nseq; j++) {
                String s = names[j];
                if (s.length() > 48) {
                    s = s.substring(0, 48);
                }
                mg.append("#").append(s).append("\n");
            }
            mg.append("\n\n");

            for (int i = 1; i < nseq; i++) {
                final double[] row = distance[i];
                for (int n = 0; n < i; n++) {
                    mg.append(" ").append(String.format(Locale.US, "%.5f", row[n]));
                }
                mg.append("\n");
                flush(w, mg);
            }
            w.write(mg.toString());
        }
    }

    /**
     * Compares the spacing signatures of every pair of sequences.
     *
     * The spacing of feature i is x_i = L / q_i, and two sequences A and B are
     * compared through the ratios of the spacings of every pair of features:
     *
     * F_ij = min(r_ij(A), r_ij(B)) / max(r_ij(A), r_ij(B)), r_ij = x_i / x_j =
     * q_j / q_i S = mean of F_ij over the pairs i < j of features seen at least
     * twice in both
     *
     * The length L cancels in r_ij, so the counts are all that is needed. With
     * d_i = ln q_i(A) - ln q_i(B), the log fold-change of feature i,
     *
     * F_ij = exp(-|d_i - d_j|)
     *
     * which depends only on the differences of d. Sorted in ascending order,
     * |d_i - d_j| = d_j - d_i for i < j, so
     *
     * sum_{i<j} F_ij = sum_j exp(-d_j) * (exp(d_1) + ... + exp(d_{j-1}))
     *
     * and the inner sums are prefix sums: S is obtained exactly in O(n log n)
     * for n features, where the pairwise definition needs O(n^2). Every pair of
     * features carries the same weight, whatever the order of the table.
     *
     * S = 1 exactly when the two count vectors are proportional (same
     * signature, any length); the report holds 100 * S and the MEGA file the
     * distance 1 - S. S measures composition, not homology: two unrelated
     * sequences of the same composition tend to S = 1 as they grow longer.
     */
    private void pairwise(int[][] m2, int nseq, int nkmers, double[][] similarity, double[][] distance) {
        // the natural log of every count once; NaN marks a feature seen fewer than twice
        final double[][] lg = new double[nseq][nkmers];
        IntStream.range(0, nseq).parallel().forEach(j -> {
            final int[] row = m2[j];
            final double[] l = lg[j];
            for (int n = 0; n < nkmers; n++) {
                l[n] = (row[n] > minkmer) ? Math.log(row[n]) : Double.NaN;
            }
        });

        IntStream.range(1, nseq).parallel().forEach(i -> {
            final double[] li = lg[i];
            final double[] d = new double[nkmers];
            for (int j = 0; j < i; j++) {
                final double[] lj = lg[j];
                int n = 0;
                for (int f = 0; f < nkmers; f++) {
                    final double x = li[f] - lj[f];
                    if (x == x) {              // defined in both sequences
                        d[n++] = x;
                    }
                }
                if (n < MIN_SHARED) {          // under-determined pair: the largest in-scale distance
                    similarity[i][j] = 0;
                    distance[i][j] = 1;
                    continue;
                }
                Arrays.sort(d, 0, n);
                final double base = d[0];      // counts are ints, so d - base stays far below overflow
                double prefix = 0, total = 0;
                for (int t = 0; t < n; t++) {
                    final double e = Math.exp(d[t] - base);
                    total += prefix / e;
                    prefix += e;
                }
                final double s = Math.min(1.0, 2.0 * total / ((double) n * (n - 1)));
                similarity[i][j] = 100 * s;
                distance[i][j] = 1 - s;
            }
        });
    }

    /**
     * A pair is compared only on this many shared features or more.
     */
    private static final int MIN_SHARED = 4;

    /**
     * The k a containment run should use, when the user names none: large
     * enough that a k-mer of the longest sequence is unlikely (p = 1%) to occur
     * in an unrelated sequence by chance, k = log4(L * 99), and never so large
     * that a diverged copy would lose all of its k-mers (a k-mer survives a
     * divergence d with probability (1-d)^k). Mash, Ondov et al. (2016), eq. 2.
     */
    public static int containmentK(long longest) {
        int k = (int) Math.ceil(Math.log(Math.max(1L, longest) * 99.0) / Math.log(4));
        return Math.max(11, Math.min(24, k | 1));   // odd k: a k-mer cannot be its own reverse complement
    }

    /**
     * How much of each sequence occurs in each other one: the CONTAINMENT
     * index, on the whole k-mer space rather than on a curated set.
     *
     * C(A,B) = |kmers(A) & kmers(B)| / |kmers(A)|
     *
     * It is asymmetric on purpose - that is the whole point. A measure of the
     * composition of a sequence as a whole, however good, cannot see a 300 bp
     * element inside a 100 kb sequence, because the element is 0.3% of what it
     * measures; containment asks a different question, "how much of A is in
     * B?", and answers 1.0 there. Reading the report, the row is A and the
     * column is B.
     *
     * The k-mers are canonical (a k-mer and its reverse complement count as
     * one), so an element inserted the other way round is still found.
     *
     * From the containment follows an estimate of the identity of the shared
     * part, p = C^(1/k), the estimator of Mash Screen (Ondov et al., Genome
     * Biology, 2019); the MEGA file holds 1 - p, taking the larger of the two
     * containments of a pair.
     */
    public void RunContainment(int k) throws IOException {
        RunContainment(k, 1);
    }

    /**
     * Containment with a FracMinHash sketch: only k-mers whose hash falls in
     * the bottom 1/{@code scaled} of the hash space are kept, so a 2 Gb
     * genome's set shrinks by that factor - from ~16 GB to ~16 MB at scaled =
     * 1000 - while the containment estimated from the kept fraction stays
     * unbiased. {@code scaled
     * = 1} keeps every k-mer, the exact measure, byte for byte as before.
     */
    public void RunContainment(int k, int scaled) throws IOException {
        final int nseq = seq.length;
        final long maxHash = maxHashForScale(scaled);
        final long[][] ks = new long[nseq][];
        final long[] lengths = new long[nseq];
        IntStream.range(0, nseq).parallel().forEach(j -> {
            ks[j] = canonical(seq[j], k, maxHash);
            lengths[j] = seq[j].length();
        });
        writeContainment(sname, lengths, ks, k, scaled);
    }

    /**
     * Containment with the k-mer sets built by streaming the files, so the
     * sequences are never held in memory - a 2 Gb genome costs only its sketch,
     * whatever the number of chromosomes. The result is byte for byte the one
     * the in-memory {@link #RunContainment(int,int)} produces on the same
     * input.
     */
    public void RunContainmentStream(String[] files, int k, int scaled) throws IOException {
        final long maxHash = maxHashForScale(scaled);
        final List<String> names = new ArrayList<>();
        final List<Long> lens = new ArrayList<>();
        final List<long[]> sk = new ArrayList<>();
        StreamingFasta.forEachSketch(files, k, maxHash, (name, length, sketch) -> {
            names.add(name);
            lens.add(length);
            sk.add(sketch);
        });
        final int nseq = names.size();
        final long[] lengths = new long[nseq];
        for (int i = 0; i < nseq; i++) {
            lengths[i] = lens.get(i);
        }
        writeContainment(names.toArray(new String[0]), lengths, sk.toArray(new long[0][]), k, scaled);
    }

    private void writeContainment(String[] names, long[] lengths, long[][] ks, int k, int scaled) throws IOException {
        final int nseq = names.length;
        final int[][] shared = new int[nseq][nseq];
        IntStream.range(0, nseq).parallel().forEach(i -> {
            for (int j = i + 1; j < nseq; j++) {
                final int s = intersection(ks[i], ks[j]);
                shared[i][j] = s;
                shared[j][i] = s;
            }
        });

        final String[] xname = new String[nseq];
        for (int j = 0; j < nseq; j++) {
            xname[j] = names[j].replace(' ', '_').replace(',', '_');
        }

        try (BufferedWriter w = newWriter(reportfile)) {
            StringBuilder sr = new StringBuilder(1 << 17);
            sr.append("containment, k=").append(k);
            if (scaled > 1) {
                sr.append(", FracMinHash scaled=").append(scaled).append(" (estimated)");
            }
            sr.append("\n\n");
            sr.append("ID\tlength\t").append(scaled > 1 ? "sketch kmers" : "distinct kmers").append("\n");
            for (int j = 0; j < nseq; j++) {
                sr.append(xname[j]).append("\t").append(lengths[j]).append("\t").append(ks[j].length).append("\n");
                flush(w, sr);
            }
            sr.append("\n\n");

            sr.append("containment %: the share of the k-mers of the ROW found in the COLUMN");
            for (int j = 0; j < nseq; j++) {
                sr.append("\t").append(xname[j]);
            }
            sr.append("\n");
            for (int i = 0; i < nseq; i++) {
                sr.append(xname[i]);
                for (int j = 0; j < nseq; j++) {
                    sr.append("\t").append(i == j ? 100 : Math.round(100 * containment(shared[i][j], ks[i].length)));
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n\n");

            sr.append("identity % of the shared part, from the larger containment of the pair");
            for (int j = 0; j < nseq; j++) {
                sr.append("\t").append(xname[j]);
            }
            sr.append("\n");
            for (int i = 0; i < nseq; i++) {
                sr.append(xname[i]);
                for (int j = 0; j < nseq; j++) {
                    sr.append("\t").append(i == j ? 100 : Math.round(100 * identity(shared[i][j], ks[i].length, ks[j].length, k, scaled)));
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n");
            w.write(sr.toString());
        }

        try (BufferedWriter w = newWriter(megafile)) {
            StringBuilder mg = new StringBuilder(1 << 17);
            mg.append("#mega\n");
            mg.append("!Format DataType=Distance DataFormat=LowerLeft NTaxa=").append(nseq).append(";\n\n");
            for (int j = 0; j < nseq; j++) {
                String s = names[j];
                if (s.length() > 48) {
                    s = s.substring(0, 48);
                }
                mg.append("#").append(s).append("\n");
            }
            mg.append("\n\n");
            for (int i = 1; i < nseq; i++) {
                for (int j = 0; j < i; j++) {
                    final double p = identity(shared[i][j], ks[i].length, ks[j].length, k, scaled);
                    mg.append(" ").append(String.format(Locale.US, "%.5f", 1 - p));
                }
                mg.append("\n");
                flush(w, mg);
            }
            w.write(mg.toString());
        }
    }

    private static double containment(int shared, int size) {
        return size == 0 ? 0 : (double) shared / size;
    }

    /**
     * the larger of the two containments of a pair, so that it does not depend
     * on which is longer
     */
    private static double maxContainment(int shared, int a, int b) {
        int m = Math.min(a, b);
        return m == 0 ? 0 : (double) shared / m;
    }

    /**
     * The identity of the shared part, from the containment: a k-mer survives a
     * divergence of d with probability (1-d)^k, so C ~ p^k and p ~ C^(1/k).
     *
     * A handful of k-mers are shared by chance alone - a given k-mer of A turns
     * up in a sequence of n distinct k-mers with probability about n/4^k - and
     * C^(1/k) turns even that into a fine-looking 60%: unless the containment
     * stands well clear of what chance gives, nothing is reported.
     */
    private static double identity(int shared, int a, int b, int k, int scaled) {
        final double c = maxContainment(shared, a, b);
        if (c <= 0 || shared < 5) {
            return 0;
        }
        // With a FracMinHash sketch the sets a, b hold ~1/scaled of the real
        // k-mers, so the real set that competes for a chance match is scaled
        // times larger. At scaled == 1 this is the exact set and the factor is 1.
        final double chance = Math.max(a, b) * (double) scaled / Math.pow(4, k);
        if (c <= 10 * chance) {
            return 0;
        }
        return Math.pow(c, 1.0 / k);
    }

    /**
     * Walks a window the size of the query along every other sequence and
     * reports where the query is found: the containment tells how much of a
     * sequence occurs in another, this tells WHERE.
     *
     * The window holds the k-mer positions of the stretch under it; a k-mer
     * entering it that belongs to the query raises the count, one leaving it
     * lowers it, so a window of any size costs the same and the sequence is
     * walked once, base by base. What the window carries is the containment of
     * the query in that very stretch.
     *
     * A few k-mers of the query are met by chance in any sequence, and a scan
     * asks the question at every base, so the best of a hundred thousand
     * windows is bound to look good. Every stretch reported is therefore
     * weighed against what chance alone gives over the whole scan: E is the
     * number of stretches this good that a sequence of no relation would be
     * expected to yield, in the manner of the E-value of BLAST, and only what
     * comes below E = 0.01 is called a hit.
     */
    public void RunScan(int k, int query) throws IOException {
        RunScan(k, query, 1);
    }

    /**
     * Scan with a FracMinHash sketch: only k-mers whose hash is in the bottom
     * 1/{@code scaled} of the space are indexed and searched. The query's index
     * shrinks by that factor and the walk does a hash lookup only for the kept
     * fraction of positions, so a large query no longer needs its whole k-mer
     * set in memory. The containment and E-value are computed in sketch space,
     * so they estimate the same quantities the exact scan reports. {@code scaled
     * = 1} indexes every k-mer - the exact scan, byte for byte as before.
     */
    public void RunScan(int k, int query, int scaled) throws IOException {
        final int nseq = seq.length;
        final long maxHash = maxHashForScale(scaled);
        final boolean sketch = (maxHash != -1L);

        final long[] qk = canonical(seq[query], k, maxHash);
        final HashMap<Long, Integer> qi = queryIndex(qk);
        final int nq = qk.length;
        final int window = seq[query].length();

        System.out.println("Query: " + sname[query] + " (" + window + " bp, " + nq + " distinct k-mers of length " + k + ")");
        if (nq == 0) {
            emptyScan(k, name(query), window);
            return;
        }
        final int[] qComp = compositionClasses(qk, k);

        final List<String> lines = new ArrayList<>();
        final List<String> hits = new ArrayList<>();
        final List<int[]> profiles = new ArrayList<>();
        final List<String> profileOf = new ArrayList<>();

        for (int j = 0; j < nseq; j++) {
            if (j != query) {
                scanSubject(sname[j], seq[j], k, scaled, nq, qi, window, qComp, maxHash, sketch,
                        lines, hits, profiles, profileOf);
            }
        }
        reportScan(name(query), window, nq, k, scaled, lines, hits, profiles, profileOf);
    }

    /**
     * Scan by streaming the files, so neither the subjects nor the query are
     * ever all in memory at once: the query is read once for its k-mer index,
     * then each subject is streamed in, scanned and dropped. Byte for byte the
     * result of the in-memory {@link #RunScan(int,int,int)} on the same input.
     */
    public void RunScanStream(String[] files, int query, int k, int scaled) throws IOException {
        final long maxHash = maxHashForScale(scaled);
        final boolean sketch = (maxHash != -1L);

        // pass 1: the query record's own k-mer index (one record held, then dropped)
        final String[] qName = {""};
        final long[][] qkHold = {new long[0]};
        final int[] cur = {0};
        StreamingFasta.forEachRecord(files, (name, s) -> {
            if (cur[0] == query) {
                qName[0] = name;
                qkHold[0] = canonical(s, k, maxHash);
            }
            cur[0]++;
        });
        final long[] qk = qkHold[0];
        final HashMap<Long, Integer> qi = queryIndex(qk);
        final int nq = qk.length;
        final int window = (int) recordLength(files, query);
        final String qx = qName[0].replace(' ', '_').replace(',', '_');

        System.out.println("Query: " + qName[0] + " (" + window + " bp, " + nq + " distinct k-mers of length " + k + ")");
        if (nq == 0) {
            emptyScan(k, qx, window);
            return;
        }
        final int[] qComp = compositionClasses(qk, k);

        final List<String> lines = new ArrayList<>();
        final List<String> hits = new ArrayList<>();
        final List<int[]> profiles = new ArrayList<>();
        final List<String> profileOf = new ArrayList<>();

        // pass 2: each subject streamed in, scanned, dropped
        final int[] jx = {0};
        StreamingFasta.forEachRecord(files, (name, s) -> {
            if (jx[0] != query) {
                scanSubject(name, s, k, scaled, nq, qi, window, qComp, maxHash, sketch,
                        lines, hits, profiles, profileOf);
            }
            jx[0]++;
        });
        reportScan(qx, window, nq, k, scaled, lines, hits, profiles, profileOf);
    }

    private static HashMap<Long, Integer> queryIndex(long[] qk) {
        final HashMap<Long, Integer> qi = new HashMap<>(2 * qk.length);
        for (int i = 0; i < qk.length; i++) {
            qi.put(qk[i], i);
        }
        return qi;
    }

    /**
     * The k-mers of the query, counted by base composition: entry a*(k+1)^2 +
     * c*(k+1) + g is the number of them made of a adenines, c cytosines and g
     * guanines (the rest thymines). How often chance puts a k-mer in a subject
     * depends on the subject's base composition and on that of the k-mer only,
     * so this histogram, built once, gives the chance level against any subject
     * in a few thousand operations.
     */
    private static int[] compositionClasses(long[] qk, int k) {
        final int m = k + 1;
        final int[] h = new int[m * m * m];
        for (long code : qk) {
            int a = 0, c = 0, g = 0;
            long x = code;
            for (int i = 0; i < k; i++) {
                final int b = (int) (x & 3);
                if (b == 0) {
                    a++;
                } else if (b == 1) {
                    c++;
                } else if (b == 2) {
                    g++;
                }
                x >>>= 2;
            }
            h[(a * m + c) * m + g]++;
        }
        return h;
    }

    /**
     * The number of distinct k-mers of the query that chance alone puts under a
     * window of {@code wp} k-mer positions of a subject whose base frequencies
     * are {@code f}: the sum, over the k-mers of the query, of 1 - (1 - p)^wp,
     * p being the probability that a position of the subject starts that k-mer
     * or its reverse complement. Taking p = 4^-k for every k-mer, as was done
     * before, is wrong by orders of magnitude for an AT-rich k-mer in AT-rich
     * DNA, and made unrelated sequences look significant.
     */
    private static double chanceLevel(int[] qComp, int k, double[] f, int wp) {
        final int m = k + 1;
        final double la = Math.log(f[0]), lc = Math.log(f[1]), lgg = Math.log(f[2]), lt = Math.log(f[3]);
        double lambda = 0;
        for (int a = 0; a <= k; a++) {
            for (int c = 0; a + c <= k; c++) {
                for (int g = 0; a + c + g <= k; g++) {
                    final int n = qComp[(a * m + c) * m + g];
                    if (n == 0) {
                        continue;
                    }
                    final int t = k - a - c - g;
                    // the k-mer itself, and its reverse complement (a<->t, c<->g)
                    final double pw = Math.exp(a * la + c * lc + g * lgg + t * lt)
                            + Math.exp(t * la + g * lc + c * lgg + a * lt);
                    lambda += n * -Math.expm1(wp * Math.log1p(-Math.min(pw, 0.5)));
                }
            }
        }
        return lambda;
    }

    /**
     * Base frequencies of a sequence (a, c, g, t), with a pseudo-count so that
     * none is zero.
     */
    private static double[] baseFrequencies(String s) {
        final long[] n = new long[4];
        for (int i = 0, e = s.length(); i < e; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b >= 0) {
                n[b]++;
            }
        }
        final double tot = n[0] + n[1] + n[2] + n[3] + 4.0;
        return new double[]{(n[0] + 1) / tot, (n[1] + 1) / tot, (n[2] + 1) / tot, (n[3] + 1) / tot};
    }

    private static long recordLength(String[] files, int index) throws IOException {
        final long[] lens = StreamingFasta.recordLengths(files);
        return (index >= 0 && index < lens.length) ? lens[index] : 0;
    }

    private void emptyScan(int k, String queryName, int window) throws IOException {
        System.out.println("The query holds no k-mer of that length: nothing to look for.");
        try (BufferedWriter w = newWriter(reportfile)) {
            w.write("scan, k=" + k + ", query=" + queryName + ", window=" + window
                    + " bp\n\nThe query holds no k-mer of length " + k + ": nothing to look for.\n");
        }
    }

    /**
     * Walks the window along one subject and appends its best window, its hit
     * regions and its profile. Identical whether the subject came from memory
     * or a stream - the only inputs are the subject's name and bases.
     */
    private void scanSubject(String rawName, String s, int k, int scaled, int nq,
            HashMap<Long, Integer> qi, int window, int[] qComp, long maxHash, boolean sketch,
            List<String> lines, List<String> hits, List<int[]> profiles, List<String> profileOf) {
        final int len = s.length();
        if (len < k) {
            return;
        }
        final String xn = rawName.replace(' ', '_').replace(',', '_');
        final int win = Math.min(window, len);      // a subject shorter than the query is taken whole
        final int wp = win - k + 1;                 // k-mer positions under the window
        final int nw = len - win + 1;               // windows walked

        // what chance alone puts under a window, and the count a window must reach.
        // The chance level is that of THIS subject: its base composition decides how
        // often each k-mer of the query turns up in it. (With a sketch the query set
        // is the sketched one, and every position of the window can still hold any of
        // its k-mers, so the same expression stands.)
        final double[] freq = baseFrequencies(s);
        final double lambda = Math.max(1e-12, chanceLevel(qComp, k, freq, wp));
        final int floor = Math.max(5, (int) Math.ceil(10 * lambda));
        // The k-mers that chance puts under a window are not independent: a chance match
        // longer than k holds several consecutive k-mers of the query. Each k-mer met is
        // followed by another with the probability that one more base agrees, sum of f^2,
        // which the E-value takes into account (see evalue). The number of separate
        // matches - a k-mer of the query met where the position before held none - is
        // reported beside it: one long exact match and many short ones are told apart.
        final double extend = freq[0] * freq[0] + freq[1] * freq[1] + freq[2] * freq[2] + freq[3] * freq[3];

        final int gap = Math.max(100, win / 10);    // k-mers farther apart than this are not one thing
        final int[] cnt = new int[nq];
        final int[] ring = new int[wp];
        final boolean[] ringRun = new boolean[wp];  // a match starts on this position
        final int stride = Math.max(1, nw / 2000);
        final int[] profile = new int[(nw + stride - 1) / stride];

        int distinct = 0, best = 0, bestAt = 0;
        int runs = 0, bestRuns = 0, peakRuns = 0;   // separate matches under the window; the most seen, overall and in the stretch
        boolean prevHit = false;
        int start = -1, peak = 0, peakAt = 0;       // the stretch being followed
        final List<int[]> found = new ArrayList<>();
        int[] hit = new int[64];                    // where a k-mer of the query is met, in order
        int nhit = 0;

        long fwd = 0, rev = 0;
        int valid = 0;
        final long mask = (k == 32) ? -1L : ((1L << (2 * k)) - 1);
        final int top = 2 * (k - 1);

        for (int i = 0; i < len; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b < 0) {
                valid = 0;
                fwd = 0;
                rev = 0;
            } else {
                fwd = ((fwd << 2) | b) & mask;
                rev = (rev >>> 2) | ((long) (3 - b) << top);
                if (valid < k) {
                    valid++;
                }
            }
            if (i < k - 1) {
                continue;
            }
            final int p = i - k + 1;                // the k-mer starting here
            int idx = -1;
            if (valid == k) {
                final long can = Math.min(fwd, rev);
                // a sketched k-mer of the subject can only match the query sketch
                // if it passes the same hash filter; skip the lookup otherwise
                if (!sketch || Long.compareUnsigned(mix64(can), maxHash) <= 0) {
                    final Integer of = qi.get(can);
                    if (of != null) {
                        idx = of;
                    }
                }
            }

            final int slot = p % wp;
            if (p >= wp) {                          // this one leaves the window
                final int old = ring[slot];
                if (old >= 0 && --cnt[old] == 0) {
                    distinct--;
                }
                if (ringRun[slot]) {
                    runs--;
                }
            }
            ring[slot] = idx;
            final boolean runStart = idx >= 0 && !prevHit;
            ringRun[slot] = runStart;
            if (runStart) {
                runs++;
            }
            prevHit = idx >= 0;
            if (idx >= 0) {
                if (cnt[idx]++ == 0) {
                    distinct++;
                }
                if (nhit == hit.length) {
                    hit = Arrays.copyOf(hit, 2 * nhit);
                }
                hit[nhit++] = p;
            }
            if (p < wp - 1) {
                continue;
            }

            final int w = p - wp + 1;               // the window starts on this base
            if (w % stride == 0) {
                profile[w / stride] = distinct;
            }
            if (distinct > best) {
                best = distinct;
                bestAt = w;
            }
            if (runs > bestRuns) {
                bestRuns = runs;
            }
            if (distinct >= floor) {                // follow the stretch while it lasts
                if (start < 0) {
                    start = w;
                    peak = 0;
                    peakRuns = 0;
                }
                if (distinct > peak) {
                    peak = distinct;
                    peakAt = w;
                }
                if (runs > peakRuns) {
                    peakRuns = runs;
                }
            } else if (start >= 0) {
                found.add(withRuns(region(hit, nhit, start, w - 2 + wp, k, len, gap, peak, peakAt), peakRuns));
                start = -1;
            }
        }
        if (start >= 0) {
            found.add(withRuns(region(hit, nhit, start, len - k, k, len, gap, peak, peakAt), peakRuns));
        }

        final double e = evalue(best, lambda, extend, nw);
        lines.add(xn + "\t" + len + "\t" + nw + "\t" + (bestAt + 1) + "\t"
                + Math.round(100.0 * best / nq) + "\t" + Math.round(100 * identityOf(best, nq, k, lambda)) + "\t"
                + best + "\t" + bestRuns + "\t" + String.format(Locale.US, "%.3g", lambda) + "\t" + String.format(Locale.US, "%.2g", e));

        for (int[] r : found) {
            final double ev = evalue(r[2], lambda, extend, nw);
            if (ev < 0.01) {
                hits.add(xn + "\t" + (r[0] + 1) + "\t" + (r[1] + 1) + "\t" + (r[1] - r[0] + 1) + "\t"
                        + (r[3] + 1) + "\t" + Math.round(100.0 * r[2] / nq) + "\t"
                        + Math.round(100 * identityOf(r[2], nq, k, lambda)) + "\t" + r[4] + "\t" + String.format(Locale.US, "%.2g", ev));
            }
        }
        profiles.add(profile);
        profileOf.add(xn + "\t" + stride);
    }

    private void reportScan(String queryName, int window, int nq, int k, int scaled,
            List<String> lines, List<String> hits, List<int[]> profiles, List<String> profileOf) throws IOException {
        System.out.println(hits.isEmpty() ? "No region of any sequence holds the query (E < 0.01)."
                : (hits.size() + " region(s) hold the query. See " + reportfile));
        for (int i = 0; i < Math.min(10, hits.size()); i++) {
            final String[] f = hits.get(i).split("\t");
            System.out.println("   " + f[0] + "  " + f[1] + ".." + f[2] + " (" + f[3] + " bp), best window at " + f[4]
                    + "  containment " + f[5] + "%  identity " + f[6] + "%  matches " + f[7] + "  E=" + f[8]);
        }

        try (BufferedWriter w = newWriter(reportfile)) {
            StringBuilder sr = new StringBuilder(1 << 17);
            sr.append("scan, k=").append(k);
            if (scaled > 1) {
                sr.append(", FracMinHash scaled=").append(scaled).append(" (estimated)");
            }
            sr.append(", query=").append(queryName)
                    .append(", window=").append(window).append(scaled > 1 ? " bp, sketch k-mers of the query=" : " bp, distinct k-mers of the query=").append(nq).append("\n\n");

            sr.append("best window of each sequence\n");
            sr.append("ID\tlength\twindows\tbest at\tcontainment %\tidentity %\tkmers found\tmatches\tby chance\tE\n");
            for (String l : lines) {
                sr.append(l).append("\n");
            }
            sr.append("\n\n");

            sr.append("the regions that hold the query (E < 0.01)\n");
            sr.append("ID\tstart\tend\tlength\tpeak at\tcontainment %\tidentity %\tmatches\tE\n");
            for (String l : hits) {
                sr.append(l).append("\n");
            }
            if (hits.isEmpty()) {
                sr.append("(none)\n");
            }
            sr.append("\n\n");

            sr.append("profile: containment % of the query, by the position of the window\n");
            for (int i = 0; i < profiles.size(); i++) {
                final String[] f = profileOf.get(i).split("\t");
                final int stride = Integer.parseInt(f[1]);
                final int[] p = profiles.get(i);
                sr.append("\n").append(f[0]).append("\nposition\tcontainment %\n");
                for (int n = 0; n < p.length; n++) {
                    sr.append(n * stride + 1).append("\t").append(Math.round(100.0 * p[n] / nq)).append("\n");
                    flush(w, sr);
                }
            }
            w.write(sr.toString());
        }
    }

    /**
     * The bounds of a region.
     *
     * Not those of the windows that hold the query: a window overlapping the
     * query by a little already carries enough of it to count, so those bounds
     * run wide by nearly a window on either side. The bounds wanted are those
     * of the k-mers of the query actually met - but a k-mer or two of it turn
     * up anywhere by chance, and one of those, a thousand bases off, would drag
     * a bound along with it. So the k-mers met are gathered into runs, no more
     * than gap apart, and the bounds are those of the thickest run: the thing
     * itself, and not what chance scattered around it.
     */
    private static int[] region(int[] hit, int n, int lo, int hi, int k, int len, int gap, int peak, int peakAt) {
        int a = 0, b = n;
        while (a < b) {                      // the first k-mer met at or after lo
            final int m = (a + b) >>> 1;
            if (hit[m] < lo) {
                a = m + 1;
            } else {
                b = m;
            }
        }
        int first = -1, last = -1, count = 0;
        int cf = -1, cl = -1, cc = 0;        // the run being walked
        for (int i = a; i < n && hit[i] <= hi; i++) {
            if (cf < 0 || hit[i] - cl > gap) {
                if (cc > count) {
                    first = cf;
                    last = cl;
                    count = cc;
                }
                cf = hit[i];
                cc = 0;
            }
            cl = hit[i];
            cc++;
        }
        if (cc > count) {
            first = cf;
            last = cl;
        }
        if (first < 0) {                     // nothing met: keep the windows
            first = lo;
            last = hi;
        }
        return new int[]{first, Math.min(len - 1, last + k - 1), peak, peakAt};
    }

    private String name(int j) {
        return sname[j].replace(' ', '_').replace(',', '_');
    }

    /**
     * identity of the found part, from the share of the query k-mers that are
     * there
     */
    private static double identityOf(int found, int nq, int k, double lambda) {
        if (nq == 0 || found < 5 || found <= 10 * lambda) {
            return 0;
        }
        return Math.pow((double) found / nq, 1.0 / k);
    }

    /**
     * How many windows holding this many k-mers of the query a sequence of no
     * relation would be expected to yield: the chance of {@code found} or more
     * under one window, times the number of windows walked.
     *
     * The k-mers met by chance come in runs. Chance matches start at a Poisson
     * rate, and a match of k bases is extended by one more base - one more
     * k-mer of the query - with probability {@code extend}, so the length of a
     * run is geometric and the number of k-mers under a window follows the
     * geometric Poisson (Polya-Aeppli) law of mean {@code lambda}. Its upper
     * tail is far heavier than that of a Poisson law of the same mean, which is
     * what a single chance match of a dozen bases beyond k calls for; treating
     * the k-mers as independent made such a match look significant. The law is
     * evaluated by its three-term recurrence, n p(n) = (r (1-t) + 2 t (n-1))
     * p(n-1) - t^2 (n-2) p(n-2), r = lambda (1-t), t = extend. Windows overlap,
     * so this is still the usual first approximation.
     */
    private static double evalue(int found, double lambda, double extend, int windows) {
        if (found <= 0) {
            return windows;
        }
        final double t = extend;
        final double r = lambda * (1 - t);     // the rate at which chance matches start
        final double c0 = r * (1 - t);
        double tail;
        if (found == 1) {
            tail = -Math.expm1(-r);
        } else {
            double a = Math.exp(-r);           // p(n-2), from p(0)
            double b = c0 * a;                 // p(n-1), from p(1)
            for (int n = 2; n <= found; n++) {
                final double c = Math.max(0, ((c0 + 2 * t * (n - 1)) * b - t * t * (n - 2) * a) / n);
                a = b;
                b = c;
            }
            tail = 0;                          // b is p(found): sum it and what follows
            for (int n = found + 1; n < found + 100000; n++) {
                tail += b;
                final double c = Math.max(0, ((c0 + 2 * t * (n - 1)) * b - t * t * (n - 2) * a) / n);
                a = b;
                b = c;
                if (b <= 1e-17 * tail) {
                    break;
                }
            }
        }
        return windows * Math.min(1.0, tail);
    }

    /**
     * A region and the largest number of separate matches a window of it held.
     */
    private static int[] withRuns(int[] region, int runs) {
        final int[] r = Arrays.copyOf(region, region.length + 1);
        r[region.length] = runs;
        return r;
    }

    /**
     * number of k-mers the two sorted sets have in common
     */
    private static int intersection(long[] a, long[] b) {
        int i = 0, j = 0, n = 0;
        while (i < a.length && j < b.length) {
            if (a[i] < b[j]) {
                i++;
            } else if (a[i] > b[j]) {
                j++;
            } else {
                n++;
                i++;
                j++;
            }
        }
        return n;
    }

    /**
     * FracMinHash: the largest hash a k-mer may have and still be kept, for a
     * downsampling factor of {@code scaled}. About one k-mer in {@code scaled}
     * passes, so a genome's k-mer set shrinks by that factor while the
     * containment estimated from the kept fraction stays unbiased (Irber et
     * al., sourmash; Hera, Pierce-Ward & Koslicki, Genome Research, 2023).
     *
     * <p>
     * {@code scaled <= 1} returns -1L (every hash passes) - the exact set, bit
     * for bit what the tool produced before sketching existed.
     */
    static long maxHashForScale(int scaled) {
        return (scaled <= 1) ? -1L : Long.divideUnsigned(-1L, scaled);
    }

    /**
     * SplitMix64 finalizer: spreads a packed k-mer code uniformly over 64 bits.
     */
    private static long mix64(long z) {
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }

    /**
     * The distinct canonical k-mers of a sequence - min(k-mer, its reverse
     * complement) - sorted, so that two of these sets can be intersected by a
     * single walk. A window carrying an ambiguity symbol holds no k-mer.
     */
    private static long[] canonical(String s, int k) {
        return canonical(s, k, -1L);
    }

    /**
     * As above, but keeps only the k-mers whose hash is at most {@code maxHash}
     * (a FracMinHash sketch), so the set of a 2 Gb genome fits in memory. With
     * {@code maxHash == -1L} nothing is filtered and the result is the exact,
     * byte-for-byte set of {@link #canonical(String,int)}.
     */
    private static long[] canonical(String s, int k, long maxHash) {
        final int n = s.length();
        if (n < k || k < 1 || k > 31) {
            return new long[0];
        }
        final boolean sketch = (maxHash != -1L);
        final long mask = (k == 32) ? -1L : ((1L << (2 * k)) - 1);
        final int top = 2 * (k - 1);
        // The exact set needs one slot per position; a sketch keeps only ~1/scaled
        // of them, so it starts small and grows - never allocating the whole-length
        // buffer, which for a 2 Gb genome would be the very 16 GB the sketch avoids.
        long[] buf = sketch ? new long[1024] : new long[n - k + 1];
        int m = 0;
        long fwd = 0, rev = 0;
        int valid = 0;
        for (int i = 0; i < n; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b < 0) {
                valid = 0;
                fwd = 0;
                rev = 0;
                continue;
            }
            fwd = ((fwd << 2) | b) & mask;
            rev = (rev >>> 2) | ((long) (3 - b) << top);   // 3-b complements a<->t and c<->g
            if (valid < k) {
                valid++;
            }
            if (valid == k) {
                final long can = Math.min(fwd, rev);
                if (!sketch) {
                    buf[m++] = can;
                } else if (Long.compareUnsigned(mix64(can), maxHash) <= 0) {
                    if (m == buf.length) {
                        buf = Arrays.copyOf(buf, buf.length * 2);
                    }
                    buf[m++] = can;
                }
            }
        }
        long[] r = Arrays.copyOf(buf, m);
        Arrays.sort(r);
        int d = 0;
        for (int i = 0; i < r.length; i++) {
            if (i == 0 || r[i] != r[i - 1]) {
                r[d++] = r[i];
            }
        }
        return Arrays.copyOf(r, d);
    }

    public static final int COSINE = 1;
    public static final int D2STAR = 2;

    /**
     * Compares the sequences on their k-mer frequencies, which is what the
     * spacing measure of Run() cannot do when the sequences differ widely in
     * length: there, a k-mer only counts once it has been seen TWICE, so a
     * short sequence arrives with a nearly empty vector and a pair of them is
     * declared maximally distant for want of data.
     *
     * Here every occurrence counts. Two statistics are offered:
     *
     * COSINE cos(X,Y) = sum(Xw*Yw) / (sqrt(sum Xw^2) * sqrt(sum Yw^2)) The
     * vectors are compared by their direction alone, so the length of the
     * sequences drops out. But every stretch of DNA has a broadly similar k-mer
     * composition, so the measure keeps a high floor: two unrelated sequences
     * still score around 0.7.
     *
     * D2STAR the count of each k-mer is first centred on what the base
     * composition of that very sequence leads one to expect, Xw~ = (Xw - Ew) /
     * sqrt(Ew), Ew = (L-k+1) * prod p(base), and the centred vectors are then
     * compared as above. What is left is the departure from chance, so
     * unrelated sequences fall to nearly zero and the measure separates far
     * better. This is the d2* of Reinert, Chew, Sun & Waterman (J Comput Biol,
     * 2009) and Song et al. (Brief Bioinform, 2014), written for exactly this
     * problem - sequences of unequal length.
     *
     * The report holds the similarity in percent; the MEGA file holds the
     * angular distance acos(s)/pi, which is a proper metric in [0,1].
     */
    public void RunFrequency(int statistic) throws IOException {
        final int nseq = seq.length;
        final int[][] m2 = countKmers();
        final double[][] e = (statistic == D2STAR) ? new double[nseq][] : null;
        if (e != null) {
            IntStream.range(0, nseq).parallel().forEach(j -> e[j] = expected(seq[j]));
        }
        writeFrequency(sname, m2, e, statistic);
    }

    /**
     * As {@link #RunFrequency(int)}, but the count and expected vectors are
     * built by streaming the files one record at a time, so the sequences are
     * never all in memory - only the small per-record vectors (palette-sized).
     * Byte for byte the same output on the same input.
     */
    public void RunFrequencyStream(String[] files, int statistic) throws IOException {
        final boolean d2 = (statistic == D2STAR);
        final List<String> names = new ArrayList<>();
        final List<int[]> counts = new ArrayList<>();
        final List<double[]> exps = new ArrayList<>();
        StreamingFasta.forEachRecord(files, (name, s) -> {
            names.add(name);
            counts.add(countKmersOne(s));
            if (d2) {
                exps.add(expected(s));
            }
        });
        writeFrequency(names.toArray(new String[0]), counts.toArray(new int[0][]),
                d2 ? exps.toArray(new double[0][]) : null, statistic);
    }

    private void writeFrequency(String[] names, int[][] m2, double[][] eArr, int statistic) throws IOException {
        final int nseq = names.length;
        final int nkmers = pt.size();

        // one unit vector per sequence: the score is then their dot product
        final double[][] v = new double[nseq][];
        IntStream.range(0, nseq).parallel().forEach(j -> {
            final double[] x = new double[nkmers];
            final int[] row = m2[j];
            if (statistic == D2STAR) {
                final double[] e = eArr[j];
                for (int i = 0; i < nkmers; i++) {
                    x[i] = (e[i] > 0) ? (row[i] - e[i]) / Math.sqrt(e[i]) : 0;
                }
            } else {
                for (int i = 0; i < nkmers; i++) {
                    x[i] = row[i];
                }
            }
            double n = 0;
            for (int i = 0; i < nkmers; i++) {
                n += x[i] * x[i];
            }
            n = Math.sqrt(n);
            if (n > 0) {
                for (int i = 0; i < nkmers; i++) {
                    x[i] /= n;
                }
            }
            v[j] = x;
        });

        final double[][] score = new double[nseq][];   // lower triangle
        for (int i = 1; i < nseq; i++) {
            score[i] = new double[i];
        }
        IntStream.range(1, nseq).parallel().forEach(i -> {
            final double[] y = v[i];
            for (int j = 0; j < i; j++) {
                final double[] x = v[j];
                double s = 0;
                for (int n = 0; n < nkmers; n++) {
                    s += x[n] * y[n];
                }
                score[i][j] = s;
            }
        });

        final String[] xname = new String[nseq];
        for (int j = 0; j < nseq; j++) {
            xname[j] = names[j].replace(' ', '_').replace(',', '_');
        }

        try (BufferedWriter w = newWriter(reportfile)) {
            StringBuilder sr = new StringBuilder(1 << 17);

            sr.append("ID/kmer");
            for (int n1 = 0; n1 < nkmers; n1++) {
                sr.append("\t").append(keys[n1]);
            }
            sr.append("\n");
            for (int j = 0; j < nseq; j++) {
                sr.append(xname[j]);
                final int[] row = m2[j];
                for (int n1 = 0; n1 < nkmers; n1++) {
                    sr.append("\t").append(row[n1]);
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n\n");

            sr.append(statistic == D2STAR ? "d2*, %" : "cosine, %");
            for (int j = 0; j < nseq; j++) {
                sr.append("\t").append(xname[j]);
            }
            sr.append("\n");
            for (int i = 0; i < nseq; i++) {
                sr.append(xname[i]);
                for (int n = 0; n < nseq; n++) {
                    if (i == n) {
                        sr.append("\t").append("100");
                    } else {
                        sr.append("\t").append(Math.round(100 * get(score, i, n)));
                    }
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n");
            w.write(sr.toString());
        }

        try (BufferedWriter w = newWriter(megafile)) {
            StringBuilder mg = new StringBuilder(1 << 17);
            mg.append("#mega\n");
            mg.append("!Format DataType=Distance DataFormat=LowerLeft NTaxa=").append(nseq).append(";\n\n");
            for (int j = 0; j < nseq; j++) {
                String s = names[j];
                if (s.length() > 48) {
                    s = s.substring(0, 48);
                }
                mg.append("#").append(s).append("\n");
            }
            mg.append("\n\n");
            for (int i = 1; i < nseq; i++) {
                for (int n = 0; n < i; n++) {
                    mg.append(" ").append(String.format(Locale.US, "%.5f", angular(score[i][n])));
                }
                mg.append("\n");
                flush(w, mg);
            }
            w.write(mg.toString());
        }
    }

    public static final int VECTOR = 3;

    // ── -vector tuning ──────────────────────────────────────────────────────────
    // Ported from TotalRepeats' -vector: a composition measure that sees homology
    // between fragments of very different length. Three coupled parts —
    //   * scale-matched seed windows: a short sequence is compared not to a long
    //     one as a whole (which averages a shared region away) but to the best
    //     WINDOW of the long one at the short one's own scale;
    //   * both strands: the short is also compared reverse-complemented, so an
    //     element inserted the other way round is found;
    //   * a per-length-scale chance floor calibrated from THIS run's own content,
    //     because the null of a composition score depends on the data and on the
    //     scale, and no flat threshold is right for both.
    // The statistic is d2* (counts centred on the base composition), whose null is
    // low and whose floor stays meaningful even on the few-sequence inputs a
    // distance matrix is usually built from; the calibrated floor tightens it
    // further when the input is large enough to measure one.
    private static final int V_CLASS_BASE = 128;   // shortest length class
    private static final int V_CLASS_COUNT = 18;    // sqrt(2)-spaced: 128 … ~46 kb
    private static final double V_WINDOW_MIN_RATIO = 1.5; // window the long only past this ratio
    private static final int V_NULL_WINDOWS = 256;  // windows drawn per class to measure its null
    private static final int V_NULL_DRAW_LIMIT = 80;// draws allowed per usable window
    private static final int V_MIN_NULL_WINDOWS = 32;
    private static final int V_MIN_NULL_PAIRS = 256;// too few screened pairs measures nothing
    private static final int V_NULL_MAX_LEN = 6000; // above this the null is ~0 anyway; skip the costly draw
    private static final int V_NULL_SCREEN_K = 25;  // a shared exact 25-mer keeps a pair OUT of the null
    private static final int V_MIN_NNZ = 8;         // a vector with fewer non-zeros is degenerate
    private static final long V_NULL_SEED = 0x9E3779B97F4A7C15L; // fixed: the run stays deterministic
    private static final double V_LN2 = Math.log(2.0);

    /**
     * Compares the sequences on their k-mer composition in the manner of
     * TotalRepeats' -vector, so that homology is found even between fragments
     * of very different length. Unlike RunFrequency, which compares whole
     * vectors, here a short sequence is matched against the best scale-matched
     * WINDOW of a longer one, on either strand, and a score is believed only
     * when it stands clear of the chance level measured for that length scale.
     *
     * The report holds the similarity in percent; the MEGA file holds the
     * angular distance, and a pair judged no better than chance is set to 0 /
     * distance 1.
     */
    /**
     * -vector reading from a {@link SeqStore}, so the input is streamed into a
     * segmented store rather than an array of Strings: the corpus (and a single
     * record) may exceed 2 GB, and the reverse-complement of each record is
     * built and dropped instead of held for every sequence at once. Byte for
     * byte the result of the array-backed version on the same input; the
     * windows the store slices are the same substrings the old code took from
     * the Strings.
     */
    public void RunVector(SeqStore st) throws IOException {
        final int nseq = st.count();

        final double[] floor = calibrateVectorFloors(st);

        // whole-sequence d2* unit vectors, forward and reverse-complement, reused as
        // the "short side" of every pair and for the comparable-length case
        final double[][] wF = new double[nseq][];
        final double[][] wR = new double[nseq][];
        IntStream.range(0, nseq).parallel().forEach(j -> {
            final String r = st.record(j);
            wF[j] = d2sUnit(r, 0, r.length());
            wR[j] = d2sUnit(revcomp(r), 0, r.length());   // the reverse complement is dropped after
        });

        final double[][] score = new double[nseq][];   // lower triangle
        for (int i = 1; i < nseq; i++) {
            score[i] = new double[i];
        }
        IntStream.range(1, nseq).parallel().forEach(i -> {
            for (int j = 0; j < i; j++) {
                score[i][j] = vectorPair(st, i, j, wF, wR, floor);
            }
        });

        final String[] xname = new String[nseq];
        for (int j = 0; j < nseq; j++) {
            xname[j] = st.name(j).replace(' ', '_').replace(',', '_');
        }

        try (BufferedWriter w = newWriter(reportfile)) {
            StringBuilder sr = new StringBuilder(1 << 16);
            sr.append("vector: scale-matched windows, both strands, chance level calibrated from the input; per-window kernel = d2*, k=").append(kmer).append("\n");
            sr.append("a pair scoring no better than the chance level of its length scale is set to 0\n\n");
            sr.append("similarity %");
            for (int j = 0; j < nseq; j++) {
                sr.append("\t").append(xname[j]);
            }
            sr.append("\n");
            for (int i = 0; i < nseq; i++) {
                sr.append(xname[i]);
                for (int n = 0; n < nseq; n++) {
                    sr.append("\t").append(i == n ? 100 : Math.round(100 * get(score, i, n)));
                }
                sr.append("\n");
                flush(w, sr);
            }
            sr.append("\n");
            w.write(sr.toString());
        }

        try (BufferedWriter w = newWriter(megafile)) {
            StringBuilder mg = new StringBuilder(1 << 16);
            mg.append("#mega\n");
            mg.append("!Format DataType=Distance DataFormat=LowerLeft NTaxa=").append(nseq).append(";\n\n");
            for (int j = 0; j < nseq; j++) {
                String s = st.name(j);
                if (s.length() > 48) {
                    s = s.substring(0, 48);
                }
                mg.append("#").append(s).append("\n");
            }
            mg.append("\n\n");
            for (int i = 1; i < nseq; i++) {
                for (int n = 0; n < i; n++) {
                    final double s = score[i][n];
                    mg.append(" ").append(String.format(Locale.US, "%.5f", s <= 0 ? 1.0 : angular(s)));
                }
                mg.append("\n");
                flush(w, mg);
            }
            w.write(mg.toString());
        }
    }

    /**
     * The score of one pair: the shorter sequence against the best window of
     * the longer at the shorter's scale (or whole-against-whole when their
     * lengths are within {@link #V_WINDOW_MIN_RATIO}), on the better of the two
     * strands. A result at or below the calibrated chance level of the
     * shorter's length scale is returned as 0.
     */
    private double vectorPair(SeqStore st, int i, int j, double[][] wF, double[][] wR, double[] floor) {
        final int shortIdx = st.length(i) <= st.length(j) ? i : j;
        final int longIdx = shortIdx == i ? j : i;
        final int shortLen = (int) st.length(shortIdx);
        final int longLen = (int) st.length(longIdx);
        final int cls = vClass(shortLen);
        final double f = floor[cls];   // -1 when chance was not measurable at this scale

        double best;
        if (longLen >= V_WINDOW_MIN_RATIO * shortLen) {
            final double[] sF = wF[shortIdx], sR = wR[shortIdx];
            final int wlen = shortLen, step = Math.max(1, wlen / 2);
            best = -1;
            for (int off = 0; off + wlen <= longLen; off += step) {
                final double[] win = d2sUnit(st.window(longIdx, off, wlen), 0, wlen);
                final double a = dot(sF, win), b = dot(sR, win);
                final double s = a > b ? a : b;
                if (s > best) {
                    best = s;
                }
            }
        } else {   // comparable length: whole against whole, either strand
            final double a = dot(wF[i], wF[j]), b = dot(wR[i], wF[j]);
            best = a > b ? a : b;
        }
        if (best < 0) {
            best = 0;
        }
        return (f >= 0 && best <= f) ? 0 : best;
    }

    private static double dot(double[] a, double[] b) {
        double s = 0;
        for (int i = 0; i < a.length; i++) {
            s += a[i] * b[i];
        }
        return s;
    }

    /**
     * The sqrt(2)-spaced length class of a sequence; the null scales with log
     * length, not length.
     */
    private static int vClass(int len) {
        if (len <= V_CLASS_BASE) {
            return 0;
        }
        final int c = (int) Math.floor(2.0 * (Math.log(len / (double) V_CLASS_BASE) / V_LN2));
        return Math.max(0, Math.min(V_CLASS_COUNT - 1, c));
    }

    /**
     * The representative (shortest) length of a class.
     */
    private static int vClassLen(int c) {
        return (int) Math.round(V_CLASS_BASE * Math.pow(2, c / 2.0));
    }

    /**
     * Reverse complement; a non-acgt symbol passes through and is never
     * counted.
     */
    private static String revcomp(String s) {
        final int n = s.length();
        final char[] r = new char[n];
        for (int i = 0; i < n; i++) {
            final char c = s.charAt(n - 1 - i);
            r[i] = c == 'a' ? 't' : c == 't' ? 'a' : c == 'c' ? 'g' : c == 'g' ? 'c' : c;
        }
        return new String(r);
    }

    /**
     * The d2* unit vector of s[from,to): the palette counts centred on what the
     * base composition of that same window leads one to expect, (count -
     * E)/sqrt(E), then scaled to unit length so a dot product is the cosine of
     * the centred vectors.
     */
    private double[] d2sUnit(String s, int from, int to) {
        final int nk = pt.size();
        final int[] cnt = new int[nk];
        final int[] tbl = codeIndex;
        final int mask = (1 << (2 * kmer)) - 1;
        int code = 0, valid = 0;
        for (int i = from; i < to; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b < 0) {
                valid = 0;
                continue;
            }
            code = ((code << 2) | b) & mask;
            if (valid < kmer) {
                valid++;
            }
            if (valid == kmer) {
                final int x = tbl[code];
                if (x >= 0) {
                    cnt[x]++;
                }
            }
        }
        final double[] p = new double[4];
        long t = 0;
        for (int i = from; i < to; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b >= 0) {
                p[b]++;
                t++;
            }
        }
        final double[] x = new double[nk];
        if (t == 0) {
            return x;
        }
        for (int i = 0; i < 4; i++) {
            p[i] /= t;
        }
        final double windows = Math.max(0, (to - from) - kmer + 1);
        double norm = 0;
        for (int i = 0; i < nk; i++) {
            final String key = keys[i];
            double e = windows;
            for (int jj = 0; jj < kmer; jj++) {
                e *= p[BASE[key.charAt(jj)]];
            }
            final double xi = (e > 0) ? (cnt[i] - e) / Math.sqrt(e) : 0;
            x[i] = xi;
            norm += xi * xi;
        }
        norm = Math.sqrt(norm);
        if (norm > 0) {
            for (int i = 0; i < nk; i++) {
                x[i] /= norm;
            }
        }
        return x;
    }

    /**
     * Per length class, the highest d2* score that non-homologous windows of
     * that length reach against each other in this run's own content — the
     * chance level of the scale. -1 for a class where chance could not be
     * measured (too little content, or too repetitive to draw enough
     * non-homologous pairs); such a class is left without a cutoff rather than
     * made unusable, so a small input still gets its raw, well-separated
     * scores.
     */
    private double[] calibrateVectorFloors(SeqStore st) {
        final double[] floors = new double[V_CLASS_COUNT];
        Arrays.fill(floors, -1.0);
        final boolean[] need = new boolean[V_CLASS_COUNT];
        for (int j = 0; j < st.count(); j++) {
            need[vClass((int) st.length(j))] = true;
        }
        // sequences long enough to draw a window of each class from
        IntStream.range(0, V_CLASS_COUNT).parallel().forEach(c -> {
            if (!need[c]) {
                return;
            }
            final int wlen = vClassLen(c);
            if (wlen < kmer || wlen > V_NULL_MAX_LEN) {
                return;   // above the cap chance is ~0; leave the class uncapped (raw scores are already low)
            }
            floors[c] = vNullFloor(st, wlen, V_NULL_SEED + c);
        });
        return floors;
    }

    /**
     * The chance level at window length {@code wlen}: draw windows of that
     * length from the input content, keep only pairs sharing no exact 25-mer,
     * and place the level where the largest of that many chance scores is
     * expected, from the centre and the lower half of their distribution (see
     * below; never above the highest score any such pair reaches). Returns -1
     * when too few windows or too few screened pairs can be drawn to measure
     * it.
     */
    private double vNullFloor(SeqStore st, int wlen, long seed) {
        final int[] eligible = eligibleFor(st, wlen);
        if (eligible.length == 0) {
            return -1;
        }
        final java.util.Random rnd = new java.util.Random(seed);
        final double[][] vec = new double[V_NULL_WINDOWS][];
        final long[][] key = new long[V_NULL_WINDOWS][];
        int built = 0;
        final int maxDraws = V_NULL_WINDOWS * V_NULL_DRAW_LIMIT;
        for (int draw = 0; draw < maxDraws && built < V_NULL_WINDOWS; draw++) {
            final int rec = eligible[rnd.nextInt(eligible.length)];
            final int slen = (int) st.length(rec);
            final int off = (slen == wlen) ? 0 : rnd.nextInt(slen - wlen + 1);
            final String win = st.window(rec, off, wlen);   // == the record's substring(off, off+wlen)
            final double[] v = d2sUnit(win, 0, wlen);
            if (nnz(v) < V_MIN_NNZ) {
                continue;
            }
            vec[built] = v;
            key[built] = canonical(win, V_NULL_SCREEN_K);
            built++;
        }
        if (built < V_MIN_NULL_WINDOWS) {
            return -1;
        }
        final double[] sc = new double[built * (built - 1) / 2];
        double best = 0;
        int pairs = 0;
        for (int a = 0; a < built; a++) {
            for (int b = a + 1; b < built; b++) {
                if (intersection(key[a], key[b]) > 0) {
                    continue;   // may be homologous: not part of the null
                }
                final double s = dot(vec[a], vec[b]);
                sc[pairs++] = s;
                if (s > best) {
                    best = s;
                }
            }
        }
        if (pairs < V_MIN_NULL_PAIRS) {
            return -1;
        }
        // The 25-mer screen only removes near-identical windows: copies of one element
        // that have diverged by more than ~10% share no 25-mer, pass it, and would set
        // the "chance" level to their own score. So the level is not taken from the
        // top of the sample, where those homologues sit, but from its centre and its
        // lower half, which they cannot reach: the median plus z times the spread
        // below the median (median - 15.87th percentile is one sigma of a normal
        // sample), z being where the largest of this many chance scores is expected.
        // It never exceeds the largest score seen, so a clean sample is cut as before.
        Arrays.sort(sc, 0, pairs);
        final double med = sc[pairs / 2];
        final double sigma = med - sc[(int) (0.1587 * pairs)];
        final double l2 = 2 * Math.log(pairs);
        final double z = Math.sqrt(l2) - (Math.log(Math.log(pairs)) + Math.log(4 * Math.PI)) / (2 * Math.sqrt(l2));
        return Math.max(0, Math.min(best, med + z * sigma));
    }

    private int[] eligibleFor(SeqStore st, int wlen) {
        int n = 0;
        for (int j = 0; j < st.count(); j++) {
            if (st.length(j) >= wlen) {
                n++;
            }
        }
        final int[] e = new int[n];
        int p = 0;
        for (int j = 0; j < st.count(); j++) {
            if (st.length(j) >= wlen) {
                e[p++] = j;
            }
        }
        return e;
    }

    private static int nnz(double[] v) {
        int n = 0;
        for (double x : v) {
            if (x != 0) {
                n++;
            }
        }
        return n;
    }

    private static double get(double[][] tri, int i, int j) {
        return i > j ? tri[i][j] : tri[j][i];
    }

    /**
     * angular distance: a metric in [0,1], 0 for identical directions
     */
    private static double angular(double s) {
        return Math.acos(Math.max(-1, Math.min(1, s))) / Math.PI;
    }

    /**
     * How often each k-mer of the set is to be expected in this sequence, from
     * its own base composition alone: E = (L-k+1) * p(b1)*...*p(bk).
     */
    private double[] expected(String s) {
        final double[] p = new double[4];
        long t = 0;
        for (int i = 0, e = s.length(); i < e; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b >= 0) {
                p[b]++;
                t++;
            }
        }
        final double[] e = new double[pt.size()];
        if (t == 0) {
            return e;
        }
        for (int i = 0; i < 4; i++) {
            p[i] /= t;
        }
        final double windows = Math.max(0, s.length() - kmer + 1);
        for (int i = 0; i < e.length; i++) {
            final String key = keys[i];
            double q = windows;
            for (int j = 0; j < kmer; j++) {
                q *= p[BASE[key.charAt(j)]];
            }
            e[i] = q;
        }
        return e;
    }

    public void RunKmerCounter(int d) throws IOException {
        final int nseq = seq.length;
        final int nkmers = pt.size();

        final int[][] m2 = countKmers();
        for (int j = 0; j < nseq; j++) {
            final int[] row = m2[j];
            final int len = seq[j].length();
            for (int n1 = 0; n1 < nkmers; n1++) {
                row[n1] = (len / (1 + row[n1]));
            }
        }

        try (BufferedWriter w = newWriter(reportfile)) {
            StringBuilder sr = new StringBuilder(1 << 17);
            sr.append("kmerID/AvDist");
            if (d == 1) {
                for (int j = 0; j < nseq; j++) {
                    sr.append("\t").append(sname[j]);
                }
                sr.append("\n");

                for (int n1 = 0; n1 < nkmers; n1++) {
                    sr.append(keys[n1]);
                    for (int j = 0; j < nseq; j++) {
                        sr.append("\t").append(m2[j][n1]);
                    }
                    sr.append("\n");
                    flush(w, sr);
                }
            } else {
                sr.append("\t").append(fPath).append("\n");
                for (int n1 = 0; n1 < nkmers; n1++) {
                    sr.append(keys[n1]);
                    long s = 0;
                    for (int j = 0; j < nseq; j++) {
                        s += m2[j][n1];
                    }
                    sr.append("\t").append(s / nseq).append("\n");
                    flush(w, sr);
                }
            }
            sr.append("\n");
            w.write(sr.toString());
        }
    }

    /**
     * Counts how often each k-mer of the set occurs in each sequence.
     *
     * The sets hold nothing but a/c/g/t, so a window carrying any other symbol
     * (n, r, y, ...) can never be one of them: two bits per base and a rolling
     * window turn the lookup into a single array read, in place of a
     * substring() and a hash lookup at every position.
     */
    /**
     * The palette k-mer counts of one sequence - one row of
     * {@link #countKmers()}.
     */
    private int[] countKmersOne(String s) {
        final int nkmers = pt.size();
        final int[] row = new int[nkmers];
        final int[] tbl = codeIndex;
        if (tbl == null) {
            for (int i = 0, e = s.length() - kmer + 1; i < e; i++) {
                final Integer x = pt.get(s.substring(i, i + kmer));
                if (x != null) {
                    row[x]++;
                }
            }
            return row;
        }
        final int mask = (1 << (2 * kmer)) - 1;
        final int k = kmer;
        int code = 0, valid = 0;
        for (int i = 0, e = s.length(); i < e; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b < 0) {
                valid = 0;
                continue;
            }
            code = ((code << 2) | b) & mask;
            if (valid < k) {
                valid++;
            }
            if (valid == k) {
                final int x = tbl[code];
                if (x >= 0) {
                    row[x]++;
                }
            }
        }
        return row;
    }

    /**
     * The counts of the spacing features of one sequence: as
     * {@link #countKmersOne}, but a word and its reverse complement fall in the
     * same column.
     */
    private int[] countSpacingOne(String s) {
        final int[] tbl = spIndex;
        if (tbl == null) {
            return countKmersOne(s);
        }
        final int[] row = new int[spKeys.length];
        final int mask = (1 << (2 * kmer)) - 1;
        final int k = kmer;
        int code = 0, valid = 0;
        for (int i = 0, e = s.length(); i < e; i++) {
            final char c = s.charAt(i);
            final int b = (c < 128) ? BASE[c] : -1;
            if (b < 0) {
                valid = 0;
                continue;
            }
            code = ((code << 2) | b) & mask;
            if (valid < k) {
                valid++;
            }
            if (valid == k) {
                final int x = tbl[code];
                if (x >= 0) {
                    row[x]++;
                }
            }
        }
        return row;
    }

    private int[][] countKmers() {
        final int nseq = seq.length;
        final int nkmers = pt.size();
        final int[][] m2 = new int[nseq][nkmers];
        final int[] tbl = codeIndex;

        if (tbl == null) {   // a set outside acgt: keep the general way
            IntStream.range(0, nseq).parallel().forEach(j -> {
                final String s = seq[j];
                final int[] row = m2[j];
                for (int i = 0, e = s.length() - kmer + 1; i < e; i++) {
                    final Integer x = pt.get(s.substring(i, i + kmer));
                    if (x != null) {
                        row[x]++;
                    }
                }
            });
            return m2;
        }

        final int mask = (1 << (2 * kmer)) - 1;
        final int k = kmer;
        IntStream.range(0, nseq).parallel().forEach(j -> {
            final String s = seq[j];
            final int[] row = m2[j];
            int code = 0;
            int valid = 0;
            for (int i = 0, e = s.length(); i < e; i++) {
                final char c = s.charAt(i);
                final int b = (c < 128) ? BASE[c] : -1;
                if (b < 0) {           // an ambiguity symbol: no k-mer of the set can match here
                    valid = 0;
                    continue;
                }
                code = ((code << 2) | b) & mask;
                if (valid < k) {
                    valid++;
                }
                if (valid == k) {
                    final int x = tbl[code];
                    if (x >= 0) {
                        row[x]++;
                    }
                }
            }
        });
        return m2;
    }

    /**
     * The column of every k-mer of the set, by its 2-bit code; -1 everywhere
     * else.
     */
    private int[] buildCodeIndex() {
        if (kmer < 1 || kmer > 15) {
            return null;
        }
        final int[] tbl = new int[1 << (2 * kmer)];
        Arrays.fill(tbl, -1);
        for (Map.Entry<String, Integer> entry : pt.entrySet()) {
            final String s = entry.getKey();
            if (s.length() != kmer) {
                return null;
            }
            int code = 0;
            for (int i = 0; i < kmer; i++) {
                final char c = s.charAt(i);
                final int b = (c < 128) ? BASE[c] : -1;
                if (b < 0) {
                    return null;   // not a plain acgt set
                }
                code = (code << 2) | b;
            }
            tbl[code] = entry.getValue();
        }
        return tbl;
    }

    private static BufferedWriter newWriter(String file) throws IOException {
        return Files.newBufferedWriter(Path.of(file), StandardCharsets.UTF_8);
    }

    /**
     * Hands the report to the writer as it grows, so that no full copy of it is
     * ever held.
     */
    private static void flush(BufferedWriter w, StringBuilder sb) throws IOException {
        if (sb.length() >= (1 << 16)) {
            w.write(sb.toString());
            sb.setLength(0);
        }
    }

    /**
     * a/c/g/t on two bits, -1 for every other symbol.
     */
    private static final int[] BASE = new int[128];

    static {
        Arrays.fill(BASE, -1);
        BASE['a'] = 0;
        BASE['c'] = 1;
        BASE['g'] = 2;
        BASE['t'] = 3;
    }

    private final HashMap<String, Integer> pt;
    private final String[] keys;
    private final int[] codeIndex;
    private int[] spIndex;      // 2-bit code -> spacing feature (a word and its reverse complement share one)
    private String[] spKeys;    // the names of the spacing features
    private final int kmer;
    private final int mod;
    private String fPath = "";
    private String tag = "";
    private String suffix;
    private String basePath;
    private String reportfile;
    private String megafile;
    private String[] seq;
    private String[] sname;
    private final int minkmer = 1;

    private static final String[] RE4MERS = {"aatt", "ttaa", "acgt", "agct", "tgca", "tcga", "catg", "ctag", "gatc", "gtac", "ccgg", "ggcc"};
    private static final String[] RE6MERS = {"aaattt", "tttaaa", "aatatt", "taatta", "attaat", "ttataa", "atgcat", "tcatga", "ttgcaa", "acatgt", "aagctt", "catatg", "tgatca", "gaattc", "tctaga", "agatct", "tgtaca", "caattg", "cttaag", "agtact", "ccatgg", "tagcta", "actagt", "gatatc", "ctatag", "tggcca", "ctgcag", "gtatac", "gcatgc", "gttaac", "cagctg", "aggcct", "gagctc", "ggatcc", "cctagg", "gtgcac", "ggtacc", "cacgtg", "gctagc", "ttcgaa", "aacgtt", "tacgta", "ctcgag", "atcgat", "gggccc", "tccgga", "gacgtc", "agcgct", "tgcgca", "accggt", "ggcgcc", "cccggg", "gtcgac", "gccggc", "gcgcgc", "ccgcgg", "tcgcga", "cggccg", "acgcgt", "cgatcg", "cgtacg"};

    private static final String[] RE8MERS = {"cctcgagg", "ggcgcgcc", "ggccggcc", "cgcgcgcg", "gcggccgc", "ttaattaa", "gtttaaac", "gcgatcgc", "cgtcgacg", "gcccgggc", "cgccggcg", "cctgcagg", "atttaaat", "tcgcgcga", "tcgtacga", "cgatatcg", "tacgcgta", "ccgtacgg", "gcgtacgc", "ccgatcgg", "cggtaccg", "atcgcgat", "ttcgcgaa", "gtcgcgac", "cgttaacg", "cgaattcg", "cgtatacg", "tcgatcga", "acgcgcgt", "cgctagcg", "acgatcgt", "tcggccga", "accgcggt", "cggatccg", "aacgcgtt", "gacgcgtc", "ctcgcgag", "agcgcgct", "acggccgt", "ggtcgacc", "cgagctcg", "tccgcgga", "gtccggac", "acgtacgt", "ccgcgcgg", "gaccggtc", "taccggta", "cgggcccg", "tgcgcgca", "cgcatgcg", "tagcgcta", "cccgcggg", "gatcgatc", "ttgcgcaa", "gggcgccc", "cggcgccg", "cgtgcacg", "atccggat", "gtacgtac", "aaccggtt", "ctacgtag", "gagcgctc", "ggacgtcc", "gctcgagc", "gttcgaac", "atgcgcat", "agtcgact", "gccgcggc", "cacgcgtg", "ccggccgg", "gtgcgcac", "ctgcgcag", "caccggtg", "agccggct", "tatcgata", "gaacgttc", "aggcgcct", "tgtcgaca", "tgccggca", "acccgggt", "taacgtta", "tggcgcca", "caacgttg", "catcgatg", "agacgtct", "cagcgctg", "tctcgaga", "ggggcccc", "tcccggga", "actcgagt", "ctccggag", "aatcgatt", "ttccggaa", "ttacgtaa", "cttcgaag", "aagcgctt", "ccccgggg", "tgacgtca", "gcacgtgc", "ggctagcc", "gcctaggc", "gactagtc", "ggttaacc", "gggtaccc", "tttcgaaa", "attcgaat", "ggtgcacc", "gggatccc", "ggtatacc", "gtctagac", "tcacgtga", "gagtactc", "agggccct", "gtagctac", "ctggccag", "taggccta", "atacgtat", "ccaattgg", "caggcctg", "ccctaggg", "cctatagg", "ggcatgcc", "ggatatcc", "gctatagc", "acctaggt", "gcttaagc", "cagtactg", "gttgcaac", "ccacgtgg", "ctgtacag", "aggtacct", "ctagctag", "aaacgttt", "tgggccca", "ggagctcc", "ccagctgg", "gtcatgac", "gtgtacac", "cccatggg", "gaggcctc", "ccttaagg", "agtgcact", "gatgcatc", "gacatgtc", "gcaattgc", "cactagtg", "gttataac", "acacgtgt", "ccatatgg", "gattaatc", "tagtacta", "gcatatgc", "acagctgt", "agctagct", "tactagta", "gtaattac", "ctgatcag", "gtggccac", "gccatggc", "tcctagga", "accatggt", "taagctta", "tggtacca", "caagcttg", "tagatcta", "aggatcct", "tgctagca", "cttgcaag", "ctaattag", "gtgatcac", "cagatctg", "ctctagag", "agcatgct", "gagatctc", "ttagctaa", "tggatcca", "gatatatc", "agagctct", "ctatatag", "cttataag", "actatagt", "ggaattcc", "atagctat", "tcagctga", "tcaattga", "agtatact", "ctcatgag", "gctgcagc", "atggccat", "tctataga", "gcagctgc", "acttaagt", "tgagctca", "aaggcctt", "agttaact", "atctagat", "gtatatac", "actgcagt", "gaagcttc", "aactagtt", "tgtgcaca", "acaattgt", "agatatct", "cattaatg", "ttggccaa", "ctttaaag", "ttgtacaa", "tgttaaca", "aagtactt", "tcttaaga", "tctgcaga", "attgcaat", "ttgatcaa", "tgatatca", "cacatgtg", "tcatatga", "atcatgat", "tgtataca", "atgatcat", "atgtacat", "gaatattc", "tccatgga", "tatgcata", "gaaatttc", "tacatgta", "caatattg", "aagatctt", "aacatgtt", "tgcatgca", "tgaattca", "catgcatg", "aatgcatt", "acatatgt", "catatatg", "ttcatgaa", "ttctagaa", "agaattct", "aaagcttt", "taatatta", "caaatttg", "attataat", "ttatataa", "tattaata", "tttgcaaa", "ataattat", "taaattta", "aattaatt", "aatatatt", "tttataaa", "aaaatttt", "aaatattt", "ttttaaaa"};
    private static final String[] RE10MERS = {"aaattaattt", "aaatcgattt", "aaatgcattt", "aaacatgttt", "aaactagttt", "aaaccggttt", "aaacgcgttt", "aaagatcttt", "aaagtacttt", "aaagcgcttt", "aaaggccttt", "aattataatt", "aatttaaatt", "aattcgaatt", "aattgcaatt", "aataattatt", "aatacgtatt", "aatagctatt", "aatgatcatt", "aatgtacatt", "aatgcgcatt", "aatggccatt", "aatcatgatt", "aatctagatt", "aatccggatt", "aatcgcgatt", "aaccatggtt", "aacctaggtt", "aacccgggtt", "aaccgcggtt", "aacaattgtt", "aacatatgtt", "aacacgtgtt", "aacagctgtt", "aactatagtt", "aacttaagtt", "aactcgagtt", "aactgcagtt", "aacgatcgtt", "aacgtacgtt", "aacggccgtt", "aaggatcctt", "aaggtacctt", "aaggcgcctt", "aagggccctt", "aagaattctt", "aagatatctt", "aagacgtctt", "aagagctctt", "aagtatactt", "aagttaactt", "aagtcgactt", "aagtgcactt", "aagcatgctt", "aagctagctt", "aagccggctt", "ataaatttat", "ataatattat", "ataacgttat", "ataagcttat", "atacatgtat", "atactagtat", "ataccggtat", "atacgcgtat", "atagatctat", "atagtactat", "atagcgctat", "ataggcctat", "atttataaat", "attttaaaat", "atttcgaaat", "atttgcaaat", "attaattaat", "attacgtaat", "attagctaat", "attgatcaat", "attgtacaat", "attgcgcaat", "attggccaat", "attcatgaat", "attctagaat", "attccggaat", "attcgcgaat", "atccatggat", "atcctaggat", "atcccgggat", "atccgcggat", "atcaattgat", "atcatatgat", "atcacgtgat", "atcagctgat", "atctatagat", "atcttaagat", "atctcgagat", "atctgcagat", "atcgatcgat", "atcgtacgat", "atcggccgat", "atggatccat", "atggtaccat", "atggcgccat", "atgggcccat", "atgaattcat", "atgatatcat", "atgacgtcat", "atgagctcat", "atgtatacat", "atgttaacat", "atgtcgacat", "atgtgcacat", "atgcatgcat", "atgctagcat", "atgccggcat", "acaaatttgt", "acaatattgt", "acaacgttgt", "acaagcttgt", "acattaatgt", "acatcgatgt", "acatgcatgt", "acagatctgt", "acagtactgt", "acagcgctgt", "acaggcctgt", "acttataagt", "actttaaagt", "acttcgaagt", "acttgcaagt", "actaattagt", "actacgtagt", "actagctagt", "actgatcagt", "actgtacagt", "actgcgcagt", "actggccagt", "actcatgagt", "actctagagt", "actccggagt", "actcgcgagt", "acccatgggt", "accctagggt", "accccggggt", "acccgcgggt", "accaattggt", "accatatggt", "accacgtggt", "accagctggt", "acctataggt", "accttaaggt", "acctcgaggt", "acctgcaggt", "accgatcggt", "accgtacggt", "accggccggt", "acggatccgt", "acggtaccgt", "acggcgccgt", "acgggcccgt", "acgaattcgt", "acgatatcgt", "acgacgtcgt", "acgagctcgt", "acgtatacgt", "acgttaacgt", "acgtcgacgt", "acgtgcacgt", "acgcatgcgt", "acgctagcgt", "acgccggcgt", "agaaatttct", "agaatattct", "agaacgttct", "agaagcttct", "agattaatct", "agatcgatct", "agatgcatct", "agacatgtct", "agactagtct", "agaccggtct", "agacgcgtct", "agttataact", "agtttaaact", "agttcgaact", "agttgcaact", "agtaattact", "agtacgtact", "agtagctact", "agtgatcact", "agtgtacact", "agtgcgcact", "agtggccact", "agtcatgact", "agtctagact", "agtccggact", "agtcgcgact", "agccatggct", "agcctaggct", "agcccgggct", "agccgcggct", "agcaattgct", "agcatatgct", "agcacgtgct", "agcagctgct", "agctatagct", "agcttaagct", "agctcgagct", "agctgcagct", "agcgatcgct", "agcgtacgct", "agcggccgct", "agggatccct", "agggtaccct", "agggcgccct", "aggggcccct", "aggaattcct", "aggatatcct", "aggacgtcct", "aggagctcct", "aggtatacct", "aggttaacct", "aggtcgacct", "aggtgcacct", "aggcatgcct", "aggctagcct", "aggccggcct", "ttaaatttaa", "ttaatattaa", "ttaacgttaa", "ttaagcttaa", "ttattaataa", "ttatcgataa", "ttatgcataa", "ttacatgtaa", "ttactagtaa", "ttaccggtaa", "ttacgcgtaa", "ttagatctaa", "ttagtactaa", "ttagcgctaa", "ttaggcctaa", "tttaattaaa", "tttacgtaaa", "tttagctaaa", "tttgatcaaa", "tttgtacaaa", "tttgcgcaaa", "tttggccaaa", "tttcatgaaa", "tttctagaaa", "tttccggaaa", "tttcgcgaaa", "ttccatggaa", "ttcctaggaa", "ttcccgggaa", "ttccgcggaa", "ttcaattgaa", "ttcatatgaa", "ttcacgtgaa", "ttcagctgaa", "ttctatagaa", "ttcttaagaa", "ttctcgagaa", "ttctgcagaa", "ttcgatcgaa", "ttcgtacgaa", "ttcggccgaa", "ttggatccaa", "ttggtaccaa", "ttggcgccaa", "ttgggcccaa", "ttgaattcaa", "ttgatatcaa", "ttgacgtcaa", "ttgagctcaa", "ttgtatacaa", "ttgttaacaa", "ttgtcgacaa", "ttgtgcacaa", "ttgcatgcaa", "ttgctagcaa", "ttgccggcaa", "taaaatttta", "taaatattta", "taaacgttta", "taaagcttta", "taattaatta", "taatcgatta", "taatgcatta", "taacatgtta", "taactagtta", "taaccggtta", "taacgcgtta", "taagatctta", "taagtactta", "taagcgctta", "taaggcctta", "tattataata", "tatttaaata", "tattcgaata", "tattgcaata", "tatgatcata", "tatgtacata", "tatgcgcata", "tatggccata", "tatcatgata", "tatctagata", "tatccggata", "tatcgcgata", "taccatggta", "tacctaggta", "tacccgggta", "taccgcggta", "tacaattgta", "tacatatgta", "tacacgtgta", "tacagctgta", "tactatagta", "tacttaagta", "tactcgagta", "tactgcagta", "tacgatcgta", "tacgtacgta", "tacggccgta", "taggatccta", "taggtaccta", "taggcgccta", "tagggcccta", "tagaattcta", "tagatatcta", "tagacgtcta", "tagagctcta", "tagtatacta", "tagttaacta", "tagtcgacta", "tagtgcacta", "tagcatgcta", "tagctagcta", "tagccggcta", "tgaaatttca", "tgaatattca", "tgaacgttca", "tgaagcttca", "tgattaatca", "tgatcgatca", "tgatgcatca", "tgacatgtca", "tgactagtca", "tgaccggtca", "tgacgcgtca", "tgagatctca", "tgagtactca", "tgagcgctca", "tgaggcctca", "tgttataaca", "tgtttaaaca", "tgttcgaaca", "tgttgcaaca", "tgtaattaca", "tgtacgtaca", "tgtagctaca", "tgtcatgaca", "tgtctagaca", "tgtccggaca", "tgtcgcgaca", "tgccatggca", "tgcctaggca", "tgcccgggca", "tgccgcggca", "tgcaattgca", "tgcatatgca", "tgcacgtgca", "tgcagctgca", "tgctatagca", "tgcttaagca", "tgctcgagca", "tgctgcagca", "tgcgatcgca", "tgcgtacgca", "tgcggccgca", "tgggatccca", "tgggtaccca", "tgggcgccca", "tggggcccca", "tggaattcca", "tggatatcca", "tggacgtcca", "tggagctcca", "tggtatacca", "tggttaacca", "tggtcgacca", "tggtgcacca", "tggcatgcca", "tggctagcca", "tggccggcca", "tcaaatttga", "tcaatattga", "tcaacgttga", "tcaagcttga", "tcattaatga", "tcatcgatga", "tcatgcatga", "tcacatgtga", "tcactagtga", "tcaccggtga", "tcacgcgtga", "tcagatctga", "tcagtactga", "tcagcgctga", "tcaggcctga", "tcttataaga", "tctttaaaga", "tcttcgaaga", "tcttgcaaga", "tctaattaga", "tctacgtaga", "tctagctaga", "tctgatcaga", "tctgtacaga", "tctgcgcaga", "tctggccaga", "tcccatggga", "tccctaggga", "tccccgggga", "tcccgcggga", "tccaattgga", "tccatatgga", "tccacgtgga", "tccagctgga", "tcctatagga", "tccttaagga", "tcctcgagga", "tcctgcagga", "tccgatcgga", "tccgtacgga", "tccggccgga", "tcggatccga", "tcggtaccga", "tcggcgccga", "tcgggcccga", "tcgaattcga", "tcgatatcga", "tcgacgtcga", "tcgagctcga", "tcgtatacga", "tcgttaacga", "tcgtcgacga", "tcgtgcacga", "tcgcatgcga", "tcgctagcga", "tcgccggcga", "ccaaatttgg", "ccaatattgg", "ccaacgttgg", "ccaagcttgg", "ccattaatgg", "ccatcgatgg", "ccatgcatgg", "ccacatgtgg", "ccactagtgg", "ccaccggtgg", "ccacgcgtgg", "ccagatctgg", "ccagtactgg", "ccagcgctgg", "ccaggcctgg", "ccttataagg", "cctttaaagg", "ccttcgaagg", "ccttgcaagg", "cctaattagg", "cctacgtagg", "cctagctagg", "cctgatcagg", "cctgtacagg", "cctgcgcagg", "cctggccagg", "cctcatgagg", "cctctagagg", "cctccggagg", "cctcgcgagg", "cccaattggg", "cccatatggg", "cccacgtggg", "cccagctggg", "ccctataggg", "cccttaaggg", "ccctcgaggg", "ccctgcaggg", "cccgatcggg", "cccgtacggg", "cccggccggg", "ccggatccgg", "ccggtaccgg", "ccggcgccgg", "ccgggcccgg", "ccgaattcgg", "ccgatatcgg", "ccgacgtcgg", "ccgagctcgg", "ccgtatacgg", "ccgttaacgg", "ccgtcgacgg", "ccgtgcacgg", "ccgcatgcgg", "ccgctagcgg", "ccgccggcgg", "caaaattttg", "caaatatttg", "caaacgtttg", "caaagctttg", "caattaattg", "caatcgattg", "caatgcattg", "caacatgttg", "caactagttg", "caaccggttg", "caacgcgttg", "caagatcttg", "caagtacttg", "caagcgcttg", "caaggccttg", "cattataatg", "catttaaatg", "cattcgaatg", "cattgcaatg", "cataattatg", "catacgtatg", "catagctatg", "catgatcatg", "catgtacatg", "catgcgcatg", "catggccatg", "catcatgatg", "catctagatg", "catccggatg", "catcgcgatg", "caccatggtg", "cacctaggtg", "cacccgggtg", "caccgcggtg", "cactatagtg", "cacttaagtg", "cactcgagtg", "cactgcagtg", "cacgatcgtg", "cacgtacgtg", "cacggccgtg", "caggatcctg", "caggtacctg", "caggcgcctg", "cagggccctg", "cagaattctg", "cagatatctg", "cagacgtctg", "cagagctctg", "cagtatactg", "cagttaactg", "cagtcgactg", "cagtgcactg", "cagcatgctg", "cagctagctg", "cagccggctg", "ctaaatttag", "ctaatattag", "ctaacgttag", "ctaagcttag", "ctattaatag", "ctatcgatag", "ctatgcatag", "ctacatgtag", "ctactagtag", "ctaccggtag", "ctacgcgtag", "ctagatctag", "ctagtactag", "ctagcgctag", "ctaggcctag", "ctttataaag", "cttttaaaag", "ctttcgaaag", "ctttgcaaag", "cttaattaag", "cttacgtaag", "cttagctaag", "cttgatcaag", "cttgtacaag", "cttgcgcaag", "cttggccaag", "cttcatgaag", "cttctagaag", "cttccggaag", "cttcgcgaag", "ctccatggag", "ctcctaggag", "ctcccgggag", "ctccgcggag", "ctcaattgag", "ctcatatgag", "ctcacgtgag", "ctcagctgag", "ctcgatcgag", "ctcgtacgag", "ctcggccgag", "ctggatccag", "ctggtaccag", "ctggcgccag", "ctgggcccag", "ctgaattcag", "ctgatatcag", "ctgacgtcag", "ctgagctcag", "ctgtatacag", "ctgttaacag", "ctgtcgacag", "ctgtgcacag", "ctgcatgcag", "ctgctagcag", "ctgccggcag", "cgaaatttcg", "cgaatattcg", "cgaacgttcg", "cgaagcttcg", "cgattaatcg", "cgatcgatcg", "cgatgcatcg", "cgacatgtcg", "cgactagtcg", "cgaccggtcg", "cgacgcgtcg", "cgagatctcg", "cgagtactcg", "cgagcgctcg", "cgaggcctcg", "cgttataacg", "cgtttaaacg", "cgttcgaacg", "cgttgcaacg", "cgtaattacg", "cgtacgtacg", "cgtagctacg", "cgtgatcacg", "cgtgtacacg", "cgtgcgcacg", "cgtggccacg", "cgtcatgacg", "cgtctagacg", "cgtccggacg", "cgtcgcgacg", "cgccatggcg", "cgcctaggcg", "cgcccgggcg", "cgccgcggcg", "cgcaattgcg", "cgcatatgcg", "cgcacgtgcg", "cgcagctgcg", "cgctatagcg", "cgcttaagcg", "cgctcgagcg", "cgctgcagcg", "cgggatcccg", "cgggtacccg", "cgggcgcccg", "cggggccccg", "cggaattccg", "cggatatccg", "cggacgtccg", "cggagctccg", "cggtataccg", "cggttaaccg", "cggtcgaccg", "cggtgcaccg", "cggcatgccg", "cggctagccg", "cggccggccg", "ggaaatttcc", "ggaatattcc", "ggaacgttcc", "ggaagcttcc", "ggattaatcc", "ggatcgatcc", "ggatgcatcc", "ggacatgtcc", "ggactagtcc", "ggaccggtcc", "ggacgcgtcc", "ggagatctcc", "ggagtactcc", "ggagcgctcc", "ggaggcctcc", "ggttataacc", "ggtttaaacc", "ggttcgaacc", "ggttgcaacc", "ggtaattacc", "ggtacgtacc", "ggtagctacc", "ggtgatcacc", "ggtgtacacc", "ggtgcgcacc", "ggtggccacc", "ggtcatgacc", "ggtctagacc", "ggtccggacc", "ggtcgcgacc", "ggccatggcc", "ggcctaggcc", "ggcccgggcc", "ggccgcggcc", "ggcaattgcc", "ggcatatgcc", "ggcacgtgcc", "ggcagctgcc", "ggctatagcc", "ggcttaagcc", "ggctcgagcc", "ggctgcagcc", "ggcgatcgcc", "ggcgtacgcc", "ggcggccgcc", "gggaattccc", "gggatatccc", "gggacgtccc", "gggagctccc", "gggtataccc", "gggttaaccc", "gggtcgaccc", "gggtgcaccc", "gggcatgccc", "gggctagccc", "gggccggccc", "gaaaattttc", "gaaatatttc", "gaaacgtttc", "gaaagctttc", "gaattaattc", "gaatcgattc", "gaatgcattc", "gaacatgttc", "gaactagttc", "gaaccggttc", "gaacgcgttc", "gaagatcttc", "gaagtacttc", "gaagcgcttc", "gaaggccttc", "gattataatc", "gatttaaatc", "gattcgaatc", "gattgcaatc", "gataattatc", "gatacgtatc", "gatagctatc", "gatgatcatc", "gatgtacatc", "gatgcgcatc", "gatggccatc", "gatcatgatc", "gatctagatc", "gatccggatc", "gatcgcgatc", "gaccatggtc", "gacctaggtc", "gacccgggtc", "gaccgcggtc", "gacaattgtc", "gacatatgtc", "gacacgtgtc", "gacagctgtc", "gactatagtc", "gacttaagtc", "gactcgagtc", "gactgcagtc", "gacgatcgtc", "gacgtacgtc", "gacggccgtc", "gaggatcctc", "gaggtacctc", "gaggcgcctc", "gagggccctc", "gagtatactc", "gagttaactc", "gagtcgactc", "gagtgcactc", "gagcatgctc", "gagctagctc", "gagccggctc", "gtaaatttac", "gtaatattac", "gtaacgttac", "gtaagcttac", "gtattaatac", "gtatcgatac", "gtatgcatac", "gtacatgtac", "gtactagtac", "gtaccggtac", "gtacgcgtac", "gtagatctac", "gtagtactac", "gtagcgctac", "gtaggcctac", "gtttataaac", "gttttaaaac", "gtttcgaaac", "gtttgcaaac", "gttaattaac", "gttacgtaac", "gttagctaac", "gttgatcaac", "gttgtacaac", "gttgcgcaac", "gttggccaac", "gttcatgaac", "gttctagaac", "gttccggaac", "gttcgcgaac", "gtccatggac", "gtcctaggac", "gtcccgggac", "gtccgcggac", "gtcaattgac", "gtcatatgac", "gtcacgtgac", "gtcagctgac", "gtctatagac", "gtcttaagac", "gtctcgagac", "gtctgcagac", "gtcgatcgac", "gtcgtacgac", "gtcggccgac", "gtggatccac", "gtggtaccac", "gtggcgccac", "gtgggcccac", "gtgaattcac", "gtgatatcac", "gtgacgtcac", "gtgagctcac", "gtgcatgcac", "gtgctagcac", "gtgccggcac", "gcaaatttgc", "gcaatattgc", "gcaacgttgc", "gcaagcttgc", "gcattaatgc", "gcatcgatgc", "gcatgcatgc", "gcacatgtgc", "gcactagtgc", "gcaccggtgc", "gcacgcgtgc", "gcagatctgc", "gcagtactgc", "gcagcgctgc", "gcaggcctgc", "gcttataagc", "gctttaaagc", "gcttcgaagc", "gcttgcaagc", "gctaattagc", "gctacgtagc", "gctagctagc", "gctgatcagc", "gctgtacagc", "gctgcgcagc", "gctggccagc", "gctcatgagc", "gctctagagc", "gctccggagc", "gctcgcgagc", "gcccatgggc", "gccctagggc", "gccccggggc", "gcccgcgggc", "gccaattggc", "gccatatggc", "gccacgtggc", "gccagctggc", "gcctataggc", "gccttaaggc", "gcctcgaggc", "gcctgcaggc", "gccgatcggc", "gccgtacggc", "gccggccggc", "gcggatccgc", "gcggtaccgc", "gcggcgccgc", "gcgggcccgc", "gcgaattcgc", "gcgatatcgc", "gcgacgtcgc", "gcgagctcgc", "gcgtatacgc", "gcgttaacgc", "gcgtcgacgc", "gcgtgcacgc"};

    private static final String[] K4MERS = {"tttc", "gaaa", "attt", "aaat", "tttg", "caaa", "tgaa", "aatt", "ttct", "aaca", "ttca", "agaa", "tgtt", "tctt", "ttga", "aaga", "tcaa", "ttgt", "aaac", "gttt", "acaa", "aaag", "tgat", "cttt", "atca", "tcat", "atga", "catt", "acat", "attg", "caat", "atgt", "aatg", "attc", "gaat", "gatt", "aatc", "tatt", "aata", "atct", "ttcc", "agtt", "aact", "gaag", "cttc", "agat", "aagt", "ggaa", "actt", "tgca", "ccaa", "tcca", "ttgg", "tgga", "ttgc", "gcaa", "gatg", "catc", "tcct", "caac", "gttg", "agga", "caag", "tgta", "aagg", "taca", "agca", "tatc", "cctt", "gata", "tgct", "acca", "atgg", "tggt", "atac", "gtat", "cttg", "ccat", "cata", "gtga", "tatg", "tctg", "gctt", "ttat", "catg", "gcat", "gttc", "ggat", "agag", "ctct", "ataa", "atcc", "aagc", "ggtt", "gaac", "agct", "caga", "aacc", "atgc", "ttac", "gtaa", "ctga", "tcag", "ctgt", "atag", "ctat", "agta", "aggt", "tact", "ctcc", "atta", "taat", "tgtc", "tcac", "tgag", "acct", "acag", "ctca", "gatc", "ggag", "gtca", "agtg", "taga", "tcta", "cctc", "taaa", "tgac", "ttta", "gaca", "gagg", "actc", "gagt", "actg", "tggc", "cagt", "taac", "gtta", "gtgg", "ggga", "tagt", "ccac", "gctg", "acta", "cacc", "cagc", "tccc", "ctgc", "gcag", "ggta", "gtct", "tacc", "ccag", "cact", "gtag", "ctgg", "ggtg", "agac", "tgcc", "tggg", "agtc", "cagg", "gact", "ggca", "cctg", "aggg", "taag", "ctaa", "gcca", "tagc", "ccct", "ctta", "gctc", "gcta", "ttag", "gagc", "ttaa", "ggct", "ctac", "ttcg", "gcac", "tcga", "tagg", "gtac", "cgaa", "ccta", "gtgc", "gggt", "tcgt", "aggc", "cgtt", "acga", "gcct", "aacg", "atcg", "cgat", "gacc", "ggcc", "ggtc", "acgt", "ggac", "accc", "ccca", "gtcc", "ctag", "ctcg", "agcc", "ccga", "cgag", "acgg", "tcgc", "cgct", "cgga", "tccg", "cgtc", "cacg", "gacg", "tcgg", "agcg", "gcga", "cgcc", "gccc", "accg", "cgta", "gccg", "tacg", "cgtg", "ccgc", "cgac", "ggcg", "gggc", "cggt", "ccgg", "cggc", "gtcg", "ccgt", "cgca", "tgcg", "gcgg", "cggg", "acgc", "gcgt", "cccg"};
    private static final String[] RE5MERS = {"aaatt", "aagtt", "acagt", "acggt", "agact", "aggct", "ataat", "atgat", "caatg", "cagtg", "ccagg", "ccggg", "cgacg", "cggcg", "ctaag", "ctgag", "gaatc", "gagtc", "gcagc", "gcggc", "ggacc", "gggcc", "gtaac", "gtgac", "taata", "tagta", "tcaga", "tcgga", "tgaca", "tggca", "ttaaa", "ttgaa", "aactt", "aattt", "accgt", "actgt", "agcct", "agtct", "atcat", "attat", "cactg", "cattg", "cccgg", "cctgg", "cgccg", "cgtcg", "ctcag", "cttag", "gactc", "gattc", "gccgc", "gctgc", "ggccc", "ggtcc", "gtcac", "gttac", "tacta", "tatta", "tccga", "tctga", "tgcca", "tgtca", "ttcaa", "tttaa"};
    private static final String[] R6MERS = {"aaattt", "tttaaa", "aatatt", "taatta", "attaat", "ttataa", "atgcat", "tcatga", "ttgcaa", "acatgt", "aagctt", "catatg", "tgatca", "gaattc", "tctaga", "agatct", "tgtaca", "caattg", "cttaag", "agtact", "ccatgg", "tagcta", "actagt", "gatatc", "ctatag", "tggcca", "ctgcag", "gtatac", "gcatgc", "gttaac", "cagctg", "aggcct", "gagctc", "ggatcc", "cctagg", "gtgcac", "ggtacc", "cacgtg", "gctagc", "ttcgaa", "aacgtt", "tacgta", "ctcgag", "atcgat", "gggccc", "tccgga", "gacgtc", "agcgct", "tgcgca", "accggt", "ggcgcc", "cccggg", "gtcgac", "gccggc", "gcgcgc", "ccgcgg", "tcgcga", "cggccg", "acgcgt", "cgatcg", "cgtacg", "atttta", "aattaa", "ttaatt", "taaaat", "atccta", "atggta", "agttga", "aaccaa", "aaggaa", "ttcctt", "ttggtt", "tcaact", "tgaagt", "taccat", "taggat", "ctaatc", "cttttc", "caaaac", "cattac", "gtaatg", "gttttg", "gaaaag", "gattag", "acttca", "acccca", "acggca", "agccga", "agggga", "tcccct", "tcggct", "tgccgt", "tggggt", "ctggtc", "ccaacc", "ccttcc", "cgaagc", "cgttgc", "caggac", "gtcctg", "gcaacg", "gcttcg", "ggaagg", "ggttgg", "gaccag", "ccggcc", "cggggc", "gccccg", "ggccgg"};
    private static final String[] R8MERS = {"aaattaaa", "aaaccaaa", "aaaggaaa", "aataataa", "aatggtaa", "aatcctaa", "aacaacaa", "aacttcaa", "aacggcaa", "aagaagaa", "aagttgaa", "aagccgaa", "ataccata", "ataggata", "attaatta", "attggtta", "attcctta", "atcaacta", "atcttcta", "atcggcta", "atgaagta", "atgttgta", "atgccgta", "acattaca", "acaggaca", "actaatca", "actggtca", "actcctca", "accaacca", "accttcca", "accggcca", "acgaagca", "acgttgca", "acgccgca", "agattaga", "agaccaga", "agtaatga", "agtggtga", "agtcctga", "agcaacga", "agcttcga", "agcggcga", "aggaagga", "aggttgga", "aggccgga", "ttattatt", "ttaccatt", "ttaggatt", "tttaattt", "tttggttt", "tttccttt", "ttcaactt", "ttcttctt", "ttcggctt", "ttgaagtt", "ttgttgtt", "ttgccgtt", "taattaat", "taaccaat", "taaggaat", "tatggtat", "tatcctat", "tacaacat", "tacttcat", "tacggcat", "tagaagat", "tagttgat", "tagccgat", "tgattagt", "tgaccagt", "tgaggagt", "tgtaatgt", "tgtcctgt", "tgcaacgt", "tgcttcgt", "tgcggcgt", "tggaaggt", "tggttggt", "tggccggt", "tcattact", "tcaccact", "tcaggact", "tctaatct", "tctggtct", "tccaacct", "tccttcct", "tccggcct", "tcgaagct", "tcgttgct", "tcgccgct", "ccattacc", "ccaccacc", "ccaggacc", "cctaatcc", "cctggtcc", "cctcctcc", "cccaaccc", "cccttccc", "cccggccc", "ccgaagcc", "ccgttgcc", "ccgccgcc", "caattaac", "caaccaac", "caaggaac", "cataatac", "catggtac", "catcctac", "cacttcac", "cacggcac", "cagaagac", "cagttgac", "cagccgac", "ctattatc", "ctaccatc", "ctaggatc", "cttaattc", "cttggttc", "cttccttc", "ctcaactc", "ctcggctc", "ctgaagtc", "ctgttgtc", "ctgccgtc", "cgattagc", "cgaccagc", "cgaggagc", "cgtaatgc", "cgtggtgc", "cgtcctgc", "cgcaacgc", "cgcttcgc", "cggaaggc", "cggttggc", "cggccggc", "ggattagg", "ggaccagg", "ggaggagg", "ggtaatgg", "ggtggtgg", "ggtcctgg", "ggcaacgg", "ggcttcgg", "ggcggcgg", "gggaaggg", "gggttggg", "gggccggg", "gaattaag", "gaaccaag", "gaaggaag", "gataatag", "gatggtag", "gatcctag", "gacaacag", "gacttcag", "gacggcag", "gagttgag", "gagccgag", "gtattatg", "gtaccatg", "gtaggatg", "gttaattg", "gttggttg", "gttccttg", "gtcaactg", "gtcttctg", "gtcggctg", "gtgaagtg", "gtgccgtg", "gcattacg", "gcaccacg", "gcaggacg", "gctaatcg", "gctggtcg", "gctcctcg", "gccaaccg", "gccttccg", "gccggccg", "gcgaagcg", "gcgttgcg"};
    private static final String[] R10MERS = {"aaataataaa", "aaattttaaa", "aaatcctaaa", "aaatggtaaa", "aaacaacaaa", "aaacttcaaa", "aaaccccaaa", "aaacggcaaa", "aaagaagaaa", "aaagttgaaa", "aaagccgaaa", "aaaggggaaa", "aattaattaa", "aattccttaa", "aattggttaa", "aataaaataa", "aatattataa", "aataccataa", "aataggataa", "aatgaagtaa", "aatgttgtaa", "aatgccgtaa", "aatggggtaa", "aatcaactaa", "aatcttctaa", "aatcccctaa", "aatcggctaa", "aaccaaccaa", "aaccttccaa", "aaccggccaa", "aacaaaacaa", "aacattacaa", "aacaccacaa", "aacaggacaa", "aactaatcaa", "aacttttcaa", "aactcctcaa", "aactggtcaa", "aacgaagcaa", "aacgttgcaa", "aacgccgcaa", "aacggggcaa", "aaggaaggaa", "aaggttggaa", "aaggccggaa", "aagaaaagaa", "aagattagaa", "aagaccagaa", "aagaggagaa", "aagtaatgaa", "aagttttgaa", "aagtcctgaa", "aagtggtgaa", "aagcaacgaa", "aagcttcgaa", "aagccccgaa", "aagcggcgaa", "ataattaata", "ataaccaata", "ataaggaata", "atacaacata", "atacttcata", "ataccccata", "atacggcata", "atagaagata", "atagttgata", "atagccgata", "ataggggata", "atttaattta", "atttccttta", "atttggttta", "attaaaatta", "attattatta", "attaccatta", "attaggatta", "attgaagtta", "attgttgtta", "attgccgtta", "attggggtta", "attcaactta", "attcttctta", "attcccctta", "attcggctta", "atccaaccta", "atccttccta", "atccggccta", "atcaaaacta", "atcattacta", "atcaccacta", "atcaggacta", "atctaatcta", "atcttttcta", "atctcctcta", "atctggtcta", "atcgaagcta", "atcgttgcta", "atcgccgcta", "atcggggcta", "atggaaggta", "atggttggta", "atggccggta", "atgaaaagta", "atgattagta", "atgaccagta", "atgaggagta", "atgtaatgta", "atgttttgta", "atgtcctgta", "atgtggtgta", "atgcaacgta", "atgcttcgta", "atgccccgta", "atgcggcgta", "acaattaaca", "acaaccaaca", "acaaggaaca", "acataataca", "acattttaca", "acatcctaca", "acatggtaca", "acagaagaca", "acagttgaca", "acagccgaca", "acaggggaca", "acttaattca", "acttccttca", "acttggttca", "actaaaatca", "actattatca", "actaccatca", "actaggatca", "actgaagtca", "actgttgtca", "actgccgtca", "actggggtca", "actcaactca", "actcttctca", "actcccctca", "actcggctca", "acccaaccca", "acccttccca", "acccggccca", "accaaaacca", "accattacca", "accaccacca", "accaggacca", "acctaatcca", "accttttcca", "acctcctcca", "acctggtcca", "accgaagcca", "accgttgcca", "accgccgcca", "accggggcca", "acggaaggca", "acggttggca", "acggccggca", "acgaaaagca", "acgattagca", "acgaccagca", "acgaggagca", "acgtaatgca", "acgttttgca", "acgtcctgca", "acgtggtgca", "acgcaacgca", "acgcttcgca", "acgccccgca", "acgcggcgca", "agaattaaga", "agaaccaaga", "agaaggaaga", "agataataga", "agattttaga", "agatcctaga", "agatggtaga", "agacaacaga", "agacttcaga", "agaccccaga", "agacggcaga", "agttaattga", "agttccttga", "agttggttga", "agtaaaatga", "agtattatga", "agtaccatga", "agtaggatga", "agtgaagtga", "agtgttgtga", "agtgccgtga", "agtggggtga", "agtcaactga", "agtcttctga", "agtcccctga", "agtcggctga", "agccaaccga", "agccttccga", "agccggccga", "agcaaaacga", "agcattacga", "agcaccacga", "agcaggacga", "agctaatcga", "agcttttcga", "agctcctcga", "agctggtcga", "agcgaagcga", "agcgttgcga", "agcgccgcga", "agcggggcga", "agggaaggga", "agggttggga", "agggccggga", "aggaaaagga", "aggattagga", "aggaccagga", "aggaggagga", "aggtaatgga", "aggttttgga", "aggtcctgga", "aggtggtgga", "aggcaacgga", "aggcttcgga", "aggccccgga", "aggcggcgga", "ttaattaatt", "ttaaccaatt", "ttaaggaatt", "ttataatatt", "ttattttatt", "ttatcctatt", "ttatggtatt", "ttacaacatt", "ttacttcatt", "ttaccccatt", "ttacggcatt", "ttagaagatt", "ttagttgatt", "ttagccgatt", "ttaggggatt", "tttaaaattt", "tttattattt", "tttaccattt", "tttaggattt", "tttgaagttt", "tttgttgttt", "tttgccgttt", "tttggggttt", "tttcaacttt", "tttcttcttt", "tttccccttt", "tttcggcttt", "ttccaacctt", "ttccttcctt", "ttccggcctt", "ttcaaaactt", "ttcattactt", "ttcaccactt", "ttcaggactt", "ttctaatctt", "ttcttttctt", "ttctcctctt", "ttctggtctt", "ttcgaagctt", "ttcgttgctt", "ttcgccgctt", "ttcggggctt", "ttggaaggtt", "ttggttggtt", "ttggccggtt", "ttgaaaagtt", "ttgattagtt", "ttgaccagtt", "ttgaggagtt", "ttgtaatgtt", "ttgttttgtt", "ttgtcctgtt", "ttgtggtgtt", "ttgcaacgtt", "ttgcttcgtt", "ttgccccgtt", "ttgcggcgtt", "taaattaaat", "taaaccaaat", "taaaggaaat", "taataataat", "taattttaat", "taatcctaat", "taatggtaat", "taacaacaat", "taacttcaat", "taaccccaat", "taacggcaat", "taagaagaat", "taagttgaat", "taagccgaat", "taaggggaat", "tattaattat", "tattccttat", "tattggttat", "tatgaagtat", "tatgttgtat", "tatgccgtat", "tatggggtat", "tatcaactat", "tatcttctat", "tatcccctat", "tatcggctat", "taccaaccat", "taccttccat", "taccggccat", "tacaaaacat", "tacattacat", "tacaccacat", "tacaggacat", "tactaatcat", "tacttttcat", "tactcctcat", "tactggtcat", "tacgaagcat", "tacgttgcat", "tacgccgcat", "tacggggcat", "taggaaggat", "taggttggat", "taggccggat", "tagaaaagat", "tagattagat", "tagaccagat", "tagaggagat", "tagtaatgat", "tagttttgat", "tagtcctgat", "tagtggtgat", "tagcaacgat", "tagcttcgat", "tagccccgat", "tagcggcgat", "tgaattaagt", "tgaaccaagt", "tgaaggaagt", "tgataatagt", "tgattttagt", "tgatcctagt", "tgatggtagt", "tgacaacagt", "tgacttcagt", "tgaccccagt", "tgacggcagt", "tgagaagagt", "tgagttgagt", "tgagccgagt", "tgaggggagt", "tgttaattgt", "tgttccttgt", "tgttggttgt", "tgtaaaatgt", "tgtattatgt", "tgtaccatgt", "tgtaggatgt", "tgtcaactgt", "tgtcttctgt", "tgtcccctgt", "tgtcggctgt", "tgccaaccgt", "tgccttccgt", "tgccggccgt", "tgcaaaacgt", "tgcattacgt", "tgcaccacgt", "tgcaggacgt", "tgctaatcgt", "tgcttttcgt", "tgctcctcgt", "tgctggtcgt", "tgcgaagcgt", "tgcgttgcgt", "tgcgccgcgt", "tgcggggcgt", "tgggaagggt", "tgggttgggt", "tgggccgggt", "tggaaaaggt", "tggattaggt", "tggaccaggt", "tggaggaggt", "tggtaatggt", "tggttttggt", "tggtcctggt", "tggtggtggt", "tggcaacggt", "tggcttcggt", "tggccccggt", "tggcggcggt", "tcaattaact", "tcaaccaact", "tcaaggaact", "tcataatact", "tcattttact", "tcatcctact", "tcatggtact", "tcacaacact", "tcacttcact", "tcaccccact", "tcacggcact", "tcagaagact", "tcagttgact", "tcagccgact", "tcaggggact", "tcttaattct", "tcttccttct", "tcttggttct", "tctaaaatct", "tctattatct", "tctaccatct", "tctaggatct", "tctgaagtct", "tctgttgtct", "tctgccgtct", "tctggggtct", "tcccaaccct", "tcccttccct", "tcccggccct", "tccaaaacct", "tccattacct", "tccaccacct", "tccaggacct", "tcctaatcct", "tccttttcct", "tcctcctcct", "tcctggtcct", "tccgaagcct", "tccgttgcct", "tccgccgcct", "tccggggcct", "tcggaaggct", "tcggttggct", "tcggccggct", "tcgaaaagct", "tcgattagct", "tcgaccagct", "tcgaggagct", "tcgtaatgct", "tcgttttgct", "tcgtcctgct", "tcgtggtgct", "tcgcaacgct", "tcgcttcgct", "tcgccccgct", "tcgcggcgct", "ccaattaacc", "ccaaccaacc", "ccaaggaacc", "ccataatacc", "ccattttacc", "ccatcctacc", "ccatggtacc", "ccacaacacc", "ccacttcacc", "ccaccccacc", "ccacggcacc", "ccagaagacc", "ccagttgacc", "ccagccgacc", "ccaggggacc", "ccttaattcc", "ccttccttcc", "ccttggttcc", "cctaaaatcc", "cctattatcc", "cctaccatcc", "cctaggatcc", "cctgaagtcc", "cctgttgtcc", "cctgccgtcc", "cctggggtcc", "cctcaactcc", "cctcttctcc", "cctcccctcc", "cctcggctcc", "cccaaaaccc", "cccattaccc", "cccaccaccc", "cccaggaccc", "ccctaatccc", "cccttttccc", "ccctcctccc", "ccctggtccc", "cccgaagccc", "cccgttgccc", "cccgccgccc", "cccggggccc", "ccggaaggcc", "ccggttggcc", "ccggccggcc", "ccgaaaagcc", "ccgattagcc", "ccgaccagcc", "ccgaggagcc", "ccgtaatgcc", "ccgttttgcc", "ccgtcctgcc", "ccgtggtgcc", "ccgcaacgcc", "ccgcttcgcc", "ccgccccgcc", "ccgcggcgcc", "caaattaaac", "caaaccaaac", "caaaggaaac", "caataataac", "caattttaac", "caatcctaac", "caatggtaac", "caacaacaac", "caacttcaac", "caaccccaac", "caacggcaac", "caagaagaac", "caagttgaac", "caagccgaac", "caaggggaac", "cattaattac", "cattccttac", "cattggttac", "cataaaatac", "catattatac", "cataccatac", "cataggatac", "catgaagtac", "catgttgtac", "catgccgtac", "catggggtac", "catcaactac", "catcttctac", "catcccctac", "catcggctac", "caccaaccac", "caccttccac", "caccggccac", "cactaatcac", "cacttttcac", "cactcctcac", "cactggtcac", "cacgaagcac", "cacgttgcac", "cacgccgcac", "cacggggcac", "caggaaggac", "caggttggac", "caggccggac", "cagaaaagac", "cagattagac", "cagaccagac", "cagaggagac", "cagtaatgac", "cagttttgac", "cagtcctgac", "cagtggtgac", "cagcaacgac", "cagcttcgac", "cagccccgac", "cagcggcgac", "ctaattaatc", "ctaaccaatc", "ctaaggaatc", "ctataatatc", "ctattttatc", "ctatcctatc", "ctatggtatc", "ctacaacatc", "ctacttcatc", "ctaccccatc", "ctacggcatc", "ctagaagatc", "ctagttgatc", "ctagccgatc", "ctaggggatc", "ctttaatttc", "ctttcctttc", "ctttggtttc", "cttaaaattc", "cttattattc", "cttaccattc", "cttaggattc", "cttgaagttc", "cttgttgttc", "cttgccgttc", "cttggggttc", "cttcaacttc", "cttcttcttc", "cttccccttc", "cttcggcttc", "ctccaacctc", "ctccttcctc", "ctccggcctc", "ctcaaaactc", "ctcattactc", "ctcaccactc", "ctcaggactc", "ctcgaagctc", "ctcgttgctc", "ctcgccgctc", "ctcggggctc", "ctggaaggtc", "ctggttggtc", "ctggccggtc", "ctgaaaagtc", "ctgattagtc", "ctgaccagtc", "ctgaggagtc", "ctgtaatgtc", "ctgttttgtc", "ctgtcctgtc", "ctgtggtgtc", "ctgcaacgtc", "ctgcttcgtc", "ctgccccgtc", "ctgcggcgtc", "cgaattaagc", "cgaaccaagc", "cgaaggaagc", "cgataatagc", "cgattttagc", "cgatcctagc", "cgatggtagc", "cgacaacagc", "cgacttcagc", "cgaccccagc", "cgacggcagc", "cgagaagagc", "cgagttgagc", "cgagccgagc", "cgaggggagc", "cgttaattgc", "cgttccttgc", "cgttggttgc", "cgtaaaatgc", "cgtattatgc", "cgtaccatgc", "cgtaggatgc", "cgtgaagtgc", "cgtgttgtgc", "cgtgccgtgc", "cgtggggtgc", "cgtcaactgc", "cgtcttctgc", "cgtcccctgc", "cgtcggctgc", "cgccaaccgc", "cgccttccgc", "cgccggccgc", "cgcaaaacgc", "cgcattacgc", "cgcaccacgc", "cgcaggacgc", "cgctaatcgc", "cgcttttcgc", "cgctcctcgc", "cgctggtcgc", "cgggaagggc", "cgggttgggc", "cgggccgggc", "cggaaaaggc", "cggattaggc", "cggaccaggc", "cggaggaggc", "cggtaatggc", "cggttttggc", "cggtcctggc", "cggtggtggc", "cggcaacggc", "cggcttcggc", "cggccccggc", "cggcggcggc", "ggaattaagg", "ggaaccaagg", "ggaaggaagg", "ggataatagg", "ggattttagg", "ggatcctagg", "ggatggtagg", "ggacaacagg", "ggacttcagg", "ggaccccagg", "ggacggcagg", "ggagaagagg", "ggagttgagg", "ggagccgagg", "ggaggggagg", "ggttaattgg", "ggttccttgg", "ggttggttgg", "ggtaaaatgg", "ggtattatgg", "ggtaccatgg", "ggtaggatgg", "ggtgaagtgg", "ggtgttgtgg", "ggtgccgtgg", "ggtggggtgg", "ggtcaactgg", "ggtcttctgg", "ggtcccctgg", "ggtcggctgg", "ggccaaccgg", "ggccttccgg", "ggccggccgg", "ggcaaaacgg", "ggcattacgg", "ggcaccacgg", "ggcaggacgg", "ggctaatcgg", "ggcttttcgg", "ggctcctcgg", "ggctggtcgg", "ggcgaagcgg", "ggcgttgcgg", "ggcgccgcgg", "ggcggggcgg", "gggaaaaggg", "gggattaggg", "gggaccaggg", "gggaggaggg", "gggtaatggg", "gggttttggg", "gggtcctggg", "gggtggtggg", "gggcaacggg", "gggcttcggg", "gggccccggg", "gggcggcggg", "gaaattaaag", "gaaaccaaag", "gaaaggaaag", "gaataataag", "gaattttaag", "gaatcctaag", "gaatggtaag", "gaacaacaag", "gaacttcaag", "gaaccccaag", "gaacggcaag", "gaagaagaag", "gaagttgaag", "gaagccgaag", "gaaggggaag", "gattaattag", "gattccttag", "gattggttag", "gataaaatag", "gatattatag", "gataccatag", "gataggatag", "gatgaagtag", "gatgttgtag", "gatgccgtag", "gatggggtag", "gatcaactag", "gatcttctag", "gatcccctag", "gatcggctag", "gaccaaccag", "gaccttccag", "gaccggccag", "gacaaaacag", "gacattacag", "gacaccacag", "gacaggacag", "gactaatcag", "gacttttcag", "gactcctcag", "gactggtcag", "gacgaagcag", "gacgttgcag", "gacgccgcag", "gacggggcag", "gaggaaggag", "gaggttggag", "gaggccggag", "gagtaatgag", "gagttttgag", "gagtcctgag", "gagtggtgag", "gagcaacgag", "gagcttcgag", "gagccccgag", "gagcggcgag", "gtaattaatg", "gtaaccaatg", "gtaaggaatg", "gtataatatg", "gtattttatg", "gtatcctatg", "gtatggtatg", "gtacaacatg", "gtacttcatg", "gtaccccatg", "gtacggcatg", "gtagaagatg", "gtagttgatg", "gtagccgatg", "gtaggggatg", "gtttaatttg", "gtttcctttg", "gtttggtttg", "gttaaaattg", "gttattattg", "gttaccattg", "gttaggattg", "gttgaagttg", "gttgttgttg", "gttgccgttg", "gttggggttg", "gttcaacttg", "gttcttcttg", "gttccccttg", "gttcggcttg", "gtccaacctg", "gtccttcctg", "gtccggcctg", "gtcaaaactg", "gtcattactg", "gtcaccactg", "gtcaggactg", "gtctaatctg", "gtcttttctg", "gtctcctctg", "gtctggtctg", "gtcgaagctg", "gtcgttgctg", "gtcgccgctg", "gtcggggctg", "gtggaaggtg", "gtggttggtg", "gtggccggtg", "gtgaaaagtg", "gtgattagtg", "gtgaccagtg", "gtgaggagtg", "gtgcaacgtg", "gtgcttcgtg", "gtgccccgtg", "gtgcggcgtg", "gcaattaacg", "gcaaccaacg", "gcaaggaacg", "gcataatacg", "gcattttacg", "gcatcctacg", "gcatggtacg", "gcacaacacg", "gcacttcacg", "gcaccccacg", "gcacggcacg", "gcagaagacg", "gcagttgacg", "gcagccgacg", "gcaggggacg", "gcttaattcg", "gcttccttcg", "gcttggttcg", "gctaaaatcg", "gctattatcg", "gctaccatcg", "gctaggatcg", "gctgaagtcg", "gctgttgtcg", "gctgccgtcg", "gctggggtcg", "gctcaactcg", "gctcttctcg", "gctcccctcg", "gctcggctcg", "gcccaacccg", "gcccttcccg", "gcccggcccg", "gccaaaaccg", "gccattaccg", "gccaccaccg", "gccaggaccg", "gcctaatccg", "gccttttccg", "gcctcctccg", "gcctggtccg", "gccgaagccg", "gccgttgccg", "gccgccgccg", "gccggggccg", "gcggaaggcg", "gcggttggcg", "gcggccggcg", "gcgaaaagcg", "gcgattagcg", "gcgaccagcg", "gcgaggagcg", "gcgtaatgcg", "gcgttttgcg", "gcgtcctgcg", "gcgtggtgcg"};
}
