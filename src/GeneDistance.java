import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class GeneDistance {

    private static final String TITLE = "GeneDistance (2024-2026) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/GeneDistance";

    private static final int OK = 0;
    private static final int USAGE_ERROR = 1;
    private static final int IO_ERROR = 2;

    public static void main(String[] args) {
        if (args.length == 0) {
            printHelp();
            return;
        }

        String infile = null;
        int kmer = 4;
        int KmerCounter = 0;
        int statistic = 0;   // 0 = the k-mer spacing measure, otherwise COSINE / D2STAR / VECTOR
        boolean spacing = false;  // -spacing asks for the k-mer spacing measure (the former default)
        int contain = 0;     // 1 for a containment run
        int scan = 0;        // the sequence to look for, 1-based, in a scan
        int ksize = AUTO_K;  // the k of the exact k-mers those two modes use
        int scaled = 1;      // FracMinHash downsampling for -contain/-scan (1 = exact); AUTO_SCALE = choose it

        for (String arg : args) {
            String a = arg.trim();
            if (a.isEmpty()) {
                continue;
            }
            String low = a.toLowerCase(Locale.ROOT);

            // whatever is not an option is the target file or folder
            if (!(a.startsWith("-") || a.equals("/?") || low.startsWith("kmer="))) {
                if (infile != null) {
                    fail("Unexpected argument: " + a + " (the target is already " + infile + ")");
                }
                infile = a;
                continue;
            }

            String opt = low.equals("/?") ? "help" : low;
            while (opt.startsWith("-")) {
                opt = opt.substring(1);
            }

            if (opt.equals("help") || opt.equals("h") || opt.equals("?")) {
                printHelp();
                return;
            } else if (opt.equals("version")) {
                System.out.println(TITLE);
                return;
            } else if (opt.equals("kmerstat")) {
                KmerCounter = 1;
            } else if (opt.equals("kmer2stat")) {
                KmerCounter = 2;
            } else if (opt.equals("cosine")) {
                statistic = SequencesSimilarity.COSINE;
            } else if (opt.equals("d2star") || opt.equals("d2*")) {
                statistic = SequencesSimilarity.D2STAR;
            } else if (opt.equals("vector")) {
                statistic = SequencesSimilarity.VECTOR;
            } else if (opt.equals("spacing") || opt.equals("space")) {
                spacing = true;
            } else if (opt.equals("contain")) {
                contain = 1;
            } else if (opt.startsWith("contain=")) {
                contain = 1;
                ksize = exactK(opt.substring(8).trim(), a);
            } else if (opt.equals("scan")) {
                scan = 1;
            } else if (opt.startsWith("scan=")) {
                scan = queryIndex(opt.substring(5).trim(), a);
            } else if (opt.startsWith("ksize=")) {
                ksize = exactK(opt.substring(6).trim(), a);
            } else if (opt.equals("scaled")) {
                scaled = AUTO_SCALE;
            } else if (opt.startsWith("scaled=")) {
                scaled = scaledOf(opt.substring(7).trim(), a);
            } else if (opt.startsWith("kmer=")) {
                kmer = model(opt.substring(5).trim(), a);
            } else {
                //fail("Unknown option: " + a);
            }
        }

        if (infile == null) {
            fail("No target file or folder given.");
        }

        // -vector is the default measure: with no measure and no other mode named,
        // run it. -spacing asks for the former default (the k-mer spacing measure).
        if (statistic == 0 && !spacing && KmerCounter == 0 && contain == 0 && scan == 0) {
            statistic = SequencesSimilarity.VECTOR;
        }

        File target = new File(infile);
        if (!target.exists()) {
            fail("No such file or folder: " + infile);
        }

        String tag = switch (statistic) {
            case SequencesSimilarity.COSINE -> "_cos";
            case SequencesSimilarity.D2STAR -> "_d2s";
            case SequencesSimilarity.VECTOR -> "_vec";
            default -> "";
        };

        boolean exact = contain != 0 || scan != 0;

        System.out.println("Current Directory: " + System.getProperty("user.dir"));
        System.out.println("Target file or Folder: " + infile);
        if (!exact) {
            System.out.println("K-mer set: " + kmer + " (" + SequencesSimilarity.modelInfo(kmer) + ")");
        }
        System.out.println("Measure: " + (scan != 0 ? "scan (where each sequence holds the query)"
                : contain != 0 ? "containment (how much of each sequence occurs in each other one)"
                : switch (statistic) {
                    case SequencesSimilarity.COSINE -> "cosine of the k-mer frequencies";
                    case SequencesSimilarity.D2STAR -> "d2* (k-mer frequencies, centred on the base composition)";
                    case SequencesSimilarity.VECTOR -> "vector (d2*, scale-matched windows, both strands) [default]";
                    default -> "k-mer spacing ratios";
                }));

        String baseName = "";
        String outBase;
        List<String> files = new ArrayList<>();

        if (target.isDirectory()) {
            baseName = target.getName();
            outBase = target.toPath().toString() + File.separator + "result";

            File[] found = target.listFiles();
            if (found == null) {
                fail("Cannot read the folder: " + infile);
                return;
            }
            // the reports are all written into this same folder as result*.xls / result*.meg,
            // whatever the measure and its k (which -contain and -scan only settle later): never
            // read one of them back in as an input sequence
            for (File file : found) {
                if (file.isFile()) {
                    String name = file.getName();
                    boolean report = name.startsWith("result")
                            && (name.endsWith(".xls") || name.endsWith(".meg"));
                    if (!report) {
                        files.add(file.getAbsolutePath());
                    }
                }
            }
            if (files.isEmpty()) {
                System.err.println("No files: " + target.toString());
                System.exit(USAGE_ERROR);
            }
        } else {
            outBase = target.toPath().toString();
            files.add(infile);
        }

        System.exit(SaveResult(KmerCounter, statistic, contain, scan, ksize, scaled, tag, kmer, baseName, files.toArray(new String[0]), outBase));
    }

    private static int SaveResult(int KmerCounter, int statistic, int contain, int scan, int ksize, int scaled,
            String tag, int model, String baseName, String[] infiles, String folder) {
        try {
            long startTime = System.nanoTime();
            System.out.println("Running...");
            SequencesSimilarity s1 = new SequencesSimilarity(model);
            s1.SetTag(tag);

            // -contain streams the files and keeps only the sketches, never the
            // sequences, so its memory is bounded whatever the genome size.
            if (contain != 0) {
                long[] lens = StreamingFasta.recordLengths(infiles);
                int nseq = lens.length;
                if (nseq == 0) {
                    System.out.println("No sequences found to process.");
                    return USAGE_ERROR;
                }
                long longest = 0;
                for (long L : lens) {
                    if (L > longest) {
                        longest = L;
                    }
                }
                int k = (ksize == AUTO_K) ? SequencesSimilarity.containmentK(longest) : ksize;
                int sc = (scaled == AUTO_SCALE) ? autoScaled(longest) : scaled;
                System.out.println("Sequences: " + nseq + ", k=" + k
                        + (ksize == AUTO_K ? " (chosen for the longest sequence, " + longest + " bp)" : "")
                        + (sc > 1 ? ", FracMinHash scaled=" + sc + (scaled == AUTO_SCALE ? " (chosen for size)" : "") : ""));
                s1.SetPaths(folder, "_c" + k + (sc > 1 ? "s" + sc : ""));
                s1.RunContainmentStream(infiles, k, sc);
                System.out.println("Report: " + s1.getReportFile());
                System.out.println("MEGA:   " + s1.getMegaFile());
                long t = (System.nanoTime() - startTime) / 1000000;
                System.out.println("Time taken: " + (t / 1000) + "." + String.format("%03d", t % 1000) + " seconds");
                return OK;
            }

            // -scan streams the subjects one at a time, so the corpus is never all
            // in memory: only the query index and the subject being scanned.
            if (scan != 0) {
                long[] lens = StreamingFasta.recordLengths(infiles);
                int nseq = lens.length;
                if (nseq == 0) {
                    System.out.println("No sequences found to process.");
                    return USAGE_ERROR;
                }
                if (scan > nseq) {
                    System.err.println("There is no sequence " + scan + ": the input holds " + nseq + ".");
                    return USAGE_ERROR;
                }
                if (nseq < 2) {
                    System.err.println("A scan needs the query and at least one sequence to look through.");
                    return USAGE_ERROR;
                }
                long longest = 0;
                for (long L : lens) {
                    if (L > longest) {
                        longest = L;
                    }
                }
                int k = (ksize == AUTO_K) ? SequencesSimilarity.containmentK(longest) : ksize;
                int sc = (scaled == AUTO_SCALE) ? autoScaled(longest) : scaled;
                System.out.println("Sequences: " + nseq + ", k=" + k
                        + (ksize == AUTO_K ? " (chosen for the longest sequence, " + longest + " bp)" : "")
                        + (sc > 1 ? ", FracMinHash scaled=" + sc + (scaled == AUTO_SCALE ? " (chosen for size)" : "") : ""));
                s1.SetPaths(folder, "_scan" + k + (sc > 1 ? "s" + sc : ""));
                s1.RunScanStream(infiles, scan - 1, k, sc);
                System.out.println("Report: " + s1.getReportFile());
                long t = (System.nanoTime() - startTime) / 1000000;
                System.out.println("Time taken: " + (t / 1000) + "." + String.format("%03d", t % 1000) + " seconds");
                return OK;
            }

            // -d2star and -cosine are purely sequential: they stream each record into
            // a small vector and never hold the sequences, so their memory is bounded.
            if (statistic == SequencesSimilarity.COSINE || statistic == SequencesSimilarity.D2STAR) {
                long[] lens = StreamingFasta.recordLengths(infiles);
                int nseq = lens.length;
                if (nseq == 0) {
                    System.out.println("No sequences found to process.");
                    return USAGE_ERROR;
                }
                System.out.println("Sequences: " + nseq + ", k-mers: " + s1.getKmerCount());
                s1.SetPaths(folder, "_k" + s1.getModel() + tag);
                s1.RunFrequencyStream(infiles, statistic);
                System.out.println("Report: " + s1.getReportFile());
                System.out.println("MEGA:   " + s1.getMegaFile());
                long t = (System.nanoTime() - startTime) / 1000000;
                System.out.println("Time taken: " + (t / 1000) + "." + String.format("%03d", t % 1000) + " seconds");
                return OK;
            }

            // -vector needs random access (windows of the long sequence, random null
            // draws), so it streams the input into a segmented store rather than an
            // array of Strings: the corpus and a single record may exceed 2 GB.
            if (statistic == SequencesSimilarity.VECTOR) {
                SeqStore st = SeqStore.build(infiles);
                int nseq = st.count();
                if (nseq == 0) {
                    System.out.println("No sequences found to process.");
                    return USAGE_ERROR;
                }
                System.out.println("Sequences: " + nseq + ", k-mers: " + s1.getKmerCount());
                s1.SetPaths(folder, "_k" + s1.getModel() + tag);
                s1.RunVector(st);
                System.out.println("Report: " + s1.getReportFile());
                System.out.println("MEGA:   " + s1.getMegaFile());
                long t = (System.nanoTime() - startTime) / 1000000;
                System.out.println("Time taken: " + (t / 1000) + "." + String.format("%03d", t % 1000) + " seconds");
                return OK;
            }

            // the k-mer spacing measure (-spacing): sequential, so it streams the
            // records into their small count vectors without holding the sequences.
            if (KmerCounter == 0) {
                long[] lens = StreamingFasta.recordLengths(infiles);
                int nseq = lens.length;
                if (nseq == 0) {
                    System.out.println("No sequences found to process.");
                    return USAGE_ERROR;
                }
                System.out.println("Sequences: " + nseq + ", k-mers: " + s1.getKmerCount());
                s1.SetPaths(folder, "_k" + s1.getModel() + tag);
                s1.RunStream(infiles);
                System.out.println("Report: " + s1.getReportFile());
                System.out.println("MEGA:   " + s1.getMegaFile());
                long t = (System.nanoTime() - startTime) / 1000000;
                System.out.println("Time taken: " + (t / 1000) + "." + String.format("%03d", t % 1000) + " seconds");
                return OK;
            }

            s1.SetFolder(baseName, infiles, folder);

            int nseq = s1.getSequenceCount();
            if (nseq == 0) {
                System.out.println("No sequences found to process.");
                return USAGE_ERROR;
            }

            System.out.println("Sequences: " + nseq + ", k-mers: " + s1.getKmerCount());

            s1.RunKmerCounter(KmerCounter);
            System.out.println("Report: " + s1.getReportFile());
            long ms = (System.nanoTime() - startTime) / 1000000;
            System.out.println("Time taken: " + (ms / 1000) + "." + String.format("%03d", ms % 1000) + " seconds");
            return OK;
        } catch (IOException e) {
            System.err.println("File I/O error: " + e.getMessage());
            return IO_ERROR;
        } catch (OutOfMemoryError e) {
            System.err.println("Out of memory: give the JVM more heap, e.g. java -Xmx8g -jar GeneDistance.jar ...");
            return IO_ERROR;
        }
    }

    private static final int AUTO_K = -1;
    private static final int AUTO_SCALE = -1;

    /** Downsampling that keeps the biggest sketch near a couple million k-mers. */
    private static int autoScaled(long longest) {
        long sc = Math.round(longest / 2_000_000.0);
        return (int) Math.max(1, Math.min(100_000, sc));
    }

    private static int scaledOf(String value, String arg) {
        if (value.equalsIgnoreCase("auto") || value.isEmpty()) {
            return AUTO_SCALE;
        }
        int s;
        try {
            s = Integer.parseInt(value);
        } catch (NumberFormatException e) {
            s = 0;
        }
        if (s < 1 || s > 1_000_000) {
            fail("The scaled of " + arg + " must be between 1 and 1000000 "
                    + "(1 = exact; ~1000 fits a 2 Gb genome). Or -scaled=auto to choose it for the size.");
        }
        return s;
    }

    private static int exactK(String value, String arg) {
        int k;
        try {
            k = Integer.parseInt(value);
        } catch (NumberFormatException e) {
            k = 0;
        }
        if (k < 6 || k > 31) {
            fail("The k of " + arg + " must be between 6 and 31 (11-16 suits most cases). "
                    + "Leave it out to have it chosen for the length of the sequences.");
        }
        return k;
    }

    private static int queryIndex(String value, String arg) {
        int q;
        try {
            q = Integer.parseInt(value);
        } catch (NumberFormatException e) {
            q = 0;
        }
        if (q < 1) {
            fail("The query of " + arg + " is the number of a sequence of the input, from 1. "
                    + "Leave it out, -scan, to look for the first one.");
        }
        return q;
    }

    private static int model(String value, String arg) {
        int m;
        try {
            m = Integer.parseInt(value);
        } catch (NumberFormatException e) {
            m = -1;
        }
        if (!SequencesSimilarity.isModel(m)) {
            StringBuilder sb = new StringBuilder();
            for (int v : SequencesSimilarity.MODELS) {
                sb.append(sb.length() == 0 ? "" : "/").append(v);
            }
            fail("Unsupported k-mer set in " + arg + ". Use -kmer=" + sb + ".");
        }
        return m;
    }

    private static void fail(String message) {
        System.err.println(message);
        System.err.println("Run with --help for the usage.");
        System.exit(USAGE_ERROR);
    }

    private static void printHelp() {
        System.out.println(TITLE + "\n");
        System.out.println("Identification of the genetic similarity, or the distance, between genomic sequences.");
        System.out.println();
        System.out.println("Usage:");
        System.out.println("  java -jar GeneDistance.jar <inputfile>|<inputfolder> [options]");
        System.out.println();
        System.out.println("Input:");
        System.out.println("  A FASTA file ('>' header line, then the sequence). Given a folder, every file in it");
        System.out.println("  is read and all the sequences are pooled into one analysis.");
        System.out.println();
        System.out.println("Options:");
        System.out.println("  -kmer=N      set of k-mers used in the analysis (default: -kmer=4):");
        for (int m : SequencesSimilarity.MODELS) {
            System.out.println(String.format("                 %-4d %s", m, SequencesSimilarity.modelInfo(m)));
        }
        System.out.println("  -contain     how much of each sequence occurs in each other one: the containment");
        System.out.println("               index on the whole k-mer space, |kmers(A) & kmers(B)| / |kmers(A)|.");
        System.out.println("               Asymmetric, and the only measure here that finds a short element");
        System.out.println("               inside a long sequence, where a measure of the composition as a whole");
        System.out.println("               is drowned by the rest of it. The report also gives the identity of");
        System.out.println("               the shared part;");
        System.out.println("  -scan        where each sequence holds the query: a window the size of the query");
        System.out.println("  -scan=N      is walked along them, and every region that holds it is reported with");
        System.out.println("               its position, its identity and an E-value. The query is the first");
        System.out.println("               sequence of the input, or the N-th;");
        System.out.println("  -ksize=K     the k that -contain and -scan use (6-31, 11-16 suits most cases).");
        System.out.println("               Left out, it is chosen for the length of the sequences. Note that");
        System.out.println("               -kmer= names one of the sets above and does not apply to those two;");
        System.out.println("  -scaled=N    FracMinHash for -contain / -scan: keep only ~1/N of the k-mers, so");
        System.out.println("               a 2 Gb genome's set shrinks by N (from ~16 GB to ~16 MB at N=1000)");
        System.out.println("               while the containment stays an unbiased estimate. N=1 is exact;");
        System.out.println("               -scaled (or -scaled=auto) chooses N for the size of the input;");
        System.out.println("  -d2star      compare the k-mer frequencies rather than their spacing, each count");
        System.out.println("               centred on the base composition of its own sequence (the d2* measure).");
        System.out.println("               Use it whenever the sequences differ widely in length: the default");
        System.out.println("               measure only counts a k-mer seen at least twice, so a short sequence");
        System.out.println("               arrives with an almost empty vector and its homology goes unseen;");
        System.out.println("  -cosine      the same, but on the raw frequencies, without centring. Simpler, and");
        System.out.println("               with a high floor: unrelated sequences still score around 70%;");
        System.out.println("  -vector      the DEFAULT measure (used when no other is named): d2* with");
        System.out.println("               scale-matched windows and both strands, for fragments of very different");
        System.out.println("               length - a short one is compared to the best window of a longer one, on");
        System.out.println("               either strand, and scores no better than chance are dropped (the chance");
        System.out.println("               level is measured from the input per length scale);");
        System.out.println("  -spacing     the former default: how the k-mers are spaced inside each sequence;");
        System.out.println("  -kmerstat    report the k-mer distances of each sequence, instead of comparing them;");
        System.out.println("  -kmer2stat   report the k-mer distances averaged over all the sequences;");
        System.out.println("  -h, --help   print this help and exit;");
        System.out.println("  --version    print the version and exit.");
        System.out.println();
        System.out.println("Output (written next to the input, 'result_kN.*' when a folder is given):");
        System.out.println("  <input>_kN.xls   k-mer counts, k-mer distances and the similarity matrix (%);");
        System.out.println("  <input>_kN.meg   distance matrix for MEGA (lower-left), to build a tree.");
        System.out.println("  -cosine and -d2star write to <input>_kN_cos.* and <input>_kN_d2s.*, so that the");
        System.out.println("  reports of the three measures never overwrite one another.");
        System.out.println();
        System.out.println("Examples:");
        System.out.println("  java -jar GeneDistance.jar sequences.fasta");
        System.out.println("  java -jar GeneDistance.jar sequences.fasta -kmer=6");
        System.out.println("  java -jar GeneDistance.jar E:\\Genomes\\Oryza_sativa -kmer=41");
        System.out.println("  java -jar GeneDistance.jar sequences.fasta -kmerstat");
        System.out.println();
        System.out.println("Exit codes: 0 = done, 1 = wrong usage or no sequence found, 2 = I/O error.");
    }
}
