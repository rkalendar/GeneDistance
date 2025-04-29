import java.io.IOException;
import java.io.File;

public class GeneDistance {

    public static void main(String[] args) {
        if (args.length > 0) {
            String infile = args[0]; // file path or Folder

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            String s = String.join(" ", args).toLowerCase() + " ";
            int model = 0;
            if (s.contains("mod=")) {
                int j = s.indexOf("mod=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    model = StrToInt(s.substring(j + 4, x));
                    if (model < 0) {
                        model = 0;
                    }
                    if (model > 4) {
                        model = 4;
                    }
                }
            }

            File folder = new File(infile);
            if (folder.exists() && (folder.isDirectory() || folder.isFile())) {
                if (folder.isDirectory()) {
                    File[] files = folder.listFiles();

                    if (files.length == 0) {
                        System.err.println("No files: " + folder.toString());
                        return;
                    }

                    int k = -1;
                    String[] filelist = new String[files.length];
                    for (File file : files) {
                        if (file.isFile()) {
                            filelist[++k] = file.getAbsolutePath();
                        }
                    }
                    SaveResult(model, filelist, folder.toPath().toString() + File.separator + "result");
                } else {
                    String[] filelist = new String[1];
                    filelist[0] = infile;
                    SaveResult(model, filelist, folder.toPath().toString()); //.xls
                }
            }
        } else {
            System.out.println("GeneDistance (2025) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/GeneDistance\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar GeneDistance.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("mod=1\t Different sets of kmers used in the analysis: mod=0...4 (default mod=0)");
        }
    }

    private static void SaveResult(int model, String[] infiles, String folder) {
        try {
            long startTime = System.nanoTime();
            System.out.println("Running...");
            SequencesSimilarity s1 = new SequencesSimilarity(model);
            s1.SetFolder(infiles, folder);
            s1.RunPattern();
            // s1.RunKmerCounter();
            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds");
        } catch (IOException e) {
            System.out.println("Incorrect file name.");
        }
    }

    public static int StrToInt(String str) {
        StringBuilder r = new StringBuilder();
        int z = 0;
        r.append(0);
        for (int i = 0; i < str.length(); i++) {
            char chr = str.charAt(i);
            if (chr > 47 && chr < 58) {
                r.append(chr);
                z++;
                if (z > 10) {
                    break;
                }
            }
            if (chr == '.' || chr == ',') {
                break;
            }
        }
        return (Integer.parseInt(r.toString()));
    }
}
