import java.io.IOException;
import java.io.File;

public class KmerDistance {

    public static void main(String[] args) {
        if (args.length > 0) {
            String infile = args[0]; // file path or Folder

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            File folder = new File(infile);
            if (folder.exists() && (folder.isDirectory() || folder.isFile())) {
                if (folder.isDirectory()) {
                    File[] files = folder.listFiles();
                    int k = -1;
                    String[] filelist = new String[files.length];
                    for (File file : files) {
                        if (file.isFile()) {
                            filelist[++k] = file.getAbsolutePath();
                        }
                    }
                    SaveResult(filelist, folder.toPath().toString() + File.separator + "result.xls");
                } else {
                    String[] filelist = new String[1];
                    filelist[0] = infile;
                    SaveResult(filelist, folder.toPath().toString() + ".xls"); //.xls
                }
            }
        } else {
            System.out.println("KmersPattern (2024) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/KmersPattern\n");
        }
    }

    private static void SaveResult(String[] infiles, String folder) {
        try {
            long startTime = System.nanoTime();
            System.out.println("Running...");
            SequencesSimilarity s1 = new SequencesSimilarity();
            s1.SetFolder(infiles, folder);
            s1.RunPattern();
           // s1.RunKmerCounter();
            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds");
        } catch (IOException e) {
            System.out.println("Incorrect file name.");
        }
    }
}
