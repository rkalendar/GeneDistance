import java.nio.charset.StandardCharsets;

/**
 * Reads the sequences of a FASTA source held in memory.
 *
 * A record opens on a '>' that begins a line; the rest of that line is its
 * name, and everything up to the next record is its sequence, from which only
 * the nucleotide symbols are kept, in lower case (the table in {@link tables}
 * maps them). Digits, gaps, spaces and line breaks are dropped, so a sequence
 * numbered or wrapped in any way is read all the same.
 */
public final class ReadingSequencesFiles {

    /**
     * The bytes are filtered in place, so the array given is modified and must
     * not be reused by the caller.
     */
    public ReadingSequencesFiles(byte[] s) {
        int n = 0;
        for (int i = 0; i < s.length; i++) {
            if (s[i] > 8 && s[i] < 128) {   // keeps the line breaks, drops control and non-ASCII bytes
                s[n++] = s[i];
            }
        }
        source = s;
        size = n;
        ReadingSequences();
        source = null;   // the source is of no use once the sequences are read
    }

    public ReadingSequencesFiles(String s) {
        this(s.getBytes(StandardCharsets.ISO_8859_1));   // through the same filtering as a file
    }

    /** The sequences, in the order they were read; empty when the source holds none. */
    public String[] getSequences() {
        return sequence;
    }

    /** The name of each sequence: its FASTA header, without the '>'. */
    public String[] getNames() {
        return name_seq;
    }

    public int getNseq() {
        return ns;
    }

    private void ReadingSequences() {
        final byte[] src = source;
        final int l = size;
        final byte[] cdn = tables.cdn;

        int letters = 0;
        int records = 0;
        for (int i = 0; i < l; i++) {
            final byte b = src[i];
            if (b == '>') {
                if (i == 0 || src[i - 1] == '\n' || src[i - 1] == '\r') {
                    records++;   // a header is a line of its own: a '>' inside one opens nothing
                }
            } else if (cdn[b] > 0) {
                letters++;
            }
        }
        if (letters == 0 || records == 0) {
            ns = 0;
            return;
        }

        name_seq = new String[records];
        sequence = new String[records];

        int n = -1;
        int t = 0;   // where the sequence of the record being read starts
        for (int i = 0; i < l; i++) {
            if (src[i] == '>' && (i == 0 || src[i - 1] == '\n' || src[i - 1] == '\r')) {
                if (n >= 0) {
                    sequence[n] = read(t, i);
                }
                n++;
                int j = i + 1;
                while (j < l && src[j] != '\n' && src[j] != '\r') {
                    j++;
                }
                // a header the source ends on, with no line break, still names its (empty) sequence
                name_seq[n] = new String(src, i + 1, j - i - 1, StandardCharsets.ISO_8859_1).trim();
                i = j;
                t = j + 1;
            }
        }
        sequence[n] = read(Math.min(t, l), l);
        ns = n + 1;
    }

    /**
     * Reads the sequence of source[from,to): keeps the nucleotide symbols, in
     * lower case, and counts them. They are packed back into that same range,
     * which the reading has left behind and no one reads again, so a sequence
     * costs one String and nothing else.
     */
    private String read(int from, int to) {
        final byte[] src = source;
        final byte[] cdn = tables.cdn;

        int x = from;
        for (int i = from; i < to; i++) {
            final byte c = cdn[src[i]];
            if (c > 0) {
                src[x++] = c;
            }
        }
        return new String(src, from, x - from, StandardCharsets.ISO_8859_1);
    }

    private String[] name_seq = new String[0];
    private String[] sequence = new String[0];
    private byte[] source;
    private int size;
    private int ns = 0;
}
