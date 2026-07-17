import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Streams a FASTA file record by record, so a multi-FASTA larger than 2 GB - one
 * that {@link java.nio.file.Files#readAllBytes} cannot hold in a single byte[] -
 * is read without ever materialising the whole file.
 *
 * <p>The byte filtering matches {@link ReadingSequencesFiles} exactly: bytes
 * outside 9..127 are dropped, a record opens on a '>' that begins a line, its
 * name is the rest of that line (trimmed), and its sequence keeps only the
 * nucleotide symbols, in lower case, mapped through {@link tables#cdn}. Given a
 * file with any nucleotide, every record is kept - including empty ones - as the
 * original does; a file with none yields nothing. Verified byte for byte against
 * {@link ReadingSequencesFiles} over the test corpus and its edge cases.
 *
 * <p>Three ways to consume it, each holding at most one record at a time rather
 * than the file:
 * <ul>
 * <li>{@link #forEachRecord} - the name and the sequence as a String, for the
 * modes that keep every sequence in memory (composition, -scan subjects);
 * <li>{@link #forEachSketch} - the name, length and FracMinHash sketch, rolled
 * over the bases without ever building the sequence String, for -contain at
 * genome scale in constant memory.
 * </ul>
 */
final class StreamingFasta {

    private static final int BUF = 1 << 20;

    /** 2-bit code of a lower-case base (what {@code tables.cdn} produces); -1 otherwise. */
    private static final int[] BASE = new int[128];

    static {
        Arrays.fill(BASE, -1);
        BASE['a'] = 0;
        BASE['c'] = 1;
        BASE['g'] = 2;
        BASE['t'] = 3;
    }

    private StreamingFasta() {
    }

    /** Called with every parsed record; the sequence is discarded once it returns. */
    interface RecordHandler {
        void accept(String name, String seq);
    }

    /** Called with every record's FracMinHash sketch, the sequence never built. */
    interface SketchHandler {
        void accept(String name, long length, long[] sketch);
    }

    /**
     * Streams every record of every file, one materialised (name, sequence) at a
     * time. A file with no nucleotide contributes nothing, exactly as a separate
     * {@link ReadingSequencesFiles} call would.
     */
    static void forEachRecord(String[] files, RecordHandler h) throws IOException {
        for (String f : files) {
            if (f == null) {
                continue;
            }
            RecordSink sink = new RecordSink(h);
            read(f, sink);
            sink.endFile();
        }
    }

    /**
     * The length of every record, in one streaming pass, without building any
     * sequence - so the k and scaled of a containment run can be chosen for the
     * longest genome before a single sketch is built.
     */
    static long[] recordLengths(String[] files) throws IOException {
        final List<Long> lens = new ArrayList<>();
        for (String f : files) {
            if (f == null) {
                continue;
            }
            LengthSink sink = new LengthSink(lens);
            read(f, sink);
            sink.endFile();
        }
        final long[] out = new long[lens.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = lens.get(i);
        }
        return out;
    }

    /**
     * Streams every record's canonical FracMinHash sketch without building the
     * sequence, so a 2 Gb genome costs only its sketch. {@code maxHash == -1L}
     * keeps every k-mer (the exact set).
     */
    static void forEachSketch(String[] files, int k, long maxHash, SketchHandler h) throws IOException {
        for (String f : files) {
            if (f == null) {
                continue;
            }
            SketchSink sink = new SketchSink(k, maxHash, h);
            read(f, sink);
            sink.endFile();
        }
    }

    // ── streaming core ──────────────────────────────────────────────────────────
    /**
     * Drives {@code sink} over one file. Bytes outside 9..127 are dropped before
     * anything else (as {@link ReadingSequencesFiles}'s in-place pre-filter does),
     * so line-start detection and header parsing see the same stream the original
     * did.
     */
    private static void read(String file, Sink sink) throws IOException {
        final byte[] cdn = tables.cdn;
        final byte[] buf = new byte[BUF];
        int prev = -1;          // last kept byte, for line-start detection
        boolean inName = false; // reading a header line
        boolean open = false;   // a record has begun
        try (InputStream in = new BufferedInputStream(new FileInputStream(file), BUF)) {
            int n;
            while ((n = in.read(buf, 0, BUF)) > 0) {
                for (int p = 0; p < n; p++) {
                    final int rb = buf[p] & 0xFF;
                    if (rb <= 8 || rb >= 128) {
                        continue;               // dropped by the pre-filter, does not affect line starts
                    }
                    final boolean lineStart = (prev == -1 || prev == '\n' || prev == '\r');
                    if (rb == '>' && lineStart) {
                        if (open) {
                            sink.endRecord();
                        }
                        sink.startRecord();
                        open = true;
                        inName = true;
                    } else if (inName) {
                        if (rb == '\n' || rb == '\r') {
                            inName = false;
                            sink.endName();
                        } else {
                            sink.nameByte(rb);
                        }
                    } else if (open) {
                        // bases before the first '>' belong to no record and are dropped,
                        // exactly as ReadingSequencesFiles discards everything before it
                        final int c = cdn[rb] & 0xFF;   // cdn is a byte table; a kept base is 97..121
                        if (c > 0) {
                            sink.base(c);
                        }
                    }
                    prev = rb;
                }
            }
        }
        if (open) {
            sink.endRecord();
        }
    }

    private interface Sink {
        void startRecord();

        void nameByte(int b);

        void endName();

        void base(int c);

        void endRecord();
    }

    /** Builds each record's name and sequence String, then hands it on and drops it. */
    private static final class RecordSink implements Sink {
        private final RecordHandler h;
        private final List<String> names = new ArrayList<>();
        private final List<String> seqs = new ArrayList<>();
        private byte[] nameBuf = new byte[128];
        private int nameLen;
        private byte[] seqBuf = new byte[1 << 16];
        private int seqLen;
        private String pendingName;
        private long fileLetters;

        RecordSink(RecordHandler h) {
            this.h = h;
        }

        @Override
        public void startRecord() {
            nameLen = 0;
            pendingName = "";
        }

        @Override
        public void nameByte(int b) {
            if (nameLen == nameBuf.length) {
                nameBuf = Arrays.copyOf(nameBuf, nameBuf.length * 2);
            }
            nameBuf[nameLen++] = (byte) b;
        }

        @Override
        public void endName() {
            pendingName = new String(nameBuf, 0, nameLen, StandardCharsets.ISO_8859_1).trim();
        }

        @Override
        public void base(int c) {
            if (seqLen == seqBuf.length) {
                seqBuf = Arrays.copyOf(seqBuf, seqBuf.length * 2);
            }
            seqBuf[seqLen++] = (byte) c;
        }

        @Override
        public void endRecord() {
            names.add(pendingName);
            seqs.add(new String(seqBuf, 0, seqLen, StandardCharsets.ISO_8859_1));
            fileLetters += seqLen;
            seqLen = 0;
        }

        void endFile() {
            if (fileLetters > 0) {   // a file with no nucleotide yields nothing, as the original does
                for (int i = 0; i < names.size(); i++) {
                    h.accept(names.get(i), seqs.get(i));
                }
            }
            names.clear();
            seqs.clear();
        }
    }

    /**
     * Rolls the canonical k-mer over each record and keeps the FracMinHash sketch,
     * so the sequence is never held. Reproduces {@code SequencesSimilarity.canonical}
     * bit for bit - same rolling, same {@code min(fwd, rev)}, same hash filter.
     */
    private static final class SketchSink implements Sink {
        private final int k;
        private final long maxHash;
        private final boolean sketch;
        private final SketchHandler h;
        private final long mask;
        private final int top;

        private byte[] nameBuf = new byte[128];
        private int nameLen;
        private String pendingName;
        private long fileLetters;

        private long[] buf = new long[1024];
        private int m;
        private long len;
        private long fwd, rev;
        private int valid;
        private final List<String> names = new ArrayList<>();
        private final List<long[]> sketches = new ArrayList<>();
        private final List<Long> lengths = new ArrayList<>();

        SketchSink(int k, long maxHash, SketchHandler h) {
            this.k = k;
            this.maxHash = maxHash;
            this.sketch = (maxHash != -1L);
            this.h = h;
            this.mask = (k >= 32) ? -1L : ((1L << (2 * k)) - 1);
            this.top = 2 * (k - 1);
        }

        @Override
        public void startRecord() {
            nameLen = 0;
            pendingName = "";
            resetRoll();
        }

        private void resetRoll() {
            m = 0;
            len = 0;
            fwd = 0;
            rev = 0;
            valid = 0;
        }

        @Override
        public void nameByte(int b) {
            if (nameLen == nameBuf.length) {
                nameBuf = Arrays.copyOf(nameBuf, nameBuf.length * 2);
            }
            nameBuf[nameLen++] = (byte) b;
        }

        @Override
        public void endName() {
            pendingName = new String(nameBuf, 0, nameLen, StandardCharsets.ISO_8859_1).trim();
        }

        @Override
        public void base(int c) {
            len++;
            final int b = (c < 128) ? BASE[c] : -1;
            if (b < 0) {
                valid = 0;
                fwd = 0;
                rev = 0;
                return;
            }
            fwd = ((fwd << 2) | b) & mask;
            rev = (rev >>> 2) | ((long) (3 - b) << top);
            if (valid < k) {
                valid++;
            }
            if (valid == k) {
                final long can = Math.min(fwd, rev);
                if (!sketch || Long.compareUnsigned(mix64(can), maxHash) <= 0) {
                    if (m == buf.length) {
                        buf = Arrays.copyOf(buf, buf.length * 2);
                    }
                    buf[m++] = can;
                }
            }
        }

        @Override
        public void endRecord() {
            long[] r = Arrays.copyOf(buf, m);
            Arrays.sort(r);
            int d = 0;
            for (int i = 0; i < r.length; i++) {
                if (i == 0 || r[i] != r[i - 1]) {
                    r[d++] = r[i];
                }
            }
            names.add(pendingName);
            sketches.add(Arrays.copyOf(r, d));
            lengths.add(len);
            fileLetters += len;
            resetRoll();
        }

        void endFile() {
            if (fileLetters > 0) {
                for (int i = 0; i < names.size(); i++) {
                    h.accept(names.get(i), lengths.get(i), sketches.get(i));
                }
            }
            names.clear();
            sketches.clear();
            lengths.clear();
        }
    }

    /** Counts each record's bases, per file, dropping a file that holds none. */
    private static final class LengthSink implements Sink {
        private final List<Long> out;
        private final List<Long> file = new ArrayList<>();
        private long len;
        private long fileLetters;

        LengthSink(List<Long> out) {
            this.out = out;
        }

        @Override
        public void startRecord() {
            len = 0;
        }

        @Override
        public void nameByte(int b) {
        }

        @Override
        public void endName() {
        }

        @Override
        public void base(int c) {
            len++;
        }

        @Override
        public void endRecord() {
            file.add(len);
            fileLetters += len;
        }

        void endFile() {
            if (fileLetters > 0) {
                out.addAll(file);
            }
            file.clear();
        }
    }

    /** SplitMix64 finalizer - identical to SequencesSimilarity.mix64. */
    private static long mix64(long z) {
        z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L;
        z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL;
        return z ^ (z >>> 31);
    }
}
