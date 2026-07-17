import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * A segmented, in-memory view over all the sequence data, built by streaming the
 * files. The bases are held in a list of byte[] segments addressed by a global
 * {@code long} offset, so the whole thing may exceed the ~2.1 GB a single Java
 * array or String can hold, and each record is reached by its global offset.
 *
 * <p>It exists for the modes that need RANDOM access to the data - {@code -vector}
 * compares a short sequence to arbitrary windows of a longer one and draws random
 * windows to calibrate its chance level - which a forward stream cannot serve. The
 * bytes are exactly the ones {@link StreamingFasta#forEachRecord} yields, so any
 * window sliced from the store equals the same substring of the record's String,
 * and the modes reading from it produce identical output.
 */
final class SeqStore {

    private static final int SEG = 1 << 30;   // 1 GB segments: a record may span several

    private final List<byte[]> segs = new ArrayList<>();
    private byte[] cur = new byte[Math.min(1 << 20, SEG)];
    private int curLen;
    private long total;

    private final List<String> names = new ArrayList<>();
    private final List<Long> offs = new ArrayList<>();
    private final List<Long> lens = new ArrayList<>();

    static SeqStore build(String[] files) throws IOException {
        final SeqStore st = new SeqStore();
        StreamingFasta.forEachRecord(files, st::add);
        st.seal();
        return st;
    }

    private void add(String name, String s) {
        names.add(name);
        offs.add(total);
        final int n = s.length();
        lens.add((long) n);
        for (int i = 0; i < n; i++) {
            if (curLen == cur.length) {
                grow();
            }
            cur[curLen++] = (byte) s.charAt(i);   // bases are ISO-8859-1 (a..z), one byte each
        }
        total += n;
    }

    private void grow() {
        if (cur.length < SEG) {
            byte[] bigger = new byte[Math.min(SEG, cur.length * 2)];
            System.arraycopy(cur, 0, bigger, 0, curLen);
            cur = bigger;
        } else {                       // segment full: push it and start a fresh one
            segs.add(cur);
            cur = new byte[Math.min(1 << 20, SEG)];
            curLen = 0;
        }
    }

    private void seal() {
        final byte[] last = new byte[curLen];
        System.arraycopy(cur, 0, last, 0, curLen);
        segs.add(last);
        cur = null;
    }

    int count() {
        return names.size();
    }

    String name(int j) {
        return names.get(j);
    }

    long length(int j) {
        return lens.get(j);
    }

    long offset(int j) {
        return offs.get(j);
    }

    /** The whole record j as a String (requires its length to fit ~2.1 GB). */
    String record(int j) {
        return substring(offs.get(j), offs.get(j) + lens.get(j));
    }

    /** The window of {@code len} bases at {@code pos} within record j. */
    String window(int j, long pos, int len) {
        final long from = offs.get(j) + pos;
        return substring(from, from + len);
    }

    /** The bases of the global range [from, to) as a String (the range must fit an int). */
    String substring(long from, long to) {
        final int len = (int) (to - from);
        final byte[] out = new byte[len];
        long g = from;
        int w = 0;
        while (w < len) {
            final int seg = (int) (g / SEG);
            final int within = (int) (g % SEG);
            final byte[] s = segs.get(seg);
            final int take = Math.min(len - w, s.length - within);
            System.arraycopy(s, within, out, w, take);
            w += take;
            g += take;
        }
        return new String(out, StandardCharsets.ISO_8859_1);
    }
}
