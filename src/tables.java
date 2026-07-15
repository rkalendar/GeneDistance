
public final class tables {

    static final public byte[] cdn = new byte[128];

    static {
// "M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T and I"
// maps every nucleotide symbol, upper or lower case, to its lower-case code; 0 for everything else
//a   b   c   d   g   h   i   k   m   n   r   s   t   u   v   w   y
//97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 121
        cdn[65] = 97;   // A
        cdn[66] = 98;   // B
        cdn[67] = 99;   // C
        cdn[68] = 100;  // D
        cdn[71] = 103;  // G
        cdn[72] = 104;  // H
        cdn[73] = 99;   // I
        cdn[75] = 107;  // K
        cdn[77] = 109;  // M
        cdn[78] = 110;  // N
        cdn[82] = 114;  // R
        cdn[83] = 115;  // S
        cdn[84] = 116;  // T
        cdn[85] = 116;  // U
        cdn[86] = 118;  // V
        cdn[87] = 119;  // W
        cdn[89] = 121;  // Y
        cdn[97] = 97;   // a
        cdn[98] = 98;    // b
        cdn[99] = 99;    // c
        cdn[100] = 100;  // d
        cdn[103] = 103;  // g
        cdn[104] = 104;  // h
        cdn[105] = 99;   // i
        cdn[107] = 107;  // k
        cdn[109] = 109;  // m
        cdn[110] = 110;  // n
        cdn[114] = 114;  // r
        cdn[115] = 115;  // s
        cdn[116] = 116;  // t
        cdn[117] = 116;  // u
        cdn[118] = 118;  // v
        cdn[119] = 119;  // w
        cdn[121] = 121;  // y
    }
}
