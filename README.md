# GeneDistance

**Identification of genetic similarity, or the distance, between genomic sequences (approximate string matching)**

[![Java](https://img.shields.io/badge/Java-26+-orange.svg)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-blue.svg)]()
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)

Six measures over one input, each answering a different question: are these two
sequences alike as a whole (`-vector`, `-d2star`, `-cosine`, `-spacing`), how much
of one occurs in the other (`-contain`), and where (`-scan`).

## Contents

- [Requirements](#requirements) · [Installation](#installation) · [Basic usage](#basic-usage) · [Options](#options)
- [The spacing signature: `-spacing`](#the-spacing-signature--spacing)
- [Sequences of unequal length](#sequences-of-unequal-length)
- [Fragments of very different length: `-vector`](#fragments-of-very-different-length--vector)
- [A short element inside a long sequence](#a-short-element-inside-a-long-sequence)
- [Chromosome-scale input (up to 2 Gb): `-scaled`](#chromosome-scale-input-up-to-2-gb--scaled)
- [FASTA larger than 2 GB](#fasta-larger-than-2-gb)
- [Where the element sits](#where-the-element-sits)
- [Output](#output) · [Examples](#examples)

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi


## Requirements

| Requirement | Details |
|-------------|---------|
| **Java** | Version 26 or higher |
| **OS** | Windows, Linux, or macOS |
| **RAM** | Default JVM heap is sufficient for most tasks |

- **Download Java:** <https://www.oracle.com/java/technologies/downloads/>
- **Set Java Path:** <https://www.java.com/en/download/help/path.html>

---

## Installation

### Option 1: Direct Download

1. Download `GeneDistance.jar` from the `dist` directory of this repository.
2. Place it in your preferred location.
3. Ensure Java 26+ is installed and available in your `PATH`. How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)

Verify Java is available:

```bash
java -version
```

### Option 2: Install Java via Conda

```bash
# Add conda-forge channel and set priority
conda config --add channels conda-forge
conda config --set channel_priority strict

# Create environment with OpenJDK 26
conda create -n java26 openjdk=26

# Activate environment
conda activate java26

# Verify installation
java -version
```


To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```GeneDistance.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar GeneDistance.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <GeneDistancePath>\dist\GeneDistance.jar <target_file_path> optional_commands```

The full list of the options is printed by:

```java -jar GeneDistance.jar --help```


### Options:

The **default measure is ```-vector```** — run with no measure flag and it is used. The former
default, the k-mer spacing measure, is now ```-spacing```.

| Option | Meaning |
| --- | --- |
| ```-kmer=N``` | set of k-mers used in the analysis: ```4```/```5```/```6```/```8```/```10``` or ```41```/```61```/```81```/```101```, default ```-kmer=4```; ```-vector``` uses ```-kmer=41``` unless a set is named (see [what the numbers mean](#what-the--kmer-numbers-mean)) |
| ```-vector``` | **the default** — a method for homology between fragments of very different length. A short sequence is compared to the **best window** of a longer one at its own scale, on **either strand**, and a score counts only when it clears the **chance level measured from the input itself**, wherever one is measured for that length class. The per-window comparison of two composition profiles is equivalent to the *d2\** statistic (k-mer counts centred on the base composition); the local window, the two strands and the calibrated floor are ```-vector```'s own — none of them part of *d2\** |
| ```-spacing``` | the k-mer **spacing signature**: the sequences are compared through the ratios of the average spacings of the k-mers of the set. Independent of sequence length and of strand; it tells whether two sequences have the same *composition* (the same genome), not whether they are homologous (see [the spacing signature](#the-spacing-signature--spacing)) |
| ```-d2star``` | compare the k-mer **frequencies** instead of their spacing, each count centred on the base composition of its own sequence (the *d2\** measure) |
| ```-cosine``` | the same, on the raw frequencies, without centring |
| ```-contain```, ```-contain=K``` | **containment**: how much of each sequence occurs in each other one, on the whole k-mer space. Finds a short element inside a long sequence. ```-contain=K``` sets the exact k in one go, as ```-ksize=K``` does |
| ```-scan```, ```-scan=N``` | **where** each sequence holds the query: a window the size of the query is walked along them, and every region that holds it is reported with its position, identity, number of separate matches and E-value. The query is the first sequence of the input, or the N-th |
| ```-ksize=K``` | the k of the exact k-mers ```-contain``` and ```-scan``` use (6-31, 11-16 suits most cases). Chosen for the length of the sequences unless given. It is unrelated to ```-kmer=```, which names one of the sets above |
| ```-scaled=N``` | **FracMinHash** for ```-contain``` / ```-scan```: keep only ~1/N of the k-mers so a 2 Gb genome fits in memory (N = 1…1000000). ```N=1``` is the default and is exact; ```-scaled``` or ```-scaled=auto``` chooses N for the size |
| ```-kmerstat``` | report the k-mer distances of each sequence, instead of comparing them |
| ```-kmer2stat``` | report the k-mer distances averaged over all the target sequences |
| ```-h```, ```--help```, ```-?```, ```/?``` | print the help and exit |
| ```--version``` | print the version and exit |

```-d2*``` is accepted for ```-d2star``` and ```-space``` for ```-spacing```. Option names are
case-insensitive and any number of leading dashes is accepted (```-vector```, ```--vector```), but
an option must carry at least one — ```kmer=N``` and ```/?``` are the only exceptions, anything else
without a dash is taken for the input path. **An option that is not recognised is ignored without a
message**, so a misspelt measure leaves the default one running: check the ```Measure:``` line each
run prints before trusting a result.

The first argument that is not an option is the target file or folder; giving a second one is an
error. Given a folder, every file directly in it is read — sub-folders are not descended into, and
the tool's own ```result*.xls``` / ```result*.meg``` reports are skipped so that a re-run never
takes them back in as input — and all their sequences are pooled into **one** analysis.

### What the ```-kmer=``` numbers mean

| Value | k-mer set | Value | k-mer set |
| --- | --- | --- | --- |
| ```4``` (default) | symmetric and mirror 4-mers | ```41``` | non-symmetric 4-mers |
| ```5``` | symmetric 5-mers | ```61``` | non-symmetric 6-mers |
| ```6``` | symmetric 6-mers | ```81``` | non-symmetric 8-mers |
| ```8``` | symmetric 8-mers | ```101``` | non-symmetric 10-mers |
| ```10``` | symmetric 10-mers | | |

The two-digit forms are words of the same length, but a **separate list** — not the plain set
counted a different way. The plain ```4```/```6```/```8```/```10``` sets hold only words that are
their own reverse complement (```5```, being of odd length, holds instead the 64 words whose flanks
are reverse-complementary around any centre); the two-digit sets drop that requirement. So ```41``` (242 of the 256 4-mers) and ```61``` (109) contain their plain
counterparts, while ```81``` (180) and ```101``` (900) are built of words that read the same
backwards and share no word at all with ```8``` (252) and ```10``` (900). ```-d2star```,
```-cosine``` and ```-vector``` count each word of the set on the strand given — a word outside the
set in force is simply not counted. ```-spacing``` pools every word with its reverse complement, so
that its result is the same on either strand: a palindromic set is unchanged by that, while ```41```
becomes 127 features, ```61``` 85 and ```5``` 32 (the run prints the number). The set in force is
printed at the start of every run that uses one. ```-kmer=``` does not apply to ```-contain``` and ```-scan```, which use
exact k-mers of length ```-ksize=K``` and print that ```k``` instead.


### The spacing signature: ```-spacing```

Each sequence is described by the average spacing of every k-mer of the set, *x<sub>i</sub> = L / q<sub>i</sub>*
(length over count, for k-mers seen at least twice), and two sequences are compared through the
ratios of those spacings, pair of k-mers by pair of k-mers:

```
r_ij = x_i / x_j = q_j / q_i                      inside one sequence: the length L cancels
F_ij = min(r_ij(A), r_ij(B)) / max(r_ij(A), r_ij(B))
S    = mean of F_ij over all pairs i < j          reported as 100*S; the .meg holds 1 - S
```

With *d<sub>i</sub> = ln q<sub>i</sub>(A) − ln q<sub>i</sub>(B)*, *F<sub>ij</sub> = exp(−|d<sub>i</sub> − d<sub>j</sub>|)*: the
measure is the spread of the log fold-changes of the k-mer counts between the two sequences. It
follows that

- **length does not matter**: *S* = 100% whenever the two count vectors are proportional (a sequence
  and a tandem copy of itself score 100%);
- **strand does not matter**: every word is pooled with its reverse complement (a palindromic word
  is its own), so a sequence and its reverse complement score 100%;
- every pair of k-mers weighs the same, and *S* is computed exactly in *O(n log n)* by sorting the
  *d<sub>i</sub>*, not in *O(n²)*;
- it measures **composition, not homology**. Two unrelated sequences of the same composition tend
  to 100% as they grow, and below a certain length the counts are too small to tell anything: two
  unrelated random sequences score about 60–70% at 2–10 kb. Use it on long sequences — roughly
  50 kb and more for the 4-mer sets, 1 Mb and more for the 6-mers — to ask *which genome* a
  chromosome or a scaffold belongs to. On 360 chromosomes of 24 animal and plant species,
  ```-spacing -kmer=6``` puts 93% of them next to a chromosome of their own species.

### Sequences of unequal length

```-spacing``` only counts a k-mer once it has been seen at least twice, so a short sequence arrives
with an almost empty vector, and its homology to a long one goes unseen. That is why it is not the
default: use ```-d2star``` there, or ```-vector``` (the default) when the lengths differ by more than
a few fold. On ```test/2.txt``` (18S rRNA of macaque and man, the same 18S diluted with random DNA,
and three random sequences), with ```-kmer=41```:

| pair | truth | ```-spacing``` | ```-cosine``` | ```-d2star``` |
| --- | --- | --- | --- | --- |
| 18S macaque / 18S human | orthologs | 99 | 100 | 100 |
| 18S / the 18S-28S unit holding it | homologous | 56 | 83 | 51 |
| 18S / 18S + 2 kb of random DNA | homologous | 81 | 96 | 75 |
| 18S / random DNA | **unrelated** | **64** | **87** | **6** |
| random / random | **unrelated** | **68** | **89** | **2** |

The reports of the measures are written side by side, as ```_kN.*```, ```_kN_cos.*```,
```_kN_d2s.*```, ```_kN_vec.*```, ```_cK.*``` and ```_scanK.xls```, and never overwrite one another.

*d2\**: Reinert G, Chew D, Sun F, Waterman MS (2009) *J Comput Biol* 16:1615-1634; Song K et al.
(2014) *Brief Bioinform* 15:343-353.


### Fragments of very different length: ```-vector```

```-d2star``` compares whole sequences, which averages a shared region away when one sequence is
much longer than the other. ```-vector``` fixes that the way TotalRepeats' *-vector* mode does: a
short sequence is compared not to a long one as a whole but to the **best window of the long one at
the short one's own scale** (step = half a window), on **either strand**, and, wherever that length
scale has a chance level to measure, a score is kept only when it stands clear of the **chance level
measured for that scale** from the input's own content (non-homologous windows screened by an exact
25-mer). On ```test/2.txt``` with ```-kmer=41```:

| pair | truth | ```-d2star``` | ```-vector``` |
| --- | --- | --- | --- |
| 18S / the 18S-28S unit holding it | 18S is inside it | 51 | **97** |
| 28S / the unit (reverse-complemented inside) | 28S is inside it | 45 | **98** |
| 18S / 18S + 8.6 kb random | homologous | 29 | 100 |
| any / random, barley, arabidopsis | **unrelated** | small | **0** |

```-vector``` works on the 4-mer set ```41``` unless ```-kmer=``` names another: the twelve
palindromic 4-mers of the general default are too few to compare the composition of windows.

The chance level of a length scale is not the highest score the sampled windows reach, since copies
of one element that have diverged by more than ~10% share no 25-mer, pass the screen, and would set
the level to their own score. It is placed where the largest of that many chance scores is expected,
from the median of the sample and its spread *below* the median, which homologues do not reach (and
never above the highest score seen).

```-vector``` is a measure for **near-identical fragments**. On families made from real
transposable-element consensus sequences (truncated copies, either strand) it assigns 100% of the
copies to their family at 2% divergence and 99% at 5%, but copies that have diverged by 15% or more
mostly fall below the chance level, because a 4-mer composition score of such copies is no longer
clear of what unrelated windows reach. Use ```-contain``` for those: it keeps 99.5% at 15%.

The floor is calibrated for length scales up to ~6 kb; above that the chance level of a composition
score is already ~0, so the class is deliberately left uncut. A scale that *cannot* be calibrated —
too little content, or too few non-homologous windows to draw from — is likewise left uncut. Either
way its raw, already well-separated scores are reported. For a short element that shares only part of *both* sequences (e.g. two different genomic
regions each carrying the same element plus their own flanks), ```-contain``` remains the sharper
tool — it asks "how much is shared" rather than "is the composition matched at scale".

Inside a single window, the two composition profiles are scored the way the *d2\** statistic does —
k-mer counts centred on the base composition each window would produce by chance, then compared as a
correlation (Reinert et al., *J Comput Biol*, 2009). That is the kernel only; ```-vector``` is not
*d2\**. *d2\** is a whole-sequence distance, whereas ```-vector``` wraps that kernel in the local
window, the two strands and the data-calibrated floor above — which is what lets it find partial
homology across very different lengths, where whole-sequence *d2\** cannot.

*Windowing and the per-scale chance floor follow* SequencesClustering *of* TotalRepeats.


### A short element inside a long sequence

No measure of the composition of a sequence *as a whole* - not the default one, not *d2\** - can
find a 300 bp element inside a 10 kb sequence: the element is 3% of what is being measured. That
is what ```-contain``` is for. It asks a different question, and it is asymmetric:

```
C(A,B) = |kmers(A) & kmers(B)| / |kmers(A)|      how much of A occurs in B
```

An Alu element (289 bp), and the same element planted, 5% mutated, in the middle of 10 kb of
unrelated DNA:

| | ```-d2star``` | ```-contain``` |
| --- | --- | --- |
| Alu / the 10 kb **holding it** | 42 | **75** |
| Alu / 10 kb of unrelated DNA | 39 | **0** |

The k-mers are canonical, so an element inserted the other way round is still found - in
```test/2.txt``` that is how the 28S turns up inside the 18S-28S unit, reverse complemented, at 98%.
The report also estimates the identity of the shared part, *p = C*<sup>1/k</sup>, and leaves it at 0
where the shared k-mers are no more than chance would give.

```-contain``` keeps the k-mers of every sequence in memory, about 8 bytes per base, so a folder of
whole genomes wants ```java -Xmx8g -jar ...``` — or a **FracMinHash sketch** to bound that.


### Chromosome-scale input (up to 2 Gb): ```-scaled```

The exact ```-contain``` / ```-scan``` hold every distinct k-mer, ~8 bytes each, so a 2 Gb genome
needs ~16 GB per sequence. ```-scaled=N``` keeps only the k-mers whose hash falls in the bottom
1/N of the hash space (a FracMinHash sketch, as in Mash / sourmash), so the set shrinks by N while
the containment estimated from the kept fraction stays unbiased. Measured on two related 50 Mb
sequences:

| | heap needed | time | containment | identity |
| --- | --- | --- | --- | --- |
| ```-contain``` (exact) | ~2 GB | 8.6 s | 8% | 86% |
| ```-contain -scaled=1000``` | **~200 MB** | **1.0 s** | 8% | 86% |

Same numbers, a tenth of the memory and a tenth of the time. At ```-scaled=1000``` a 2 Gb genome's
sketch is ~2 M hashes (~16 MB) instead of ~16 GB. ```-scaled=1``` (the default) is the exact
measure, bit for bit; ```-scaled``` with no number picks N from the **longest sequence** in the
input — about one hash kept per 2 Mb of it, so a 2 Gb chromosome gives N≈1000, and an input whose
longest sequence is under ~3 Mb stays exact at N=1. Reports go to
```_cKsN.*``` / ```_scanKsN.*``` so a sketch run never overwrites an exact one.

*FracMinHash*: Irber L et al. (2022) *sourmash*; Hera MR, Pierce-Ward NT, Koslicki D (2023)
*Genome Research* 33(7):1061-1068.


### FASTA larger than 2 GB

The input is read **record by record from a stream**, not loaded whole, so a multi-FASTA past the
~2.1 GB that a single Java array can hold is read all the same. The records are byte for byte the
ones the old whole-file reader produced (verified over the test corpus and its edge cases), so every
mode's output is unchanged.

What that costs in memory depends on the mode:

| mode | held in memory | a 2 Gb genome |
| --- | --- | --- |
| ```-d2star``` / ```-cosine``` | one small vector per record, streamed | **constant** — ran 2×50 Mb in a 256 MB heap |
| ```-contain -scaled=N``` | only the sketches | **~16 MB per genome** — 512 MB heap is plenty |
| ```-scan``` | the query index and one subject at a time | ~ the largest record |
| ```-vector``` (default) | a segmented ```SeqStore``` (it needs random access to windows) | ~ the data, but a record may exceed 2.1 GB |
| ```-spacing``` | one small vector per record, streamed | **constant** |

```-d2star```/```-cosine```/```-spacing``` only pass over each record once, so they stream it into a
small vector and never hold the sequences. ```-vector``` (the default) compares a short sequence to
arbitrary windows of a longer one and draws random windows to calibrate, so it needs random access; it
streams the input into a **segmented ```SeqStore```** (1 GB byte[] segments addressed by a global
```long``` offset), which holds the data once — no second copy for the reverse strand — and lets a
single chromosome exceed the 2.1 GB a String can hold. Every mode now reads the input as a stream.

Demonstrated on a 2.09 GiB, 3 × 740 Mb file: ```Files.readAllBytes``` cannot allocate an array that
large and the old reader fails with an out-of-memory error even at ```-Xmx8g```; the streaming
reader runs ```-contain -scaled=1000``` on it in a 512 MB heap in 35 s. One record may still exceed
~2.1 Gb only in the exact (unsketched) modes that hold it as a String; ```-contain```/```-scan```
with ```-scaled``` never build the String, so a single chromosome of any size is fine.


### Where the element sits

```-contain``` says how much of a sequence is there, ```-scan``` says **where**. A window the size of
the query is walked along every other sequence, base by base, and each region that holds the query
is reported with its bounds, the identity of the match, and an E-value - the number of regions this
good that an unrelated sequence of that length would be expected to yield, so that the best of a
hundred thousand windows is not mistaken for a find.

On ```test/2.txt```, looking for the 18S rRNA of the macaque (the first sequence):

```
> java -jar GeneDistance.jar test/2.txt -scan

Query: 18S_Macaca_mulatta (1868 bp, 1856 distinct k-mers of length 11)
5 region(s) hold the query.
   18S_1kMacaca_mulatta  1..1868 (1868 bp), best window at 1  containment 100%  identity 100%  matches 1  E=0.0
   18S_2kMacaca_mulatta  1..1868 (1868 bp), best window at 1  containment 100%  identity 100%  matches 1  E=0.0
   18S-28S_Macaca_mulatta  3641..5819 (2179 bp), best window at 3641  containment 100%  identity 100%  matches 5  E=0.0
   18S-28S_Homo_sapiens  3590..5664 (2075 bp), best window at 3644  containment 99%  identity 100%  matches 10  E=0.0
   18S_Homo_sapiens  1..1869 (1869 bp), best window at 1  containment 99%  identity 100%  matches 3  E=0.0
```

The 18S is found at base 3641 of the 18S-28S unit, after the 5' ETS, and nothing is reported in the
random sequences, the 28S or the unrelated genes: their best windows hold 2-7 k-mers of the query
where chance gives 1.4-4.6, and score E = 16 to 730 (chance alone would give as much that many
times over). *matches* is the number of separate exact matches under the window: one for an
identical copy, more as substitutions cut it into pieces.

The report also holds the profile - the containment of the query under the window, along each
sequence - which plots as a curve of where the query is and is not.

*E-value*: the number of windows holding that many k-mers of the query that a sequence of no relation
would be expected to yield (Altschul SF et al. (1990) *J Mol Biol* 215:403-410). Two things are
taken from the sequence searched rather than assumed. The chance of meeting each k-mer of the query
follows from the **base composition of that sequence** - an AT-rich k-mer turns up in AT-rich DNA
orders of magnitude more often than 4<sup>-k</sup>. And the k-mers met by chance are **not
independent**: a chance match longer than k holds several consecutive k-mers, so their number under
a window follows a geometric Poisson (Polya-Aeppli) law, not a Poisson one. With both, none of 560
searches of unrelated chicken DNA with Arabidopsis repeats reports a region; treating the k-mers as
independent and equally likely, as the first versions did, reported one in 13% of them. The values
are conservative. Windows
overlap, so it remains the usual first approximation, not an exact probability. A printed 0.0 is a
value below the range of a double.

*containment*: Broder AZ (1997) *SEQUENCES'97*; Koslicki D, Zabeti H (2019) *Appl Math Comput*
354:206-215; Ondov BD et al. (2019) *Mash Screen*, *Genome Biology* 20:232. *Choice of k*: Ondov BD
et al. (2016) *Mash*, *Genome Biology* 17:132.


### Output:

Written next to the input (```result_kN.*``` when a folder is given):

| File | Content |
| --- | --- |
| ```<input>_kN.xls``` | ```-spacing```: k-mer counts, average k-mer spacings, and the similarity matrix (%) |
| ```<input>_kN.meg``` | distance matrix for MEGA (lower-left), to build a tree; for ```-spacing``` it is 1 − S |
| ```<input>_cK.xls``` | ```-contain```: containment matrix (%) and the identity (%) of the shared part |
| ```<input>_cK.meg``` | ```-contain```: distance matrix (1 - identity) for MEGA |
| ```<input>_scanK.xls``` | ```-scan```: best window per sequence, the regions that hold the query (position, identity, separate matches, E-value) and the containment profile |

A ```-cosine```, ```-d2star``` or ```-vector``` run writes ```_kN_cos.*```, ```_kN_d2s.*``` or
```_kN_vec.*```; a sketched ```-contain``` / ```-scan``` writes ```_cKsN.*``` / ```_scanKsN.*```. The
comparison measures therefore never overwrite one another, and they hold less than ```_kN.xls```
above: ```_kN_cos.xls``` and ```_kN_d2s.xls``` carry the k-mer counts and the similarity matrix,
```_kN_vec.xls``` the similarity matrix alone — the k-mer-distance block is written by ```-spacing```
only. The ```.meg``` matrix is the same in all of them. ```-kmerstat``` and ```-kmer2stat``` are the
exception: both write their table to ```<input>_kN.xls```, the name ```-spacing``` uses, so either
overwrites a previous ```-spacing``` report and leaves its ```_kN.meg``` beside it, now stale.

Exit codes: 0 = done, 1 = wrong usage or no sequence found, 2 = I/O error or out of memory. Under
```-kmerstat``` / ```-kmer2stat``` a failed read is reported and the run carries on, so an input that
cannot be read there ends at 1, or at 0 on a partial report when only some of the files failed.


### Examples:
```
java -jar C:\GeneDistance\dist\GeneDistance.jar C:\KmersPattern\test\t1.txt

java -jar C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Oryza_sativa\ 

java -jar C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Oryza_sativa\ -kmer=41

java -jar C:\GeneDistance\dist\GeneDistance.jar C:\KmersPattern\test\t1.txt -kmerstat

```
