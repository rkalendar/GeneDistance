## GeneDistance
## Identifications of genetic similarity or the distance between genomic sequences (in approximate string matching)

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

| Option | Meaning |
| --- | --- |
| ```-kmer=N``` | set of k-mers used in the analysis, ```N``` = 4/5/6/8/10 (symmetric) or 41/61/81/101 (non-symmetric), default ```-kmer=4``` |
| ```-d2star``` | compare the k-mer **frequencies** instead of their spacing, each count centred on the base composition of its own sequence (the *d2\** measure) |
| ```-cosine``` | the same, on the raw frequencies, without centring |
| ```-vector``` | *d2\** with **scale-matched windows and both strands**: a short fragment is compared to the best window of a longer one, on either strand, and scores no better than chance (measured from the input per length scale) are dropped. For homology between fragments of very different length |
| ```-contain``` | **containment**: how much of each sequence occurs in each other one, on the whole k-mer space. Finds a short element inside a long sequence |
| ```-scan```, ```-scan=N``` | **where** each sequence holds the query: a window the size of the query is walked along them. The query is the first sequence of the input, or the N-th |
| ```-ksize=K``` | the k of the exact k-mers ```-contain``` and ```-scan``` use (6-31). Chosen for the length of the sequences unless given |
| ```-kmerstat``` | report the k-mer distances of each sequence, instead of comparing them |
| ```-kmer2stat``` | report the k-mer distances averaged over all the target sequences |
| ```-h```, ```--help``` | print the help and exit |
| ```--version``` | print the version and exit |


### Sequences of unequal length

The default measure compares how the k-mers are *spaced* inside each sequence, and a k-mer only
counts once it has been seen at least twice. A short sequence therefore arrives with an almost
empty vector, and its homology to a long one goes unseen. Use ```-d2star``` there. On
```test/2.txt``` (18S rRNA of macaque and man, the same 18S diluted with random DNA, and three
random sequences), with ```-kmer=41```:

| pair | truth | default | ```-cosine``` | ```-d2star``` |
| --- | --- | --- | --- | --- |
| 18S macaque / 18S human | orthologs | 99 | 100 | 100 |
| 18S / the 18S-28S unit holding it | homologous | 55 | 83 | 51 |
| 18S / 18S + 2 kb of random DNA | homologous | 77 | 96 | 75 |
| 18S / random DNA | **unrelated** | **61** | **87** | **6** |
| random / random | **unrelated** | **65** | **89** | **2** |

The reports of the measures are written side by side, as ```_kN.*```, ```_kN_cos.*```,
```_kN_d2s.*```, ```_kN_vec.*```, ```_cK.*``` and ```_scanK.xls```, and never overwrite one another.

*d2\**: Reinert G, Chew D, Sun F, Waterman MS (2009) *J Comput Biol* 16:1615-1634; Song K et al.
(2014) *Brief Bioinform* 15:343-353.


### Fragments of very different length: ```-vector```

```-d2star``` compares whole sequences, which averages a shared region away when one sequence is
much longer than the other. ```-vector``` fixes that the way TotalRepeats' *-vector* mode does: a
short sequence is compared not to a long one as a whole but to the **best window of the long one at
the short one's own scale** (step = half a window), on **either strand**, and a score is kept only
when it stands clear of the **chance level measured for that length scale** from the input's own
content (non-homologous windows screened by an exact 25-mer). On ```test/2.txt``` with ```-kmer=41```:

| pair | truth | ```-d2star``` | ```-vector``` |
| --- | --- | --- | --- |
| 18S / the 18S-28S unit holding it | 18S is inside it | 51 | **97** |
| 28S / the unit (reverse-complemented inside) | 28S is inside it | 45 | **98** |
| 18S / 18S + 8.6 kb random | homologous | 29 | 100 |
| any / random, barley, arabidopsis | **unrelated** | small | **0** |

The calibrated floor needs enough content to measure the chance level; on a small input a length
scale that cannot be calibrated is simply left uncut (its raw, already well-separated scores are
reported). For a short element that shares only part of *both* sequences (e.g. two different genomic
regions each carrying the same element plus their own flanks), ```-contain``` remains the sharper
tool — it asks "how much is shared" rather than "is the composition matched at scale".

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
whole genomes wants ```java -Xmx8g -jar ...```.


### Where the element sits

```-contain``` says how much of a sequence is there, ```-scan``` says **where**. A window the size of
the query is walked along every other sequence, base by base, and each region that holds the query
is reported with its bounds, the identity of the match, and an E-value - the number of regions this
good that an unrelated sequence of that length would be expected to yield, so that the best of a
hundred thousand windows is not mistaken for a find.

An Alu element planted, 5% mutated, at base 5000 of 10 kb of unrelated DNA:

```
> java -jar GeneDistance.jar mix.fa -scan

Query: Q_AluYb11_289 (289 bp, 279 distinct k-mers of length 11)
4 region(s) hold the query.
   M_AluYb11_5pct           1..289 (289 bp), best window at 1     containment 70%  identity 97%  E=0.0
   H_AluJo_283_real_homolog 1..274 (274 bp), best window at 1     containment 39%  identity 92%  E=0.0
   L_10kb_CONTAINS_AluYb11  5008..5288 (281 bp), best window at 5000  containment 75%  identity 97%  E=0.0
```

The element is found where it was put, and the 10 kb of unrelated DNA holds nothing (its best window
scores E = 180, which is to say chance alone would give as much 180 times over). On ```test/2.txt```
the 18S of the macaque is found at base 3641 of the 18S-28S unit, after the 5' ETS, and the 28S at
7900, reverse complemented.

The report also holds the profile - the containment of the query under the window, along each
sequence - which plots as a curve of where the query is and is not.

*E-value*: the Poisson tail of the count of shared k-mers over the windows walked, in the manner of
Altschul SF et al. (1990) *J Mol Biol* 215:403-410. The k-mers of a window overlap, so it is the
usual first approximation, not an exact probability.

*containment*: Broder AZ (1997) *SEQUENCES'97*; Koslicki D, Zabeti H (2019) *Appl Math Comput*
354:206-215; Ondov BD et al. (2019) *Mash Screen*, *Genome Biology* 20:232. *Choice of k*: Ondov BD
et al. (2016) *Mash*, *Genome Biology* 17:132.


### Output:

Written next to the input (```result_kN.*``` when a folder is given):

| File | Content |
| --- | --- |
| ```<input>_kN.xls``` | k-mer counts, k-mer distances, and the similarity matrix (%) |
| ```<input>_kN.meg``` | distance matrix for MEGA (lower-left), to build a tree |
| ```<input>_cK.xls``` | ```-contain```: containment matrix (%) and the identity (%) of the shared part |
| ```<input>_cK.meg``` | ```-contain```: distance matrix (1 - identity) for MEGA |
| ```<input>_scanK.xls``` | ```-scan```: best window per sequence, the regions that hold the query (position, identity, E-value) and the containment profile |

Exit codes: 0 = done, 1 = wrong usage or no sequence found, 2 = I/O error.


### Examples:
```
java -jar C:\GeneDistance\dist\GeneDistance.jar C:\KmersPattern\test\t1.txt

java -jar C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Oryza_sativa\ 

java -jar C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Oryza_sativa\ -kmer=41

java -jar C:\GeneDistance\dist\GeneDistance.jar C:\KmersPattern\test\t1.txt -kmerstat

```
