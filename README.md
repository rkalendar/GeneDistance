## GeneDistance
## Identifications of genetic similarity or the distance between genomic sequences (approximate matching algorithm)

Alignment-free method for calculating genetic distances between DNA sequences as a basis for similarity distance estimation and phylogenetic reconstruction. The proposed algorithm is effective for measuring homology in cases of mixed sequences and different lengths, including individual chromosomes of the same or different species can be used for analysis. The application of the proposed algorithm is not limited to determining the degree of homology between sequences, but can be used for phylogenetic analysis, species identification and sequence classification in genomic assemblies.

## Author
Ruslan Kalendar 
email: ruslan.kalendar@helsinki.fi

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 23 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)

How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)

The program generates a file for analysis in the software MEGA 12: https://www.megasoftware.net/

To run the project from the command line (CLI). Command-line options, separated by spaces. 
The executive file ```GeneDistance.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar GeneDistance.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <GeneDistancePath>\dist\GeneDistance.jar <target_file_path> optional_commands```


### Examples:
```
java -jar C:\GeneDistance\dist\GeneDistance.jar C:\GeneDistance\test\t1.txt

java -jar C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Chloroplast\ -mod=2

```

### Large genome usage (you will have to show the program to use more RAM, for example as listed here, up to 64 Gb memory: -Xms16g -Xmx64g):
```
java -jar -Xms16g -Xmx64g C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\T2T-CHM13v2.0\ -mod=3
```
For chromosomes larger than 500 Mb you will need to use more memory, 128 Gb:
```
java -jar -Xms32g -Xmx128g C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Cycas_panzhihuaensis\ -mod=4
```


## Sequence Entry:
Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of the following:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.





