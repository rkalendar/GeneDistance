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


To run the project from the command line (CLI). Command-line options, separated by spaces. 
The executive file ```GeneDistance.jar``` is in the ```dist``` directory, which can be copied to any location. 
Go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar GeneDistance.jar <target_file_path/Folder_path>```


### Basic usage:

```java -jar <GeneDistancePath>\dist\GeneDistance.jar <target_file_path> optional_commands```


### Examples:
```
java -jar C:\GeneDistance\dist\GeneDistance.jar C:\GeneDistancePath\test\t1.txt

java -jar C:\GeneDistance\dist\GeneDistance.jar E:\Genomes\Oryza_sativa\ 

```
