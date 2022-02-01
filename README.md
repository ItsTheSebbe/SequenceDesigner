# Sequence Designer for caDNAno

## Overview
This program makes it possible to sequence scaffold and staple strands, even if the cadnano file contains multiple separate scaffolds.

## Usage
The program can be run using the following command:

```python
python3 seq_designer.py <cadnano json file> <scaffold file>
```
For example:
```python
python3 seq_designer.py json_files/test_virtual.json scaffold_files/M13mp18 
```
## Input
The program will require two inputs as arguments:
- cadnano .json file
- scaffold sequence file

## Output
The program will generate three output files:
- scaffolds.txt
- staples.txt
- visualized_sequence.txt
### scaffolds.txt
scaffolds.txt contains the sequences of the scaffold strands. Moreover, it contains the start and end location, and the length of each scaffold.
### staples.txt
staples.txt contains the sequences of the staple strands. Moreover, it contains the start and end location, and the length of each staple.
### visualized_sequence.txt
visualized_sequence.txt contains a nicely formatted visualization of the scaffold and staple sequence data, analogous to the visual representation in cadnano. This might be useful for checking the final results.

## Example Output
Here are some examples for the outputs using a small caDNAno file. (json_files/small_twobreak.json and M13mp18 scaffold). 
![image](https://user-images.githubusercontent.com/28595211/152023399-08882dd2-fa9e-4b99-ab19-b349729acefc.png)


### scaffolds.txt
```
Start,End,Sequence,Length
1[6],0[5],GTGATGATT,9
0[6],1[7],AATGCTACTAC,11
```
### staples.txt
```
Start,End,Sequence,Length
1[2],1[11],ATCACGTAGT,10
0[11],0[2],AGCATTAATC,10
```
### visualized_sequence.txt
```
Scaffold 0    |--GATTAATGCT---------|
Staple 0      |--CTAATTACGA---------|

Staple 1      |--ATCACGTAGT---------|
Scaffold 1    |--TAGTGCATCA---------|
```
