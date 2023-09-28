# protein_tools.py

**protein_tools.py** - is a tool which allows the performing of various procedures for a user entered protein sequences. 

### Usage

The tool works by calling the function `run_protein_tools`, which takes arbitrary number of arguments with protein sequencies (*str*) and the name of the procedure to be performed (always the last argument, *str*, see the usage examples below). The output is the result of the procedure as *string* if one sequence is submitted or *list* if several.

**NOTE:**  For the procedure `check_mutations` a fixed number of string arguments are used: one RNA sequence, one protein sequence and the name of procedure itself.

### Procedures

- `compute_molecular_weight` — computes molecular weight of protein sequence in g/mol
- `compute_length` — computes the number of amino acids in protein sequence
- `compute_hydrophobicity` — computes the percentage of gydrophobic aminoacids in protein sequence
- `check_mutations` — 

### Examples
```python
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'GSCKRGPRT', 'compute_length') # [10, 18, 9]
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'GSCKRGPRT', 'compute_molecular_weight') # [1055.496, 1886.872, 942.482]
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'GSCKRGPRT', 'compute_hydrophobicity') # [50.0, 27.778, 11.111]

```
   
### Additional information
- The program works **only** with protein or RNA sequences. If any of the entered sequences contain inappropriate characters or cannot exist, the program will display an error. Sequences can contain characters of any case.

```python
run_protein_tools('ATA', 'DefinitelyNotDNA', 'transcribe') # ValueError: Invalid alpabet
run_protein_tools('ATGU', 'reverse') # ValueError: Invalid alpabet
```

### Contacts
Author contributions:
