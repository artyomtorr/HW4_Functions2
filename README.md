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
- 

### Examples
```python
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'GSCKRGPRT', 'compute_length') # [10, 18, 9]
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'GSCKRGPRT', 'compute_molecular_weight') # [1055.496, 1886.872, 942.482]
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'GSCKRGPRT', 'compute_hydrophobicity') # [50.0, 27.778, 11.111]
run_protein_tools('AUGGAUCAUcAAUAA', 'MDKL*', 'check_mutations') #'Mutations:K3, L4.'
```
   
### Additional information
- The program works **only** with protein and RNA sequences. If any of the entered sequences contain inappropriate characters or cannot exist, the program will display an error. Sequences can contain characters of any case.

```python
run_protein_tools('PROTEIN', 'compute_molecular_weight') # ValueError: Invalid protein sequence
run_protein_tools('AUGGAU_AUcAAUAA', 'MDKL*', 'check_mutations')# ValueError: Invalid RNA sequence
```

### Contacts
Please use contacts below to reach out with any comments, concerns, or discussions regarding **protein_tools.py.** <br>
- Artyom Toropov ([@artyomtorr](github.com/artyomtorr)) <br>
- Sofiya Vinogradova ([@sofiyaga57](github.com/sofiyaga57)) <br>
- Nikita Zherko ([@rereremin](github.com/rereremin)) <br>
![изображение](https://github.com/artyomtorr/HW4_Functions2/assets/144557024/88f1c523-711a-40d7-9134-30c6b6639037)


*Author contributions:* <br> 
Artyom Toropov (teamlead): functions *is_protein*, *molecular_weight*, *run_protein_tools* <br> 
Sofiya Vinogradova: functions ..., <br> 
Nikita Zherko: functions *compute_hydrophobicity*, *check_mutations*.
