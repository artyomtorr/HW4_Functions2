# protein_tools.py

**protein_tools.py** - is a tool which allows the performing of various procedures for a user entered protein sequences. 

### Usage

The tool works by calling the function `run_protein_tools`, which takes arbitrary number of arguments with protein sequencies (*str*) and the name of the procedure to be performed (always the last argument, *str*, see the usage examples below). The output is the result of the procedure as *string, tuple* or *dictionary* if one sequence is submitted or *list* if several.

**NOTE:**  For the procedure `check_mutations` a fixed number of string arguments are used: one RNA sequence, one protein sequence and the name of procedure itself.

### Procedures

- `compute_molecular_weight` — computes molecular weight of protein sequence in g/mol
- `compute_length` — computes the number of amino acids in protein sequence
- `compute_hydrophobicity` — computes the percentage of gydrophobic aminoacids in protein sequence
- `check_mutations` — checks missense mutations in the protein sequence after translation
- `protein_to_dna`- returns possible variants of DNAs for a given protein sequence
- `count_amino_acids` - calculates the number of each aminoacid in protein sequence

### Examples
```python
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_length')
#[('MAEGEITNLP', 10), ('tGQYLAMDTSgLLYGSQT', 18)]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_molecular_weight')
#[('MAEGEITNLP', 1055.496), ('tGQYLAMDTSgLLYGSQT', 1886.872)]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_hydrophobicity')
#[('MAEGEITNLP', 50.0), ('tGQYLAMDTSgLLYGSQT', 27.778)]

run_protein_tools('AUGGAUCAUcAAUAA', 'MDKL*', 'check_mutations')
#'Mutations: K3, L4.'

run_protein_tools('MAEGLP', 'LYGSQT','protein_to_dna')
#['ATG GCT/GCC/GCA/GCG GAA/GAG GGT/GGC/GGA/GGG TTA/TTG/CTT/CTC/CTA/CTG CCT/CCC/CCA/CCG',
#'TTA/TTG/CTT/CTC/CTA/CTG TAT/TAC GGT/GGC/GGA/GGG TCT/TCC/TCA/TCG/AGT/AGC CAA/CAG ACT/ACC/ACA/ACG']

run_protein_tools('MAEGLP', 'LYGSQT','count_amino_acids')
#[{'M': 1, 'A': 1, 'E': 1, 'G': 1, 'L': 1, 'P': 1},
#{'L': 1, 'Y': 1, 'G': 1, 'S': 1, 'Q': 1, 'T': 1}]
```
   
### Additional information
- The program works **only** with protein and RNA sequences. If any of the entered sequences contain inappropriate characters or cannot exist, the program will display an error. Sequences can contain characters of any case.

```python
run_protein_tools('PROTEIN', 'compute_molecular_weight') #ValueError: Invalid protein sequence
run_protein_tools('AUGGAU_AUcAAUAA', 'MDKL*', 'check_mutations') #ValueError: Invalid RNA sequence
```
- For the procedure `check_mutations` there are extra requirements for RNA and protein sequences: mRNA sequences must contain **start-codon** and **one of the stop-codons**, protein sequnces must start with **"M"** and ends with **"*"** (stop-codon). 
```python
run_protein_tools("AUGGUAGGGAAAUUUUGA", "MGGKF", 'check_mutations') #ValueError: Stop (*) is absent
run_protein_tools("AUGGUAGGGAAAUUUUGA", "GGKF*", 'check_mutations') #ValueError: Start (M) is absent
```
### Contacts
Please use contacts below to reach out with any comments, concerns, or discussions regarding **protein_tools.py.** <br>
- Artyom Toropov ([@artyomtorr](https://github.com/artyomtorr/)) <br>
- Sofiya Vinogradova ([@sofiyaga57](https://github.com/sofiyaga57/)) <br>
- Nikita Zherko ([@rereremin](https://github.com/rereremin/)) <br>
![изображение](https://github.com/artyomtorr/HW4_Functions2/assets/144557024/88f1c523-711a-40d7-9134-30c6b6639037)


*Author contributions:* <br> 
Artyom Toropov (teamlead): functions `is_protein`, `is_rna`, `compute_molecular_weight`, `run_protein_tools` <br> 
Sofiya Vinogradova: functions `compute_length`, `count_amino_acids`, `protein_to_dna` <br> 
Nikita Zherko: functions `compute_hydrophobicity`, `translate_rna`, `check_mutations`
