alphabet_protein = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

alphabet_rna = {'A', 'U', 'G', 'C'}

amino_acid_masses = {
    'A': 71.03711,
    'R': 156.10111,
    'N': 114.04293,
    'D': 115.02694,
    'C': 103.00919,
    'Q': 128.05858,
    'E': 129.04259,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'L': 113.08406,
    'K': 128.09496,
    'M': 131.04049,
    'F': 147.06841,
    'P': 97.05276,
    'S': 87.03203,
    'T': 101.04768,
    'W': 186.07931,
    'Y': 163.06333,
    'V': 99.06841
}

gydrophobic_aminoacids = {"A", "V", "L", "I", "P", "F", "W", "M"}

dna_codons = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'C': ['TGT', 'TGC'],
        'D': ['GAT', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['TTT', 'TTC'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'K': ['AAA', 'AAG'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'M': ['ATG'],
        'N': ['AAT', 'AAC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        '*': ["UAA", "UAG", "UGA"]}

rna_codons = {
        "F": ["UUC", "UUU"], "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "I": ["AUU", "AUC", "AUA"], "M": ["AUG"], "V": ["GUU", "GUC", "GUA", "GUG"],
        "S": ["UCU", "UCC", "UCA", "UCG"], "P": ["CCU", "CCC", "CCA", "CCG"],
        "T": ["ACU", "ACC", "ACA", "ACG"], "A": ["GCU", "GCC", "GCA", "GCG"],
        "Y": ["UAC", "UAU"], "*": ["UAA", "UAG", "UGA"], "H": ["CAU", "CAC"],
        "Q": ["CAA", "CAG"], "N": ["AAU", "AAC"],
        "K": ["AAA", "AAG"], "D": ["GAU", "GAC"], "E": ["GAA", "GAG"],
        "C": ["UGU", "UGC"], "W": ["UGG"], "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "S": ["AGU", "AGC"], "G": ["GGU", "GGC", "GGA", "GGG"]
    }


def is_protein(seq:str):
    """
    Check the existence of a protein sequence, return boolean.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= alphabet_protein


def is_rna(seq:str):
    """
    Check the existence of a RNA sequence, return boolean.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= alphabet_rna


def compute_molecular_weight(seq:str):
    """
    Compute molecular weight (g/mol) of protein sequence.
    """
    molecular_weight = 0
    for amino_acid in seq.upper():
        molecular_weight += amino_acid_masses[amino_acid]
    return round(molecular_weight, 3)


def compute_length(seq:str):
    """
    Compute the length of protein sequence.
    """
    return len(seq)    


def compute_hydrophobicity(protein:str) -> tuple:
    """
    Compute the percentage of gydrophobic aminoacids in protein sequence.

    Argument:
    - protein (str): protein sequence. Include hydrophobic 
    and hydrophilic aminoacids.

    Return:
    - tuple, result of computation percentage of gydrophobic aminoacids.
    """
    count_of_gydrophobic = 0
    for i in range(len(protein)):
        if protein[i] in gydrophobic_aminoacids:
            count_of_gydrophobic += 1

    percentage = round(count_of_gydrophobic / len(protein) * 100, 3)

    return protein, percentage


def translate_rna(seq:str) -> str:
    """
    Perform the translation of mRNA seguence into protein sequence.

    Argument:
    - seq (str): mRNA sequence. Must contain start-codon and one of 
    the stop-codons.

    Return:
    - str, protein sequence after translation. 
    Always starts with "M" and ends with "*".
    """
    triplets = [seq[i:i + 3].upper() for i in range(0, len(seq), 3)]
    protein = []
    for triplet in triplets:
        for aminoacid in rna_codons.keys():
            if triplet in rna_codons[aminoacid]:
                protein.append(aminoacid)

    if protein[-1] != "*":
        raise ValueError("Stop-codon (*) is absent in mRNA")
    if protein[0] != "M":
        raise ValueError("Start-codon (M) is absent in mRNA")

    start = protein.index("M")
    stop = protein.index("*")
    return "".join(protein[start:stop + 1])


def check_mutations(seq:str, protein:str) -> str:
    """
    Check mutations in the protein sequence after translation.

    Use additional function "translation(seq)".
    This function doesn't show mutations, which don't lead to 
    change aminoacids in protein sequence. 

    Arguments:
    - seq (str): translation sequence of mRNA with/without mutations
    - protein (str): protein for comparison with protein after translation.
    Every protein starts with "M" and ends with "*" (stop-codon). 
    Remark: is_protein(seq) doesn't see "*", but it's used in the other part of function.

    Return:
    - str, if mRNA without mutations return "Protein without mutations." 
    If some mutations in protein, return aminoacid(s) and their position(s)

    Examples:
    - "AUGGUAGGGAAAUUUUGA", "MVGKF*" ->  "Protein without mutations."
    - "AUGGUAGGGAAAUUUUGA", "MGGVF*" -> "Mutations:G2, V4."
    - "AUGGUAGGGAAAUUUUGA", "MGGKF" –> "ValueError: Stop (*) is absent"
    - "AUGGUAGGGAAAUUUUGA", "GGKF*" –> "ValueError: Start (M) is absent"
    - "AUGAAAAAAUGA", "MK*" -> "ValueError: Different length of translated protein and protein"
    """

    correct_protein = translation(seq)
    bank_of_mutations = []
    
    if is_protein(protein[:-1]) is not True:
        raise ValueError("Invalid protein sequence")
    if is_rna(seq) is not True:
        raise ValueError("Invalid RNA sequence")
    if protein[-1] != "*":
        raise ValueError("Stop (*) is absent")
    if protein[0] != "M":
        raise ValueError("Start (M) is absent")
    if len(protein) != len(seq)/3:
        raise ValueError("Different length of translated protein and protein")
    
    for i in range(len(correct_protein)):
        if correct_protein[i] != protein[i]:
            bank_of_mutations.append(f'{protein[i]}{i + 1}')

    if len(bank_of_mutations) == 0:
        return "Protein without mutations."
    else:
        return "Mutations:" + ", ".join(bank_of_mutations) + "."


def run_protein_tools(*args:str):
    """
    Function containing methods for protein analysis.
    
    Takes arbitrary number of arguments with protein sequencies
    and the name of the procedure to be performed (always the 
    last argument). Returns the result of the procedure as string 
    if one sequnce is submitted or list if several.

    If procedure 'check_mutations' is used then input must be only three
    arguments: RNA sequence, protein sequence and the name of procedure 
    itself.
    """
    *seqs, procedure = args
    results = []
    d_of_functions = {'compute_molecular_weight': compute_molecular_weight, 
                  'compute_length': compute_length,
                  'compute_hydrophobicity': compute_hydrophobicity,
                 }
    if procedure == 'check_mutations':
        results.append(check_mutations(seqs[0], seqs[1]))
    else:
        for seq in seqs:
            if is_protein(seq) is not True:
                raise ValueError("Invalid protein sequence")
            if procedure not in d_of_functions:
                raise ValueError("Wrong procedure name")
            else:
                results.append(d_of_functions[procedure](seq))
    if len(results) == 1:
        return results[0]
    else:
        return results
        
