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


def is_protein(seq):
    unique_chars = set(seq.upper())
    return unique_chars <= alphabet_protein


def is_rna(seq):
    unique_chars = set(seq.upper())
    return unique_chars <= alphabet_rna


def compute_molecular_weight(seq):
    molecular_weight = 0
    for amino_acid in seq:
        molecular_weight += amino_acid_masses[amino_acid]
    return round(molecular_weight, 3)


def compute_length(seq: str):
    """
    Compute the length of the input amino acid sequence
    
    """
    return len(seq)    


def compute_hydrophobicity(protein):

    count_of_gydrophobic = 0
    if is_protein(protein):
        for i in range(len(protein)):
            if protein[i] in gydrophobic_aminoacids:
                count_of_gydrophobic += 1

    percentage = round(count_of_gydrophobic / len(protein) * 100, 3)

    return f"Percentage of gydrophobic aminoacids in {protein} = {percentage}%."


def translation(seq):
    """
    """
    triplets = [seq[i:i + 3].upper() for i in range(0, len(seq), 3)]
    protein = []
    for triplet in triplets:
        for aminoacid in rna_codons.keys():
            if triplet in rna_codons[aminoacid]:
                protein.append(aminoacid)

    start = protein.index("M")
    stop = protein.index("*")
    return "".join(protein[start:stop + 1])


def check_mutations(seq, protein):

    if is_protein(protein[:-1]) is not True:
        raise ValueError("Invalid protein sequence")
    if is_rna(seq) is not True:
        raise ValueError("Invalid RNA sequence")

    correct_protein = translation(seq)
    bank_of_mutations = []
    
    for i in range(len(correct_protein)):
        if correct_protein[i] != protein[i]:
            bank_of_mutations.append(f'{protein[i]}{i + 1}')

    if len(bank_of_mutations) == 0:
        return "Protein without mutations."
    else:
        return "Mutations:" + ", ".join(bank_of_mutations) + "."


def run_protein_tools(*seqs_and_procedure):
    procedure = seqs_and_procedure[-1]
    seqs = seqs_and_procedure[:-1]

    results = []
    if procedure == 'check_mutations':
        results.append(check_mutations(seqs[0], seqs[1]))
        
    else:
        for seq in seqs:
            seq = seq.upper()
            if is_protein(seq) is not True:
                raise ValueError("Invalid protein sequence")
            if procedure == 'compute_molecular_weight':
                results.append(molecular_weight(seq))
            elif procedure == 'compute_length':
                results.append(compute_length(seq))
            elif procedure == 'compute_hydrophobicity':
                results.append(compute_hydrophobicity(seq))
    if len(results) == 1:
        return results[0]
    else:
        return results
      