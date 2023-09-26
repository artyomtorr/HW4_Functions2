
def level_of_hydrophobic(protein):
    gydrophobic_aminoacids = {"A", "V", "L", "I", "P", "F", "W", "M"}

    count_of_gydrophobic = 0
    if is_protein(protein):
        for i in range(len(protein)):
            if protein[i] in gydrophobic_aminoacids:
                count_of_gydrophobic += 1

    percentage = count_of_gydrophobic / len(protein) * 100

    return f"Percentage of gydrophobic aminoacids in {protein} = {percentage}%."


def translation(seq):
    """
    """
    gene_code = {
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
    triplets = [seq[i:i + 3].upper() for i in range(0, len(seq), 3)]
    protein = []
    for triplet in triplets:
        for aminoacid in gene_code.keys():
            if triplet in gene_code[aminoacid]:
                protein.append(aminoacid)

    if is_protein("".join(protein)):
        start = protein.index("M")
        stop = protein.index("*")
        return "".join(protein[start:stop + 1])
    else:
        return "This sequence doesn't include the gene."


def mutations(seq, protein):
    correct_protein = translation(seq)

    if is_protein(protein):
        bank_of_mutations = []
        for i in range(len(correct_protein)):
            if correct_protein[i] != protein[i]:
                bank_of_mutations.append(f'{protein[i]}{i + 1}')

        if len(bank_of_mutations) == 0:
            return "Protein without mutations."
        else:
            return "Mutations:" + ", ".join(bank_of_mutations) + "."
    else:
        return "It isn't a protein."


alphabet_protein = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
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


def is_protein(seq):
    unique_chars = set(seq)
    return unique_chars <= alphabet_protein


def molecular_weight(seq):
    molecular_weight = 0
    for amino_acid in seq:
        molecular_weight += amino_acid_masses[amino_acid]
    return round(molecular_weight, 3)


def run_protein_tools(*seqs_and_procedure):
    procedure = seqs_and_procedure[-1]
    seqs = seqs_and_procedure[:-1]

    results = []

    for seq in seqs:
        seq = seq.upper()
        if is_protein(seq) is not True:
            raise ValueError("Invalid alphabet")
        if procedure == 'molecular_weight':
            results.append(molecular_weight(seq))

    if len(results) == 1:
        return results[0]
    else:
        return results

