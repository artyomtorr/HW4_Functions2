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


