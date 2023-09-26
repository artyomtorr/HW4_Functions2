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
