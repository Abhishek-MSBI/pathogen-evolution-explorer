from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_sample_data(dataset_name):
    """
    Provide sample genomic data for demonstration
    
    Args:
        dataset_name: String indicating which sample dataset to return
    
    Returns:
        list: List of SeqRecord objects
    """
    sample_data = {
        "SARS-CoV-2 Variants": [
            ("Wuhan-Hu-1", "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAA"),
            ("Alpha", "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAG"),
            ("Beta", "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTAC"),
            ("Delta", "ATGTTTGTTTTTCTTGTTTTATTGCCACTAGTCTCTAGTCAGTGTGTTGA")
        ],
        "Influenza H1N1": [
            ("H1N1-2009", "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATG"),
            ("H1N1-2018", "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATA"),
            ("H1N1-2019", "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATC"),
            ("H1N1-2020", "ATGAAGGCAATACTAGTAGTTCTGCTATATACATTTGCAACCGCAAATT")
        ],
        "E. coli Strains": [
            ("K-12", "ATGAGTCTTAATCATGCCCCTTGGCTTGTTGTTTTACGCTTAAAATCGAT"),
            ("O157:H7", "ATGAGTCTTAATCATGCCCCTTGGCTTGTTGTTTTACGCTTAAAATCGAA"),
            ("B", "ATGAGTCTTAATCATGCCCCTTGGCTTGTTGTTTTACGCTTAAAATCGAC"),
            ("W", "ATGAGTCTTAATCATGCCCCTTGGCTTGTTGTTTTACGCTTAAAATCGAG")
        ]
    }

    if dataset_name not in sample_data:
        raise ValueError(f"Unknown dataset: {dataset_name}")

    sequences = []
    for name, seq in sample_data[dataset_name]:
        sequences.append(
            SeqRecord(
                Seq(seq),
                id=name,
                name=name,
                description=f"Sample sequence for {name}"
            )
        )
    
    return sequences
