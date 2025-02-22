from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from io import StringIO
import numpy as np

def process_fasta(file_content):
    """
    Process FASTA file content and return sequences
    
    Args:
        file_content: StringIO object containing FASTA data
    
    Returns:
        list: List of SeqRecord objects
    """
    try:
        sequences = list(SeqIO.parse(file_content, "fasta"))
        if not sequences:
            raise ValueError("No valid sequences found in the file")
        return sequences
    except Exception as e:
        raise ValueError(f"Invalid FASTA format: {str(e)}")

def perform_alignment(sequences, method="Multiple Sequence Alignment"):
    """
    Perform sequence alignment using specified method
    
    Args:
        sequences: List of SeqRecord objects
        method: String indicating alignment method
    
    Returns:
        MultipleSeqAlignment: Aligned sequences
    """
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required for alignment")

    if method == "Multiple Sequence Alignment":
        # Create a simple progressive alignment
        alignment = MultipleSeqAlignment(sequences)
        return alignment
    else:
        # Perform pairwise alignment
        from Bio import pairwise2
        alignments = []
        reference = sequences[0]
        
        for seq in sequences[1:]:
            alignment = pairwise2.align.globalxx(reference.seq, seq.seq)[0]
            alignments.append(alignment)
            
        return MultipleSeqAlignment(alignments)

def calculate_distance_matrix(alignment):
    """
    Calculate distance matrix from alignment
    
    Args:
        alignment: MultipleSeqAlignment object
    
    Returns:
        numpy.ndarray: Distance matrix
    """
    calculator = DistanceCalculator('identity')
    return calculator.get_distance(alignment)
