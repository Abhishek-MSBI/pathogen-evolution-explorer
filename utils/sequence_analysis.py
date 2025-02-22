from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio import Align
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

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'

    if method == "Multiple Sequence Alignment":
        # Create a progressive alignment
        alignment = MultipleSeqAlignment([])
        # Add first sequence
        alignment.append(sequences[0])

        # Progressively align remaining sequences
        for seq in sequences[1:]:
            alignment.append(seq)

        return alignment
    else:
        # Perform pairwise alignment
        alignments = []
        reference = sequences[0]

        for seq in sequences[1:]:
            alignment = aligner.align(reference.seq, seq.seq)[0]
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