from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

    try:
        # Convert sequences to proper format
        records = []
        for seq in sequences:
            if isinstance(seq.seq, str):
                seq.seq = Seq(seq.seq)
            records.append(seq)

        # Create alignment
        alignment = MultipleSeqAlignment(records)

        # Ensure all sequences are the same length
        max_length = max(len(seq) for seq in alignment)
        for record in alignment:
            if len(record.seq) < max_length:
                # Pad with gaps if necessary
                record.seq = record.seq + "-" * (max_length - len(record.seq))

        return alignment

    except Exception as e:
        raise ValueError(f"Error during sequence alignment: {str(e)}")

def calculate_distance_matrix(alignment):
    """
    Calculate distance matrix from alignment

    Args:
        alignment: MultipleSeqAlignment object

    Returns:
        numpy.ndarray: Distance matrix
    """
    try:
        calculator = DistanceCalculator('identity')
        return calculator.get_distance(alignment)
    except Exception as e:
        raise ValueError(f"Error calculating distance matrix: {str(e)}")