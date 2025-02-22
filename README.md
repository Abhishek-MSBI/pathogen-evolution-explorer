# 🧬 Pathogen Evolution Explorer

A web-based application for analyzing pathogen genomic sequences and exploring their evolutionary relationships. Built with Streamlit and Biopython.

## Features

- **Sequence Statistics**: Calculate and visualize basic sequence metrics including GC content and length statistics
- **Sequence Comparison**: Generate similarity matrices and interactive heatmaps to compare multiple sequences
- **Phylogenetic Analysis**: Create interactive phylogenetic trees to visualize evolutionary relationships
- **Sample Datasets**: Includes pre-loaded datasets for SARS-CoV-2 variants, Influenza H1N1, and E. coli strains

## Installation

1. Clone the repository:
```bash
git clone https://github.com/your-username/pathogen-evolution-explorer.git
cd pathogen-evolution-explorer
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Run the application:
```bash
streamlit run main.py
```

## Usage

1. **Using Sample Datasets**:
   - Select "Use sample dataset" from the sidebar
   - Choose from available sample datasets
   - Explore the analysis in different tabs

2. **Analyzing Your Own Data**:
   - Select "Upload FASTA file" from the sidebar
   - Upload your FASTA format file
   - View analysis results across all tabs

## Requirements

- Python 3.11+
- Streamlit
- Biopython
- Plotly
- Pandas
- NumPy

## Project Structure

```
├── .streamlit/          # Streamlit configuration
├── assets/             # Static assets and styles
├── utils/              # Utility functions
│   ├── sample_data.py  # Sample dataset provider
│   ├── sequence_analysis.py  # Sequence analysis tools
│   ├── sequence_stats.py    # Statistical calculations
│   └── visualization.py     # Visualization functions
└── main.py            # Main application file
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
