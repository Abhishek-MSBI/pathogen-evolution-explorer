<<<<<<< HEAD
---
title: Pathogen Evolution Explorer
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: streamlit
app_file: main.py
pinned: false
---

# 🧬 Pathogen Evolution Explorer

[![Hugging Face Spaces](https://img.shields.io/badge/Hugging%20Face-Spaces-blue)](https://huggingface.co/spaces/srabhishek/Pathogenevolutionexplorer)

A web-based application for analyzing pathogen genomic sequences and exploring their evolutionary relationships. Built with Streamlit and Biopython.

## Features

- **Sequence Statistics**: Calculate and visualize basic sequence metrics including GC content and length statistics
- **Sequence Comparison**: Generate similarity matrices and interactive heatmaps to compare multiple sequences
- **Sample Datasets**: Includes pre-loaded datasets for SARS-CoV-2 variants, Influenza H1N1, and E. coli strains
- **Phylogenetic Analysis**: *(Under Development)* Interactive phylogenetic tree visualization feature coming soon

## Quick Start

### Option 1: Using GitHub Desktop
1. In GitHub Desktop:
   - Click File > Clone Repository
   - Enter URL: `https://github.com/Abhishek-MSBI/pathogen-evolution-explorer`
   - Choose your local path
   - Click Clone

2. Install the required packages:
```bash
pip install streamlit biopython plotly pandas numpy


```

3. Run the application:
```bash
streamlit run main.py
```

### Option 2: Manual Download and Setup
1. Download the project files
2. Install the required packages:
```bash
pip install streamlit biopython plotly pandas numpy
```

3. Run the application:
```bash
streamlit run main.py
```

## Usage

1. **Using Sample Datasets**:
   - Select "Use sample dataset" from the sidebar
   - Choose from available sample datasets
   - Explore sequence statistics and comparisons in different tabs

2. **Analyzing Your Own Data**:
   - Select "Upload FASTA file" from the sidebar
   - Upload your FASTA format file
   - View analysis results across available tabs

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

## Future Features

- Enhanced phylogenetic tree visualization
- Integration with NCBI/EMBL APIs
- Advanced comparative analysis tools
- Export functionality for analysis results

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Pushing to GitHub

If using command line:
```bash
git init
git add .
git commit -m "Initial commit: Pathogen Evolution Explorer application"
git branch -M main
git remote add origin https://github.com/Abhishek-MSBI/pathogen-evolution-explorer.git
git push -u origin main
```

If using GitHub Desktop:
1. Open GitHub Desktop
2. Add Local Repository (File > Add Local Repository)
3. Select the project folder
4. Make initial commit
5. Publish repository
=======
---
title: Pathogenevolutionexplorer
emoji: ⚡
colorFrom: indigo
colorTo: pink
sdk: streamlit
sdk_version: 1.42.2
app_file: app.py
pinned: false
short_description: A web-based application for analyzing pathogen genomic seque
---

Check out the configuration reference at https://huggingface.co/docs/hub/spaces-config-reference
>>>>>>> space/main
