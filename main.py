import streamlit as st
import pandas as pd
from io import StringIO
from utils.sequence_analysis import process_fasta, perform_alignment
from utils.visualization import create_phylogenetic_tree
from utils.sample_data import get_sample_data

# Page configuration
st.set_page_config(
    page_title="Pathogen Evolution Explorer",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS
with open('assets/app_style.css') as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)

def main():
    st.title("🧬 Pathogen Evolution Explorer")
    
    # Add detailed instructions
    st.markdown("""
    ### How to Test the Application:

    1. **Using Sample Dataset:**
       - Select "Use sample dataset" from the sidebar
       - Choose "SARS-CoV-2 Variants" from the dropdown menu
       - The phylogenetic tree will appear automatically

    2. **Uploading Your Own Data:**
       - Select "Upload FASTA file" from the sidebar
       - Click "Browse files" to upload your FASTA format file
       - The tree will generate after successful upload

    3. **Visualization Options:**
       - Choose between "Circular Tree" or "Rectangular Tree" layouts
       - Select your preferred alignment method
    """)

    st.markdown("""
    Explore and analyze evolutionary relationships of pathogens using genomic data.
    Upload your FASTA files or use our sample datasets to generate phylogenetic trees.
    """)

    # Sidebar
    st.sidebar.header("Controls")
    
    # Data input method selection
    input_method = st.sidebar.radio(
        "Choose input method:",
        ["Upload FASTA file", "Use sample dataset"]
    )

    sequences = None
    
    if input_method == "Upload FASTA file":
        uploaded_file = st.sidebar.file_uploader(
            "Upload FASTA file",
            type=["fasta", "fa", "fna"],
            help="Upload a FASTA format file containing genomic sequences"
        )
        
        if uploaded_file:
            try:
                content = StringIO(uploaded_file.getvalue().decode())
                sequences = process_fasta(content)
                st.success("File successfully uploaded and processed!")
            except Exception as e:
                st.error(f"Error processing file: {str(e)}")
                return
    else:
        sample_options = ["SARS-CoV-2 Variants", "Influenza H1N1", "E. coli Strains"]
        selected_sample = st.sidebar.selectbox("Select sample dataset:", sample_options)
        sequences = get_sample_data(selected_sample)

    if sequences:
        # Analysis options
        st.sidebar.header("Analysis Settings")
        alignment_method = st.sidebar.selectbox(
            "Alignment Method",
            ["Multiple Sequence Alignment", "Pairwise Alignment"]
        )
        
        visualization_type = st.sidebar.selectbox(
            "Visualization Type",
            ["Circular Tree", "Rectangular Tree"]
        )

        # Main content area
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.header("Phylogenetic Tree")
            with st.spinner("Generating phylogenetic tree..."):
                try:
                    alignment = perform_alignment(sequences, method=alignment_method)
                    fig = create_phylogenetic_tree(
                        alignment,
                        tree_type=visualization_type
                    )
                    st.plotly_chart(fig, use_container_width=True)
                except Exception as e:
                    st.error(f"Error generating tree: {str(e)}")

        with col2:
            st.header("Information Panel")
            st.markdown("""
            ### About Phylogenetic Trees
            
            Phylogenetic trees show evolutionary relationships between biological species
            or other entities that share a common ancestor.
            
            **Key Components:**
            - **Nodes**: Represent species or sequences
            - **Branches**: Show evolutionary relationships
            - **Branch Length**: Indicates genetic distance
            
            ### Analysis Methods
            
            **Multiple Sequence Alignment (MSA)**
            - Compares multiple sequences simultaneously
            - Identifies conserved regions
            - Shows evolutionary patterns
            
            **Pairwise Alignment**
            - Compares sequences two at a time
            - Useful for closely related sequences
            - Faster than MSA
            """)

    # Footer
    st.markdown("---")
    st.markdown(
        "Created with ❤️ using Streamlit, Biopython, and Plotly. "
        "For educational and research purposes."
    )

if __name__ == "__main__":
    main()