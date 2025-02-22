import streamlit as st
import pandas as pd
from io import StringIO
from utils.sequence_analysis import process_fasta, perform_alignment
from utils.visualization import create_phylogenetic_tree
from utils.sample_data import get_sample_data
from utils.sequence_stats import calculate_sequence_stats, compare_sequences
import plotly.graph_objects as go # Added import for heatmap

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

    st.markdown("""
    Explore and analyze pathogen sequences through multiple analysis methods:
    - Basic sequence statistics
    - Sequence similarity comparison
    - Phylogenetic relationships
    """)

    # Sidebar controls
    st.sidebar.header("Controls")
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
        # Analysis tabs
        tab1, tab2, tab3 = st.tabs([
            "Sequence Statistics", 
            "Sequence Comparison", 
            "Phylogenetic Analysis"
        ])

        with tab1:
            st.header("Basic Sequence Statistics")
            stats_df = calculate_sequence_stats(sequences)
            st.dataframe(stats_df, use_container_width=True)

            # Display summary statistics
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Number of Sequences", len(sequences))
                st.metric("Average Length", int(stats_df['Length'].mean()))
            with col2:
                st.metric("Average GC Content", f"{stats_df['GC Content (%)'].mean():.2f}%")

        with tab2:
            st.header("Sequence Similarity Matrix")
            similarity_df = compare_sequences(sequences)
            st.dataframe(similarity_df, use_container_width=True)

            # Add heatmap visualization
            st.plotly_chart(
                create_similarity_heatmap(similarity_df),
                use_container_width=True
            )

        with tab3:
            st.header("Phylogenetic Analysis")
            try:
                alignment = perform_alignment(sequences) #removed method=alignment_method
                fig = create_phylogenetic_tree(alignment) #removed tree_type=visualization_type
                st.plotly_chart(fig, use_container_width=True)
            except Exception as e:
                st.warning("Phylogenetic tree generation is currently under maintenance. Please check back later.")
                st.info("You can still explore sequence statistics and comparisons in the other tabs.")

def create_similarity_heatmap(similarity_df):
    """Create a heatmap visualization of sequence similarities"""
    import plotly.graph_objects as go

    fig = go.Figure(data=go.Heatmap(
        z=similarity_df.values,
        x=similarity_df.columns,
        y=similarity_df.index,
        colorscale='Viridis',
        zmin=0,
        zmax=100
    ))

    fig.update_layout(
        title="Sequence Similarity Heatmap",
        xaxis_title="Sequence ID",
        yaxis_title="Sequence ID",
        height=500
    )

    return fig

if __name__ == "__main__":
    main()