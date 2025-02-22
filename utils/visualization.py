import plotly.graph_objects as go
from Bio import Phylo
from io import StringIO
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def create_phylogenetic_tree(alignment, tree_type="Circular"):
    """
    Create an interactive phylogenetic tree visualization

    Args:
        alignment: MultipleSeqAlignment object
        tree_type: String indicating tree layout type

    Returns:
        plotly.graph_objects.Figure: Interactive tree visualization
    """
    try:
        # Initialize distance calculator explicitly
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)

        # Construct tree
        constructor = DistanceTreeConstructor()
        tree = constructor.build_tree(dm)

        # Convert tree to newick format
        tree_string = StringIO()
        Phylo.write(tree, tree_string, 'newick')

        # Create plotly figure
        fig = go.Figure()

        # Extract coordinates for visualization
        coords = []
        for i, clade in enumerate(tree.get_terminals()):
            x = float(clade.branch_length or 0)
            y = i
            coords.append((x, y, clade.name))

        xs, ys, labels = zip(*coords)

        if tree_type == "Circular":
            # Convert to polar coordinates
            theta = np.linspace(0, 2*np.pi, len(xs))
            r = np.array(xs)

            fig.add_trace(go.Scatterpolar(
                r=r,
                theta=theta,
                text=labels,
                mode='markers+text',
                name='Species'
            ))

            fig.update_layout(
                polar=dict(
                    radialaxis=dict(visible=True, title="Genetic Distance"),
                    angularaxis=dict(visible=True)
                )
            )
        else:
            # Rectangular layout
            fig.add_trace(go.Scatter(
                x=xs,
                y=ys,
                text=labels,
                mode='markers+text',
                name='Species'
            ))

            fig.update_layout(
                xaxis_title="Genetic Distance",
                yaxis_title="Species",
                showlegend=False
            )

        fig.update_layout(
            title="Phylogenetic Tree",
            height=600,
            template="plotly_white"
        )

        return fig
    except Exception as e:
        raise Exception(f"Error creating phylogenetic tree: {str(e)}")