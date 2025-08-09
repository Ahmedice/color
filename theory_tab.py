import streamlit as st
import plotly.graph_objects as go
import numpy as np
from scipy.interpolate import interp1d

class TheoryTab:
    def __init__(self):
        self.render()

    def render(self):
        st.header("Theory Tab Content")

        # Sample data based on the image
        x = [230, 260, 280]
        y = [18, 36, 18]

        # Add data points to connect to zero and extend to 340
        x = [220] + x + [340]
        y = [0] + y + [0]

        # Create a spline interpolation
        f = interp1d(x, y, kind='quadratic')

        # Generate x values for the curve
        xnew = np.linspace(min(x), max(x), num=100, endpoint=True)

        # Generate y values for the curve
        ynew = f(xnew)

        # Create the plot
        fig = go.Figure()

        # Add the scatter trace
        fig.add_trace(go.Scatter(x=xnew, y=ynew, mode='lines', line=dict(color='red')))

        # Add labels for the wavelengths
        fig.add_trace(go.Scatter(x=[230, 260, 280], y=[18, 36, 18], mode='markers+text', text=["Salt", "DNA/RNA", "Protein"], textposition="top center",
                                 textfont=dict(size=16)))


        # Set the axis labels
        fig.update_layout(
            xaxis_title="Wavelength (nm)",
            yaxis_title="10mm Absorbance",
            xaxis_range=[220, 340],
            yaxis_range=[0, 40]
        )

        # Display the plot in Streamlit
        st.plotly_chart(fig, key="initial_chart")

        # Add sliders for 230, 260, and 280 wavelengths
        salt = st.slider("Salt", 0, 40, 18)
        dna_rna = st.slider("DNA/RNA", 0, 40, 36)
        protein = st.slider("Protein", 0, 40, 18)

        # Update the plot based on slider values
        y_updated = [salt, dna_rna, protein]

        # Add data points to connect to zero and extend to 340
        y_updated = [0] + y_updated + [0]

        # Create a new spline interpolation with updated values
        f_updated = interp1d(x, y_updated, kind='quadratic')
        ynew_updated = f_updated(xnew)

        fig.data[0].y = ynew_updated
        st.plotly_chart(fig, key="updated_chart")