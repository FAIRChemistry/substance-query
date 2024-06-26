"""This module contains an ipywidget for querying PubChem for a 
chemical substance and returning a research data model object containing
the substance's information.
"""

from os import wait

from IPython.display import display
from click import style
import ipywidgets as widgets
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp

from . import Substance


class SubstanceWidget:
    def __init__(self, substance_data_model: Substance = None) -> Substance:
        self.data_model = substance_data_model
        self.query_compound = ""
        self.compounds = []
        self.selected_compound = None

        # Create the label widget for the title
        title_label = widgets.Label(value="Search for a compound:")

        # Create the text input widget for the search query
        text_input = widgets.Text(placeholder="Enter compound name")

        # Create the search button
        search_button = widgets.Button(description="Search")

        # Create the compound dropdown widget for the search results
        compound_dropdown = widgets.Dropdown(
            options=[],
            description="Select a compound:",
            layout={"width": "max-content"},
            style={"description_width": "initial"},
            disabled=True,
        )

        # Create the ouput widget for compound name and structure
        output = widgets.Output()

        # Define a function to handle the search button click event
        def on_search_button_click(b):
            self.query_compound = text_input.value
            self.compounds = self._query_pubchem(self.query_compound)
            compound_dropdown.options = [
                compound.iupac_name.capitalize() for compound in self.compounds
            ]
            compound_dropdown.value = compound_dropdown.options[0]
            compound_dropdown.disabled = False

        # Attach the function to the search button click event
        search_button.on_click(on_search_button_click)

        # Define a function to handle the dropdown change event
        def on_compound_dropdown_change(change):
            with output:
                output.clear_output()
                compound_name = change.new
                self.selected_compound = self.compounds[
                    compound_dropdown.options.index(compound_name)
                ]
                self._display_compound(self.selected_compound)

        # Attach the function to the dropdown change event
        compound_dropdown.observe(on_compound_dropdown_change, names="value")

        # Create a container for the widgets
        container = widgets.VBox(
            [title_label, text_input, search_button, compound_dropdown, output]
        )

        # Display the container
        display(container)

    # Create the function to query PubChem
    @staticmethod
    def _query_pubchem(compound_name):
        """
        Get a compound by its name.

        Args:
            compound_name: The name of the compound to search for.

        Returns:
            A list of compounds that match the name.
        """
        return pcp.get_compounds(compound_name, "name")

    @staticmethod
    def _display_compound(compound):
        iupac_name = compound.iupac_name
        inchi = compound.inchi
        synonyms = [name.capitalize() for name in compound.synonyms]
        molecule = Chem.MolFromInchi(inchi)

        image = Draw.MolToImage(molecule)

        display(widgets.HTML(value=f"<b>IUPAC name:</b>"))
        display(widgets.HTML(value=f"{iupac_name.capitalize()}"))
        display(widgets.HTML(value=f"<b>Structure:</b>"))
        display(image)
        display(widgets.HTML(value=f"<b>Other names:</b>"))
        for name in synonyms:
            display(widgets.HTML(value=f"{name}"))
