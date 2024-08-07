"""This module contains an ipywidget for querying PubChem for a 
chemical substance and returning a research data model object containing
the substance's information.
"""

from os import wait
import re

from IPython.display import display
from click import style
import ipywidgets as widgets
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pcp

from . import Substance


RE_SMILES = re.compile(r"/^([^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,})$/ig")
RE_INCHI = re.compile(
    r"/^((InChI=)?[^J][0-9BCOHNSOPrIFla+\-\(\)\\\/,pqbtmsih]{6,})$/ig"
)
RE_INCHIKEY = re.compile(r"/^([0-9A-Z\-]+)$/")


class SubstanceWidget:
    def __init__(self, substance_data_model: Substance = None) -> Substance:
        if not substance_data_model:
            substance_data_model = Substance()
        self.substance = substance_data_model
        self.query = ""
        self.compounds = []
        self.selected_compound = None

        # Create the label widget for the title
        title_label = widgets.Label(value="Search for a compound:")

        # Create the text input widget for the search query
        text_input = widgets.Text(placeholder="Enter name, CID, SMILES, or InChI(Key)")

        # Create the search button
        search_button = widgets.Button(description="Search PubChem")

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
            del self.compounds[:]
            self.query = text_input.value
            match self.query:
                case self.query if self.query.isdigit():
                    self.compounds.append(pcp.Compound.from_cid(self.query))
                case self.query if RE_SMILES.match(self.query):
                    self.compounds = pcp.get_compounds(self.query, "smiles")
                case self.query if RE_INCHI.match(self.query):
                    self.compounds = pcp.get_compounds(self.query, "inchi")
                case self.query if RE_INCHIKEY.match(self.query):
                    self.compounds = pcp.get_compounds(self.query, "inchikey")
                case _:
                    self.compounds = pcp.get_compounds(self.query, "name")
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
                self.substance.id = "https://pubchem.ncbi.nlm.nih.gov/compound/" + str(
                    self.selected_compound.cid
                )
                self.substance.iupac_name = (
                    self.selected_compound.iupac_name.capitalize()
                )
                self.substance.canonical_smiles = (
                    self.selected_compound.canonical_smiles
                )
                self.substance.inchi_key = self.selected_compound.inchikey
                self.substance.molecular_weight = (
                    self.selected_compound.molecular_weight
                )

        # Attach the function to the dropdown change event
        compound_dropdown.observe(on_compound_dropdown_change, names="value")

        # Create a container for the widgets
        container = widgets.VBox(
            [title_label, text_input, search_button, compound_dropdown, output]
        )

        # Display the container
        display(container)

    def get_substance_data_model(self) -> Substance:
        return self.substance

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
        for i, name in enumerate(synonyms):
            if i >= 9:
                display(widgets.HTML(value="..."))
                break
            display(widgets.HTML(value=f"{name}"))
            i += 1
