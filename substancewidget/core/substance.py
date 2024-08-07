from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
class Substance(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Specification of a chemical substance."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    label: Optional[str] = element(
        description="Any name or identifier of the substance.",
        default=None,
        tag="label",
        json_schema_extra=dict(),
    )

    iupac_name: Optional[str] = element(
        description="IUPAC name of the substance.",
        default=None,
        tag="iupac_name",
        json_schema_extra=dict(),
    )

    canonical_smiles: Optional[str] = element(
        description=(
            "Simplified molecular-input line-entry system (SMILES) notation of the"
            " substance."
        ),
        default=None,
        tag="canonical_smiles",
        json_schema_extra=dict(
            regex="/^([^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,})$/ig",
        ),
    )

    inchi: Optional[str] = element(
        description="IUPAC International Chemical Identifier (InChI) of the substance.",
        default=None,
        tag="inchi",
        json_schema_extra=dict(
            regex="/^((InChI=)?[^J][0-9BCOHNSOPrIFla+\-\(\)\\\/,pqbtmsih]{6,})$/ig",
        ),
    )

    inchi_key: Optional[str] = element(
        description=(
            "IUPAC International Chemical Identifier (InChI) key of the substance."
        ),
        default=None,
        tag="inchi_key",
        json_schema_extra=dict(
            regex="/^([0-9A-Z\-]+)$/",
        ),
    )

    molecular_weight: Optional[float] = element(
        description="Molecular weight of the substance in g/mol.",
        default=None,
        tag="molecular_weight",
        json_schema_extra=dict(),
    )

    lot_number: Optional[str] = element(
        description="Identifier of the lot of the substance.",
        default=None,
        tag="lot_number",
        json_schema_extra=dict(),
    )

    manufacturer: Optional[str] = element(
        description="Manufacturer of the substance.",
        default=None,
        tag="manufacturer",
        json_schema_extra=dict(),
    )

    _repo: Optional[str] = PrivateAttr(
        default="https://github.com/FAIRChemistry/substance-widget"
    )
    _commit: Optional[str] = PrivateAttr(
        default="f33d876ffc611df676fd4b9759543aa38f666dcd"
    )

    _raw_xml_data: Dict = PrivateAttr(default_factory=dict)

    @model_validator(mode="after")
    def _parse_raw_xml_data(self):
        for attr, value in self:
            if isinstance(value, (ListPlus, list)) and all(
                isinstance(i, _Element) for i in value
            ):
                self._raw_xml_data[attr] = [elem2dict(i) for i in value]
            elif isinstance(value, _Element):
                self._raw_xml_data[attr] = elem2dict(value)

        return self
