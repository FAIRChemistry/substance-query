from typing import Dict, Optional
from uuid import uuid4

import sdRDM
from lxml.etree import _Element
from pydantic import PrivateAttr, model_validator
from pydantic_xml import attr, element
from sdRDM.base.datatypes import Identifier
from sdRDM.base.listplus import ListPlus
from sdRDM.base.utils import forge_signature
from sdRDM.tools.utils import elem2dict


@forge_signature
class PreparationStep(
    sdRDM.DataModel,
    search_mode="unordered",
):
    """Specification of a preparation step."""

    id: Optional[str] = attr(
        name="id",
        alias="@id",
        description="Unique identifier of the given object.",
        default_factory=lambda: str(uuid4()),
    )

    label: Optional[str] = element(
        description="Label of the step.",
        default=None,
        tag="label",
        json_schema_extra=dict(),
    )

    preparation_id: Optional[Identifier] = element(
        description="Identifier of the preparation step.",
        default=None,
        tag="preparation_id",
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
