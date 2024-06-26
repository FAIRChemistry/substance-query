# Substance research data model

Python object model specifications based on the [software-driven-rdm](https://github.com/JR-1991/software-driven-rdm) Python library.

## Core objects

### Substance

Specification of a chemical substance.

- label
  - Type: string
  - Description: Any name or identifier of the substance.
- substance_id
  - Type: Identifier
  - Description: Identifier of the substance.  
- iupac_name
  - Type: string
  - Description: IUPAC name of the substance.
- canonical_smiles
  - Type: string
  - Description: Simplified molecular-input line-entry system (SMILES) notation of the substance.
  - Regex: /^([^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,})$/ig
- inchi_key
  - Type: string
  - Description: IUPAC International Chemical Identifier (InChI) key of the substance.
  - Regex: /^([0-9A-Z\-]+)$/
- molecular_weight
  - Type: float
  - Description: Molecular weight of the substance in g/mol.
- lot_number
  - Type: string
  - Description: Identifier of the lot of the substance.
- preparation_procedure
  - Type: [PreparationProcedure](#preparationprocedure)
  - Description: Procedure used to prepare the substance.
- analytical_data
  - Type: [AnalyticalData](#analyticaldata)
  - Description: Analytical data of the substance.
  - Multiple: True

### PreparationProcedure

Specification of a preparation procedure.

- preparation_description
  - Type: string
  - Description: Description of the preparation procedure.
- preparation_steps
  - Type: [PreparationStep](#preparationstep)
  - Description: Steps of the preparation procedure.

### PreparationStep

Specification of a preparation step.

- label
  - Type: string
  - Description: Label of the step.
- preparation_id
  - Type: Identifier
  - Description: Identifier of the preparation step.

### AnalyticalData

Specification of analytical data.

- label
  - Type: string
  - Description: Label of the analytical data.
- analytical_method
  - Type: string
  - Description: Method used to obtain the analytical data.
- analytics_id
  - Type: Identifier
  - Description: Identifier of the analytical data.
