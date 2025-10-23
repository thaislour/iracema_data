# Data and Software Repository for "IRACEMA: A Database Management System for Bioactive Compounds Obtained and Characterized by Brazilian Researchers"

This repository contains the dataset and software scripts used in the manuscript submitted to the Journal of Chemical Information and Modeling (JCIM).

## 1. Reproducibility Requirements

The goal of this repository is to enable full reproducibility of the results presented in the iracema platform.

* Software / Environment: The software scripts were developed in Python. To reproduce the computational environment, use the `requirements.txt` file located at the repository root.
* Recommended minimum Python version: 3.10 or newer.
* Tools used: Molecular properties were calculated using the following free, publicly available libraries:
	* RDKit (version listed in `requirements.txt`: 2024.3.2) - https://www.rdkit.org/
	* NumPy (RDKit dependency and numerical processing) - https://numpy.org/

## Data and Software Availability

All data and scripts necessary to reproduce the results reported in the manuscript are provided in this repository as machine-readable files. Below we group the resources and provide provenance and reproducibility notes.

1) Data files (repository root)
The primary data files are located at the repository root:

| File | Description |
| :--- | :--- |
| `iracema_predicted_data.json` | Computed molecular properties to each compound shown in the platform such as InChI, InChIKey, molecular formula, TPSA, exact molecular weight, MolWt, logP and additional ADME properties|
| `iracema_database_activities.csv` | Activity-level export of the IRACEMA database (one row per activity; includes SMILES, activity value/type, DOI and assay metadata). Contains curated molecule information extracted from the literature and its reported biological activities. **Each row represents an activity.** |
| `iracema_database_molecules.csv` | Molecule-level export of the IRACEMA database (one row per unique molecule; includes IRACEMA_ID, SMILES, InChI/InChIKey, formula, RDKit-derived columns and publication metadata). **Each row represents a molecule and its computed properties.**|

### `iracema_predicted_data.json` — format and schema

`iracema_predicted_data.json` is a machine-readable JSON mapping where each top-level key is the input SMILES string (or canonical SMILES) for a molecule and the corresponding value is an object containing grouped properties. The main groups present in the file are:

- `pharmacokinetics`: pharmacokinetic and ADME-related outputs (these were generated using the SwissADME web tool). Example fields: LogKp (skin permeation), absorption flags, P-gp substrate predictions, CYP inhibition flags, etc.
- `druglikeness`: medicinal-chemistry / drug-likeness summaries (rules-of-5 flags, PAINS, synthetic accessibility score, and similar aggregated flags).
- `rdKitData`: RDKit-derived properties and identifiers (InChI, InChIKey, molecular formula, TPSA, exact molecular weight, MolWt, LogP computed by RDKit, 2D/3D coordinates, and fingerprint objects).

Example — quick Python access pattern

```python
import json

with open('predict_properties_data.json') as fh:
	data = json.load(fh)

# pick a molecule by SMILES key
smiles = next(iter(data))
mol = data[smiles]
logp = mol.get('rdKitData', {}).get('logp')
sa = mol.get('druglikeness', {}).get('SyntheticAccessibility')
swiss_kp = mol.get('pharmacokinetics', {}).get('LogKp_skin_permeation')
print(smiles, logp, sa, swiss_kp)
```

2) Scripts and calculation code (`src`)
- `src/properties_calculation_script.py` — RDKit-based Python functions used to compute many of the molecular properties reported in `predict_properties_data.json` (InChI/InChIKey, logP, TPSA, formulas/masses, 2D/3D coordinates, SDF/PDB export, fingerprints). This file documents the exact functions used and is sufficient to reproduce the RDKit-derived properties when run on the same input SMILES.

3) External tools used
- SwissADME (https://www.swissadme.ch/) — Pharmacokinetic and medicinal-chemistry properties (for example, qualitative absorption, P-gp substrate prediction, CYP inhibition flags, synthetic accessibility, and other cheminformatics flags) were obtained using the free SwissADME web tool. SwissADME is publicly accessible and was used as an external, web-based resource.

4) Environment and reproducibility
- `requirements.txt` (repository root) lists the Python packages and pinned versions used in the analysis. 

5) Limitations and confidential data
- No confidential data required for reproducing the essential results are withheld. Any confidential information (if present in other materials) is not necessary to reproduce the main results reported in the manuscript.