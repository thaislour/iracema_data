from rdkit import Chem
from rdkit.Chem import inchi, Crippen, rdMolDescriptors, Descriptors, AllChem, rdFingerprintGenerator
from rdkit.Chem import SDWriter
from io import StringIO
from rdkit.Chem.rdmolfiles import MolToPDBBlock

def validate_smiles(smiles):
    if not smiles or not isinstance(smiles, str):
        raise ValueError("Invalid SMILES string")

def smiles_to_mol(smiles):
    validate_smiles(smiles)
    rdkit_mol = Chem.MolFromSmiles(smiles)
    if rdkit_mol is None:
        raise ValueError("Invalid SMILES string")
    return rdkit_mol

def smiles_to_canon_smiles(smiles):
    try:
        canon = Chem.CanonSmiles(smiles)
    except Exception as ex:
        raise ValueError(str(ex))
    return canon


def mol_to_inchi(rdkit_mol):
    try:
        inchi_str = inchi.MolToInchi(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return inchi_str


def inchi_to_inchikey(inchi_str):
    try:
        inchikey = inchi.InchiToInchiKey(inchi_str)
        if inchikey is None:
            raise ValueError("Empty inchikey return, check inchi string validity")
    except Exception as ex:
        raise ValueError(str(ex))
    return inchikey


def mol_to_inchikey(rdkit_mol):
    try:
        inchikey = inchi.MolToInchiKey(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return inchikey


def calculate_logp(rdkit_mol):
    try:
        logp = Crippen.MolLogP(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return logp

def calculate_tpsa(rdkit_mol):
    try:
        tpsa = rdMolDescriptors.CalcTPSA(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return tpsa


def calculate_mol_formula(rdkit_mol):
    try:
        mol_formula = rdMolDescriptors.CalcMolFormula(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return mol_formula


def calculate_exact_mol_wt(rdkit_mol):
    try:
        exact_mol_wt = rdMolDescriptors.CalcExactMolWt(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return exact_mol_wt


def calculate_mol_wt(rdkit_mol):
    try:
        mol_wt = Descriptors.MolWt(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return mol_wt

def calculate_coordinates_2d(rdkit_mol):
    try:
        AllChem.Compute2DCoords(rdkit_mol)
        coordinates_2d = Chem.MolToMolBlock(rdkit_mol)
    except Exception as ex:
        raise ValueError(str(ex))
    return coordinates_2d

def calculate_coordinates_3d(rdkit_mol):
    try:
        rdkit_mol_3d = Chem.AddHs(rdkit_mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xf00d
        AllChem.EmbedMolecule(rdkit_mol_3d, params)
        AllChem.UFFOptimizeMolecule(rdkit_mol_3d)
        coordinates_3d = Chem.MolToMolBlock(rdkit_mol_3d)
    except Exception as ex:
        raise ValueError(str(ex))
    return coordinates_3d

def calculate_coordinates_3d_mol(rdkit_mol):
    try:
        rdkit_mol_3d = Chem.AddHs(rdkit_mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xf00d
        AllChem.EmbedMolecule(rdkit_mol_3d, params)
        AllChem.UFFOptimizeMolecule(rdkit_mol_3d)
    except Exception as ex:
        raise ValueError(str(ex))
    return rdkit_mol_3d

def get_morgan_fingerprint(identifier, radius=2, fpSize=1024):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize)
    try:
        if isinstance(identifier, str):
            rdkit_mol = smiles_to_mol(identifier)
        elif isinstance(identifier, Chem.Mol):
            rdkit_mol = identifier
        else:
            raise ValueError("Input must be a SMILES string or an RDKit Mol object")
        fingerprint = mfpgen.GetFingerprint(rdkit_mol)
    except Exception as ex:
        raise ValueError(f"Error generating fingerprint: {str(ex)}")
    return fingerprint


def get_pdb_from_mol(rdkit_mol):
    return MolToPDBBlock(rdkit_mol)


def generate_sdf(molecules, use_3d=False):
    sdf_strings = []
    for mol_data in molecules:
        mol_id = mol_data['id']
        smiles = mol_data['smiles']
        name = mol_data['name']
        try:
            rdkit_mol = smiles_to_mol(smiles)
            
            # If use_3d is True, calculate 3D coordinates, otherwise calculate 2D coordinates
            if use_3d:
                try:
                    rdkit_mol = calculate_coordinates_3d_mol(rdkit_mol)
                except Exception as ex:
                    raise ValueError(str(ex))
            else:
                calculate_coordinates_2d(rdkit_mol)

            rdkit_mol.SetProp("_Name", name)
            rdkit_mol.SetProp("ID", mol_id)

            # Write the molecule to a string in SDF format
            sio = StringIO()
            SDWriter(sio).write(rdkit_mol)
            sdf_strings.append(sio.getvalue())
            sio.close()

        except Exception as ex:
            raise ValueError(str(ex))

    return "\n".join(sdf_strings)

def generateData(smiles):

        # Rdkit generators
        rdkit_mol = smiles_to_mol(smiles)
        canon_smiles = smiles_to_canon_smiles(smiles)
        inchi_str = mol_to_inchi(rdkit_mol)
        inchikey = inchi_to_inchikey(inchi_str)
        logp = calculate_logp(rdkit_mol)
        tpsa = calculate_tpsa(rdkit_mol)
        mol_formula = calculate_mol_formula(rdkit_mol)
        exact_mol_wt = calculate_exact_mol_wt(rdkit_mol)
        mol_wt = calculate_mol_wt(rdkit_mol)
        coordinates_2d = calculate_coordinates_2d(rdkit_mol)
        coordinates_3d = calculate_coordinates_3d(rdkit_mol)
        pdb_file = get_pdb_from_mol(rdkit_mol)
        
        data = {"mol" : str(rdkit_mol), 
                "canon_smiles" : canon_smiles, 
                "inchi_str" : inchi_str, 
                "inchikey" : inchikey, 
                "logp" : logp, 
                "tpsa" : tpsa, 
                "mol_formula" : mol_formula, 
                "exact_mol_wt" : exact_mol_wt, 
                "mol_wt" : mol_wt, 
                "coordinates_2d" : coordinates_2d, 
                "coordinates_3d" : coordinates_3d,
                "pdb" : pdb_file}
        return data