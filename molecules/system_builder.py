"""Methods used to create systems from job statepoint."""
import mbuild as mb
from mbuild.lib.molecules.water import WaterSPC

from reproducibility_project.src.molecules.benzene_ua import BenzeneUA
from reproducibility_project.src.molecules.ethanol_aa import EthanolAA
from reproducibility_project.src.molecules.methane_ua import MethaneUA
from reproducibility_project.src.molecules.pentane_ua import PentaneUA

from reproducibility_project.src.molecules.tip3p import tip3p
from reproducibility_project.src.molecules.tip3p_ew_b import tip3p_ew_b
from reproducibility_project.src.molecules.tip3p_ew_f import tip3p_ew_f
from reproducibility_project.src.molecules.tips3p import tips3p
from reproducibility_project.src.molecules.opc3 import opc3
from reproducibility_project.src.molecules.a99SB_disp import a99SB_disp
from reproducibility_project.src.molecules.opc import opc
from reproducibility_project.src.molecules.tip4p_2005 import tip4p_2005
from reproducibility_project.src.molecules.tip4p_d import tip4p_d
from reproducibility_project.src.molecules.tip4p_ew import tip4p_ew


def construct_system(
    sp,
    scale_liq_box=1.0,
    scale_vap_box=1.0,
    constrain=False,
    fix_orientation=False,
):
    """Construct systems according to job statepoint.

    Parameters
    ----------
    sp: dict (from job.sp)
        Dictionary contains information necessary to construct a system.
        Stored as state of job. The dictionary should resemble:
        {"molecule": str,
         "engine": str,
         "replica": int,
         "temperature": float (in K),
         "pressure": float (in kPa),
         "ensemble": str,
         "N_liquid": int,
         "N_vap": int,
         "box_L_liq": int (nm),
         "box_L_vap", int (nm),
         "init_liq_den": float (g/cm3),
         "init_vap_den": float (g/cm3),
         "mass": float (g/mol),
         "forcefield_name": str,
         "cutoff_style": str,
         "r_cut": float (in nm)}
    scale_liq_box: float, optional, default=1.0
        Option to scale sizes of the liquid box.
    scale_vap_box: float, optional, default=1.0
        Option to scale sizes of the vapor box.
    constrain: boolean, optional, default=False
        Option to use constrainmol on bond lengths to make them exactly equal
        to the bond lengths in the FF file.
    fix_orientation : boolean, optional, default=False
        Option to specify that compounds should not be rotated when filling the
        box. Useful for rigid bodies in HOOMD.

    Returns
    -------
    [filled_liq_box, filled_vap_box]
        Return list of system as specified.
    """
    molecule = get_molecule(sp)
    liq_box = mb.Box([sp["box_L_liq_x"] * scale_liq_box, sp["box_L_liq_y"] * scale_liq_box, sp["box_L_liq_z"] * scale_liq_box])
    if(sp.pdbid):
        filled_liq_box = mb.fill_box(
            compound=[molecule],
            density=[sp["init_liq_den"]],
            box=liq_box,
            fix_orientation=fix_orientation,
        )
    else:
        filled_liq_box = mb.fill_box(
            compound=[molecule],
            n_compounds=[sp["N_liquid"]],
            box=liq_box,
            fix_orientation=fix_orientation,
        )

    if sp["box_L_vap"] and sp["N_vap"]:
        vap_box = mb.Box([sp["box_L_vap"] * scale_vap_box] * 3)
        filled_vap_box = mb.fill_box(
            compound=[molecule], n_compounds=[sp["N_vap"]], box=vap_box
        )
        boxes = [filled_liq_box, filled_vap_box]
    else:
        boxes = [filled_liq_box, None]

    if not constrain or sp.molecule == "methaneUA":
        return boxes

    # If we reached this far, we need constrainmol, foyer, and load_ff
    import foyer
    from constrainmol import ConstrainedMolecule

    from reproducibility_project.src.utils.forcefields import load_ff

    ff = load_ff(sp.forcefield_name)
    parmed_molecule = molecule.to_parmed()
    typed_molecule = ff.apply(parmed_molecule)
    constrain_mol = ConstrainedMolecule(typed_molecule)

    for box in boxes:
        if box is None:
            continue
        else:
            for mol in box.children:
                constrain_mol.update_xyz(mol.xyz * 10)  # nm to angstrom
                constrain_mol.solve()
                mol.xyz = constrain_mol.xyz / 10.0  # angstrom to nm

    return boxes

def construct_ion_system(
    sp,
    scale_liq_box=1.0,
    scale_vap_box=1.0,
    constrain=False,
    fix_orientation=False,
):
    """Construct systems according to job statepoint.

    Parameters
    ----------
    sp: dict (from job.sp)
        Dictionary contains information necessary to construct a system.
        Stored as state of job. The dictionary should resemble:
        {"molecule": str,
         "engine": str,
         "replica": int,
         "temperature": float (in K),
         "pressure": float (in kPa),
         "ensemble": str,
         "N_liquid": int,
         "N_vap": int,
         "box_L_liq": int (nm),
         "box_L_vap", int (nm),
         "init_liq_den": float (g/cm3),
         "init_vap_den": float (g/cm3),
         "mass": float (g/mol),
         "forcefield_name": str,
         "cutoff_style": str,
         "r_cut": float (in nm)}
    scale_liq_box: float, optional, default=1.0
        Option to scale sizes of the liquid box.
    scale_vap_box: float, optional, default=1.0
        Option to scale sizes of the vapor box.
    constrain: boolean, optional, default=False
        Option to use constrainmol on bond lengths to make them exactly equal
        to the bond lengths in the FF file.
    fix_orientation : boolean, optional, default=False
        Option to specify that compounds should not be rotated when filling the
        box. Useful for rigid bodies in HOOMD.

    Returns
    -------
    [filled_liq_box, filled_vap_box]
        Return list of system as specified.
    """

    FF_file_cation = 'oplsaa'
    FF_file_anion = 'oplsaa'

    cation = mb.load('CCO', smiles=True)
    ethanol.name = 'ETO'
    ethanol.energy_minimize(forcefield=FF_file_ethanol, steps=10**5)

    FF_dict = {water.name: FF_file_water, ethanol.name: FF_file_ethanol}

    residues_list = [ethanol.name, water.name]

    fix_bonds_angles_residues = [water.name]

    molecule = get_molecule(sp)
    liq_box = mb.Box([sp["box_L_liq_x"] * scale_liq_box, sp["box_L_liq_y"] * scale_liq_box, sp["box_L_liq_z"] * scale_liq_box])
    if(sp.pdbid):
        filled_liq_box = mb.fill_box(
            compound=[molecule],
            density=[sp["init_liq_den"]],
            box=liq_box,
            fix_orientation=fix_orientation,
        )
    else:
        filled_liq_box = mb.fill_box(
            compound=[molecule],
            n_compounds=[sp["N_liquid"]],
            box=liq_box,
            fix_orientation=fix_orientation,
        )

    if sp["box_L_vap"] and sp["N_vap"]:
        vap_box = mb.Box([sp["box_L_vap"] * scale_vap_box] * 3)
        filled_vap_box = mb.fill_box(
            compound=[molecule], n_compounds=[sp["N_vap"]], box=vap_box
        )
        boxes = [filled_liq_box, filled_vap_box]
    else:
        boxes = [filled_liq_box, None]

    if not constrain or sp.molecule == "methaneUA":
        return boxes

    # If we reached this far, we need constrainmol, foyer, and load_ff
    import foyer
    from constrainmol import ConstrainedMolecule

    from reproducibility_project.src.utils.forcefields import load_ff

    ff = load_ff(sp.forcefield_name)
    parmed_molecule = molecule.to_parmed()
    typed_molecule = ff.apply(parmed_molecule)
    constrain_mol = ConstrainedMolecule(typed_molecule)

    for box in boxes:
        if box is None:
            continue
        else:
            for mol in box.children:
                constrain_mol.update_xyz(mol.xyz * 10)  # nm to angstrom
                constrain_mol.solve()
                mol.xyz = constrain_mol.xyz / 10.0  # angstrom to nm

    return boxes


def get_molecule(sp):
    """Construct the mbuild molecule for the job statepoint.

    Parameters
    ----------
    sp: dict (from job.sp)
        Dictionary contains information necessary to construct a system.
        Stored as state of job. The dictionary should resemble:
        {"molecule": str,
         "engine": str,
         "replica": int,
         "temperature": float (in K),
         "pressure": float (in kPa),
         "ensemble": str,
         "N_liquid": int,
         "N_vap": int,
         "box_L_liq": int (nm),
         "box_L_vap", int (nm),
         "init_liq_den": float (g/cm3),
         "init_vap_den": float (g/cm3),
         "mass": float (g/mol),
         "forcefield_name": str,
         "cutoff_style": str,
         "r_cut": float (in nm)}

    Returns
    -------
    molecule
        Return mBuild molecule for the statepoint.
    """
    molecule_dict = {
        "methaneUA": MethaneUA(),
        "pentaneUA-flexible_bonds": PentaneUA(),
        "pentaneUA-constrain_bonds": PentaneUA(),
        "pentaneUA": PentaneUA(),
        "benzeneUA": BenzeneUA(),
        "waterSPCE": WaterSPC(),
        "ethanolAA": EthanolAA(),
        "spce": WaterSPC(),
        "tip3p": tip3p(),
        "tip3p_ew_b": tip3p_ew_f(),
        "tip3p_ew_f": tip3p_ew_f(),
        "tips3p": tips3p(),
        "opc3": opc3(),
        "tip4p_ew": tip4p_ew(),
        "tip4p_2005": tip4p_2005(),
        "tip4p_d": tip4p_d(),
        "a99SB_disp": a99SB_disp(),
        "opc": opc(),
    }
    molecule = molecule_dict[sp["molecule"]]
    # For now, hardcode all molecule names to WAT
    #molecule.name = sp["molecule"]
    molecule.name = "WAT"
    return molecule
