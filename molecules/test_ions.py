"""Methods used to create systems from job statepoint."""
import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
from reproducibility_project.src.utils.forcefields import get_ff_path_ion

FF_file_cation = get_ff_path_ion("custom", "sod")
FF_file_anion = get_ff_path_ion("custom", "cla")

cation = mb.load('[Na]', smiles=True)
cation.name = 'SOD'

anion = mb.load('[Cl-]', smiles=True)
anion.name = 'CLA'

num_cations = 2
num_anions = 2


FF_dict = {cation.name: FF_file_cation, anion.name: FF_file_anion}

residues_list = [cation.name, anion.name]

#liq_box = mb.Box([sp["box_L_liq_x"] * scale_liq_box, sp["box_L_liq_y"] * scale_liq_box, sp["box_L_liq_z"] * scale_liq_box])
liq_box = mb.Box([1.0, 1.0, 1.0])


filled_liq_box = mb.fill_box(
    compound=[cation, anion],
    n_compounds=[num_cations,num_anions],
    box=liq_box
)

boxes = [filled_liq_box, None]

charmm = mf_charmm.Charmm(filled_liq_box,
                          'GEMC_NVT_water_ethanol_liq',
                          structure_box_1=None,
                          filename_box_1=None,
                          ff_filename="ions",
                          forcefield_selection=FF_dict,
                          residues=residues_list,
                          bead_to_atom_name_dict=None,
                          fix_residue=None,
                          gomc_fix_bonds_angles=None,
                          reorder_res_in_pdb_psf=True
                          )

charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()






