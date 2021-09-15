"""
Use this script to update the thermo and kinetic parameters of the kinetic mechanism with the most recent version
of RMG-database. The script creates a new folder to store the new model with updated parameters.
"""
import os

import numpy as np

from rmgpy.chemkin import save_chemkin
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.kinetics.depository import DepositoryReaction
from rmgpy.kinetics import KineticsData, Arrhenius
from rmgpy.tools.loader import load_rmg_job
from rmgpy.thermo.thermoengine import submit


# path to the folder storing input.py, chem_annotated.inp, and species_dictionary.txt
modelpath = 'original'
# load in the model
rmg = load_rmg_job(input_file=os.path.join(modelpath, 'input.py'),
                   chemkin_file=os.path.join(modelpath, 'chem_annotated.inp'),
                   species_dict=os.path.join(modelpath, 'species_dictionary.txt'),
                   generate_images=False)
# load the database
rmg.load_database()

# set each species thermo to None and then obtain its updated thermo from the database
for spc in rmg.reaction_model.core.species:
    spc.thermo = None
    submit(spc, solvent_name=rmg.reaction_model.solvent_name)

# create new directory to store the output
os.mkdir(os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo'))
chemkin_path = os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo', 'chem_updated_thermo.inp')
chemkin_verbose_path = os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo', 'chem_annotated_updated_thermo.inp')
dictionary_path = os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo', 'species_dictionary.txt')
save_chemkin(rmg.reaction_model,
             chemkin_path,
             chemkin_verbose_path,
             dictionary_path,
             save_edge_species=False)

for i, reaction in enumerate(rmg.reaction_model.core.reactions):
    # get the new reaction template
    reaction.template = rmg.database.kinetics.generate_reactions_from_families(reactants=reaction.reactants,
                                                                               products=reaction.products,
                                                                               only_families=[reaction.family])[0].template
    reaction.kinetics = None
    rmg.reaction_model.apply_kinetics_to_reaction(reaction)
    if isinstance(reaction.kinetics, KineticsData):
        reaction.kinetics = reaction.kinetics.to_arrhenius()

    #  correct barrier heights of estimated kinetics
    if isinstance(reaction, TemplateReaction) or isinstance(reaction, DepositoryReaction):
        reaction.fix_barrier_height()

    if rmg.reaction_model.pressure_dependence and reaction.is_unimolecular():
        # if this is going to be run through pressure dependence code, make sure the barrier is positive
        reaction.fix_barrier_height(force_positive=True)

os.mkdir(os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo_kinetics') )
chemkin_path = os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo_kinetics', 'chem_updated_thermo_kinetics.inp')
chemkin_verbose_path = os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo_kinetics', 'chem_annotated_updated_thermo_kinetics.inp')
dictionary_path = os.path.join(os.path.dirname(rmg.output_directory), 'updated_thermo_kinetics', 'species_dictionary.txt')
save_chemkin(rmg.reaction_model,
             chemkin_path,
             chemkin_verbose_path,
             dictionary_path,
             save_edge_species=False)
