# This demonstrates one way to generate an initial alignment between two
# PDB sequences. It can later be edited by hand.

# Set Modeller environment (including search patch for model.read())
from modeller import *
env = environ()
env.io.atom_files_directory = [".", "../atom_files/"]

# Create a new empty alignment and model:
aln = alignment(env)
mdl = model(env)

# Read the whole 5z82 atom file
code='target_native_structure_5z82'
mdl.read(file=code, model_segment=('FIRST:@', 'END:'))

# Add the model sequence to the alignment
aln.append_model(mdl, align_codes=code, atom_files=code)

# Read 4i1a atom file, and add to alignment
code='template_4i1a'
mdl.read(file=code, model_segment=('FIRST:@', 'END:'))
aln.append_model(mdl, align_codes=code, atom_files=code)

# Align them by sequence
aln.malign(gap_penalties_1d=(-500, -300))
aln.write(file='fer1-seq.ali')

# Align them by structure
aln.malign3d(gap_penalties_3d=(0.0, 2.0))

# check the alignment for its suitability for modeling
aln.check()

aln.write(file='T0951.ali')
