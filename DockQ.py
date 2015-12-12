#!/usr/bin/env python

import Bio.PDB
import sys
import os
import re
from Bio.SVDSuperimposer import SVDSuperimposer

def parse_fnat(fnat_out):
    fnat=0;
    inter=[]
    for line in fnat_out:
#        print line
        line=line.rstrip('\n')
        match=re.search(r'NATIVE: (\d+)(\w) (\d+)(\w)',line)
        if(re.search(r'^Fnat',line)):
            list=line.split(' ')
            fnat=list[3]
        elif(match):
            #print line
            res1=match.group(1)
            chain1=match.group(2)
            res2=match.group(3)
            chain2=match.group(4)
           # print res1 + ' ' + chain1 + ' ' + res2 + ' ' + chain2
            inter.append(res1 + chain1)
            inter.append(res2 + chain2)
    return (fnat,inter)

exec_path=os.path.dirname(os.path.abspath(sys.argv[0]))
model=sys.argv[1]
native=sys.argv[2]

cmd=exec_path + '/fnat ' + model + ' ' + native + ' 5'


fnat_out = os.popen(cmd).readlines()
(fnat,interface)=parse_fnat(fnat_out)
#print inter_pairs
#fnat=grep(re('^Fnat',fnat_out))

          


# Select what residues numbers you wish to align
# and put them in a list
start_id = 100
end_id   = 200
atoms_to_be_aligned = range(start_id, end_id + 1)

# Start the parser
pdb_parser = Bio.PDB.PDBParser(QUIET = True)

# Get the structures
ref_structure = pdb_parser.get_structure("reference", native)
sample_structure = pdb_parser.get_structure("model", model)

# Use the first model in the pdb-files for alignment
# Change the number 0 if you want to align to another structure
ref_model    = ref_structure[0]
sample_model = sample_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
ref_atoms = []
sample_atoms = []

common_interface=[]
chain_res={}

# Do the same for the sample structure
for sample_chain in sample_model:
  chain=sample_chain.id
  if chain not in chain_res.keys():
      chain_res[chain]=[]
  for sample_res in sample_chain:
    resname=sample_res.get_id()[1]
    key=str(resname) + chain
#    print key
    chain_res[chain].append(key)
#    chain_res[chain].append(sample_res)
    if key in interface:
 #     print key  
      sample_atoms.append(sample_res['CA'])
      common_interface.append(key)

#print inter_pairs

chain_ref={}
common_residues=[]
# Iterate of all chains in the model in order to find all residues
for ref_chain in ref_model:
  # Iterate of all residues in each model in order to find proper atoms
#  print dir(ref_chain)
  chain=ref_chain.id
  if chain not in chain_ref.keys():
      chain_ref[chain]=[]
  for ref_res in ref_chain:
      resname=ref_res.get_id()[1]
      key=str(resname) + chain
      #print ref_res
#      print key
     # print chain_res.values()
      if key in chain_res[chain]: # if key is present in sample
          #print key
          chain_ref[chain].append(ref_res['CA'])
          common_residues.append(key)
          #chain_sample.append((ref_res['CA'])
      if key in common_interface:
      # Check if residue number ( .get_id() ) is in the list
      # Append CA atom to list
        #print key  
        ref_atoms.append(ref_res['CA'])

#get the ones that were present in native        
chain_sample={}
for sample_chain in sample_model:
  chain=sample_chain.id
  if chain not in chain_sample.keys():
      chain_sample[chain]=[]
  for sample_res in sample_chain:
    resname=sample_res.get_id()[1]
    key=str(resname) + chain
    if key in common_residues:
        chain_sample[chain].append(sample_res['CA'])
    #if key in common_residues:
 #     print key  
      #sample_atoms.append(sample_res['CA'])
      #common_interface.append(key)

        
#print ref_atoms
#print sample_atoms
#print sample_model.get_atoms()
# Now we initiate the superimposer:
super_imposer = Bio.PDB.Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_model.get_atoms())

# Print RMSD:
irms=super_imposer.rms

(chain1,chain2)=chain_sample.keys()

ligand_chain=chain1
receptor_chain=chain2
if(len(chain_sample[chain1]) > len(chain_sample[chain2])):
    receptor_chain=chain1
    ligand_chain=chain2


#print chain_sample.keys()

#Set to align on receptor
super_imposer.set_atoms(chain_ref[receptor_chain], chain_sample[receptor_chain])
super_imposer.apply(sample_model.get_atoms())
#print dir(super_imposer)
#print chain1_rms
#Set to calculate rms on ligand
super_imposer.set_atoms(chain_ref[ligand_chain], chain_sample[ligand_chain])
Lrms=super_imposer.rms

DockQ=(float(fnat) + 1/(1+(irms/1.5)*(irms/1.5)) + 1/(1+(Lrms/8.5)*(Lrms/8.5)))/3

print 'Fnat ' + fnat
print 'iRMS ' + str(irms)
print 'LRMS ' + str(Lrms)
print 'DockQ ' + str(DockQ)
# Save the aligned version of 1UBQ.pdb
#io = Bio.PDB.PDBIO()
#io.set_structure(sample_structure) 
#io.save("1UBQ_aligned.pdb")
