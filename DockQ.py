#!/usr/bin/env python

import Bio.PDB
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
import sys
import os
import re
import tempfile
import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer
from math import sqrt
from argparse import ArgumentParser
import itertools
import commands

def parse_fnat(fnat_out):
    fnat=-1;
    nat_correct=-1
    nat_total=-1
    fnonnat=-1
    nonnat_count=-1
    model_total=-1
    inter=[]
    for line in fnat_out.split("\n"):
#        print line
        line=line.rstrip('\n')
        match=re.search(r'NATIVE: (\d+)(\w) (\d+)(\w)',line)
        if(re.search(r'^Fnat',line)):
            list=line.split(' ')
            fnat=float(list[3])
            nat_correct=int(list[1])
            nat_total=int(list[2])
        elif(re.search(r'^Fnonnat',line)):
            list=line.split(' ')
            fnonnat=float(list[3])
            nonnat_count=int(list[1])
            model_total=int(list[2])
        elif(match):
            #print line
            res1=match.group(1)
            chain1=match.group(2)
            res2=match.group(3)
            chain2=match.group(4)
           # print res1 + ' ' + chain1 + ' ' + res2 + ' ' + chain2
            inter.append(res1 + chain1)
            inter.append(res2 + chain2)
    return (fnat,nat_correct,nat_total,fnonnat,nonnat_count,model_total,inter)

def capri_class(fnat,iRMS,LRMS):

    if(fnat < 0.1 or (LRMS > 10.0 and iRMS > 4.0)):
        return 'Incorrect'
    elif((fnat >= 0.1 and fnat < 0.3) and (LRMS <= 10.0 or iRMS <= 4.0) or (fnat >= 0.3 and LRMS > 5.0 and iRMS > 2.0)):
        return 'Acceptable'
    elif((fnat >= 0.3 and fnat < 0.5) and (LRMS <= 5.0 or iRMS <= 2.0) or (fnat >= 0.5 and LRMS > 1.0 and iRMS > 1.0)):
        return 'Medium'
    elif(fnat >= 0.5 and (LRMS <= 1.0 or iRMS <= 1.0)):
        return 'High'
    else:
        return 'Undef'

def capri_class_DockQ(DockQ):

    (c1,c2,c3)=(0.23,0.49,0.80)
    if(DockQ < c1):
        return 'Incorrect'
    elif(DockQ >= c1 and DockQ < c2):
        return 'Acceptable'
    elif(DockQ >= c2 and DockQ < c3):
        return 'Medium'
    elif(DockQ >= c3):
        return 'High'
    else:
        return 'Undef'


def calc_DockQ(model,native,use_CA_only=False):
    
    exec_path=os.path.dirname(os.path.abspath(sys.argv[0]))    
    atom_for_sup=['CA','C','N','O']
    if(use_CA_only):
        atom_for_sup=['CA']

    cmd_fnat=exec_path + '/fnat ' + model + ' ' + native + ' 5'
    #cmd_interface=exec_path + '/fnat ' + model + ' ' + native + ' 10 backbone'
    cmd_interface=exec_path + '/fnat ' + model + ' ' + native + ' 10'


    #fnat_out = os.popen(cmd_fnat).readlines()
    fnat_out = commands.getoutput(cmd_fnat)
#    sys.exit()
    (fnat,nat_correct,nat_total,fnonnat,nonnat_count,model_total,interface5A)=parse_fnat(fnat_out)
    assert fnat!=-1, "Error running cmd: %s\n" % (cmd_fnat)
#    inter_out = os.popen(cmd_interface).readlines()
    inter_out = commands.getoutput(cmd_interface)
    (fnat_bb,nat_correct_bb,nat_total_bb,fnonnat_bb,nonnat_count_bb,model_total_bb,interface)=parse_fnat(inter_out)
    assert fnat_bb!=-1, "Error running cmd: %s\n" % (cmd_interface)

    #print fnat
    #Use same interface as for fnat for iRMS
    #interface=interface5A


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


    #find atoms common in both sample and native
    atoms_def_sample=[]
    atoms_def_in_both=[]
    #first read in sample
    for sample_chain in sample_model:
#        print sample_chain
        chain=sample_chain.id
#        print chain
        for sample_res in sample_chain:
           # print sample_res
            if sample_res.get_id()[0] != ' ': #Skip hetatm.
                continue
            resname=sample_res.get_id()[1]
            key=str(resname) + chain
            for a in atom_for_sup:
                atom_key=key + '.' + a
                if a in sample_res:
                    if atom_key in atoms_def_sample:
                        print atom_key + ' already added (MODEL)!!!'
                    atoms_def_sample.append(atom_key)

    #then read in native also present in sample
    for ref_chain in ref_model:
        chain=ref_chain.id
        for ref_res in ref_chain:
            #print ref_res
            if ref_res.get_id()[0] != ' ': #Skip hetatm.
#                print ref_res.get_id()
                continue
            resname=ref_res.get_id()[1]
            key=str(resname) + chain
            for a in atom_for_sup:
                atom_key=key + '.' + a
                if a in ref_res and atom_key in atoms_def_sample:
                    if atom_key in atoms_def_in_both:
                        print atom_key + ' already added (Native)!!!' 
                    atoms_def_in_both.append(atom_key)


#    print atoms_def_in_both
    for sample_chain in sample_model:
        chain=sample_chain.id
        if chain not in chain_res.keys():
            chain_res[chain]=[]
        for sample_res in sample_chain:
            if sample_res.get_id()[0] != ' ': #Skip hetatm.
                continue
            resname=sample_res.get_id()[1]
            key=str(resname) + chain
            chain_res[chain].append(key)
            if key in interface:
                for a in atom_for_sup:
                    atom_key=key + '.' + a
                    if a in sample_res and atom_key in atoms_def_in_both:
                        sample_atoms.append(sample_res[a])
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
            if ref_res.get_id()[0] != ' ': #Skip hetatm.
                continue
            resname=ref_res.get_id()[1]
            key=str(resname) + chain

            #print ref_res
            #      print key
            # print chain_res.values()
            if key in chain_res[chain]: # if key is present in sample
                #print key
                for a in atom_for_sup:
                    atom_key=key + '.' + a
                    if a in ref_res and atom_key in atoms_def_in_both:
                        chain_ref[chain].append(ref_res[a])
                        common_residues.append(key)
                      #chain_sample.append((ref_res['CA'])
            if key in common_interface:
              # Check if residue number ( .get_id() ) is in the list
              # Append CA atom to list
                #print key  
                for a in atom_for_sup:
                    atom_key=key + '.' + a
                    #print atom_key
                    if a in ref_res and atom_key in atoms_def_in_both:
                        ref_atoms.append(ref_res[a])



    #get the ones that are present in native        
    chain_sample={}
    for sample_chain in sample_model:
        chain=sample_chain.id
        if chain not in chain_sample.keys():
            chain_sample[chain]=[]
        for sample_res in sample_chain:
            if sample_res.get_id()[0] != ' ': #Skip hetatm.
                continue
            resname=sample_res.get_id()[1]
            key=str(resname) + chain
            if key in common_residues:
                for a in atom_for_sup:
                    atom_key=key + '.' + a
                    if a in sample_res and atom_key in atoms_def_in_both:
                        chain_sample[chain].append(sample_res[a])

        #if key in common_residues:
        #     print key  
        #sample_atoms.append(sample_res['CA'])
        #common_interface.append(key)


    assert len(ref_atoms)!=0, "length of native is zero"
    assert len(sample_atoms)!=0, "length of model is zero"
    assert len(ref_atoms)==len(sample_atoms), "Different number of atoms in native and model %d %d\n" % (len(ref_atoms),len(sample_atoms))

    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_model.get_atoms())

    # Print RMSD:
    irms=super_imposer.rms

    (chain1,chain2)=chain_sample.keys()

    ligand_chain=chain1
    receptor_chain=chain2
    len1=len(chain_res[chain1])
    len2=len(chain_res[chain2])

    assert len1!=0, "%s chain has zero length!\n" % chain1
    assert len2!=0, "%s chain has zero length!\n" % chain2

    class1='ligand'
    class2='receptor'
    if(len(chain_sample[chain1]) > len(chain_sample[chain2])):
        receptor_chain=chain1
        ligand_chain=chain2
        class1='receptor'
        class2='ligand'



    #print len1
    #print len2
    #print chain_sample.keys()

    #Set to align on receptor
    assert len(chain_ref[receptor_chain])==len(chain_sample[receptor_chain]), "Different number of atoms in native and model receptor (chain %c) %d %d\n" % (receptor_chain,len(chain_ref[receptor_chain]),len(chain_sample[receptor_chain]))

    super_imposer.set_atoms(chain_ref[receptor_chain], chain_sample[receptor_chain])
    super_imposer.apply(sample_model.get_atoms())
    receptor_chain_rms=super_imposer.rms
    #print receptor_chain_rms
    #print dir(super_imposer)
    #print chain1_rms

    #Grep out the transformed ligand coords

    #print ligand_chain

    #print chain_ref[ligand_chain]
    #print chain_sample[ligand_chain]
    #l1=len(chain_ref[ligand_chain])
    #l2=len(chain_sample[ligand_chain])




    assert len(chain_ref[ligand_chain])!=0 or len(chain_sample[ligand_chain])!=0, "Zero number of equivalent atoms in native and model ligand (chain %s) %d %d.\nCheck that the residue numbers in model and native is consistent\n" % (ligand_chain,len(chain_ref[ligand_chain]),len(chain_sample[ligand_chain]))


    assert len(chain_ref[ligand_chain])==len(chain_sample[ligand_chain]), "Different number of atoms in native and model ligand (chain %c) %d %d\n" % (ligand_chain,len(chain_ref[ligand_chain]),len(chain_sample[ligand_chain]))

    coord1=np.array([atom.coord for atom in chain_ref[ligand_chain]])
    coord2=np.array([atom.coord for atom in chain_sample[ligand_chain]])

    #coord1=np.array([atom.coord for atom in chain_ref[receptor_chain]])
    #coord2=np.array([atom.coord for atom in chain_sample[receptor_chain]])

    #print len(coord1)
    #print len(coord2)

    sup=SVDSuperimposer()
    Lrms = sup._rms(coord1,coord2) #using the private _rms function which does not superimpose


    #super_imposer.set_atoms(chain_ref[ligand_chain], chain_sample[ligand_chain])
    #super_imposer.apply(sample_model.get_atoms())
    #coord1=np.array([atom.coord for atom in chain_ref[receptor_chain]])
    #coord2=np.array([atom.coord for atom in chain_sample[receptor_chain]])
    #Rrms= sup._rms(coord1,coord2)
    #should give same result as above line
    #diff = coord1-coord2
    #l = len(diff) #number of atoms
    #from math import sqrt
    #print sqrt(sum(sum(diff*diff))/l)
    #print np.sqrt(np.sum(diff**2)/l)
    DockQ=(float(fnat) + 1/(1+(irms/1.5)*(irms/1.5)) + 1/(1+(Lrms/8.5)*(Lrms/8.5)))/3
    dict={}
    dict['DockQ']=DockQ
    dict['irms']=irms
    dict['Lrms']=Lrms
    dict['fnat']=fnat
    dict['nat_correct']=nat_correct
    dict['nat_total']=nat_total

    dict['fnonnat']=fnonnat
    dict['nonnat_count']=nonnat_count
    dict['model_total']=model_total
    
    dict['chain1']=chain1
    dict['chain2']=chain2
    dict['len1']=len1
    dict['len2']=len2
    dict['class1']=class1
    dict['class2']=class2
    
    return dict

def get_pdb_chains(pdb):
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    pdb_struct = pdb_parser.get_structure("reference", pdb)[0]
    chain=[]
    for c in pdb_struct:
        chain.append(c.id)
    return chain
#ATOM   3312  CA
#ATOM   3315  CB  ALA H 126     -21.802  31.674  73.597  1.00 58.05           C  

def make_two_chain_pdb(pdb,group1,group2): #renumber from 1
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)
    pdb_struct = pdb_parser.get_structure("reference", pdb)[0]
    for c in pdb_struct:
        if c.id in group1:
            c.id='A'
        if c.id in group2:
            c.id='B'
    (code,outfile)=tempfile.mkstemp()
    io=Bio.PDB.PDBIO()
    io.set_structure(pdb_struct)
    io.save(outfile)
    exec_path=os.path.dirname(os.path.abspath(sys.argv[0]))    
    cmd=exec_path + '/scripts/renumber_pdb.pl ' + outfile
    os.system(cmd)
    os.remove(outfile)
    return outfile +'.renum'

def change_chain(pdb_string,chain):
    new_str=[];
    for line in pdb_string:
        s=list(line);
        s[21]=chain;
        new_str.append("".join(s));
    return "\n".join(new_str);

def make_two_chain_pdb_perm(pdb,group1,group2): #not yet ready
    pdb_chains={}
    f=open(pdb);
    for line in f.readlines():
        if line[0:4] == "ATOM":
            #       s=list(line);
            #print line
            chain=line[21]
            atom=line[13:16]
            resnum=int(line[22:26])
            #        print atom + ':' + str(resnum) +':'
            if chain not in pdb_chains:
                pdb_chains[chain]=[]
            pdb_chains[chain].append(line)
 #       print line
 #       s[21]='B'
 #       print "".join(s)
#        print chain


    f.close()
    #sys.exit()
    (code,outfile)=tempfile.mkstemp()
    f=open(outfile,'w')
    for c in group1:
     #   print pdb_chains[c]
        f.write(change_chain(pdb_chains[c],"A"))
    f.write("TER\n");
    for c in group2:
        f.write(change_chain(pdb_chains[c],"B"))
    f.close();
    #print outfile
    exec_path=os.path.dirname(os.path.abspath(sys.argv[0]))    
    cmd=exec_path + '/scripts/renumber_pdb.pl ' + outfile
    os.system(cmd)
    os.remove(outfile)
    return outfile +'.renum'
    
def main():

    parser=ArgumentParser(description="DockQ - Quality measure for protein-protein docking models")
    parser.add_argument('model',metavar='<model>',type=str,nargs=1,help='path to model file')
    parser.add_argument('native',metavar='<native>',type=str,nargs=1,help='path to native file')
    parser.add_argument('-short',default=False,action='store_true',help='short output')
    parser.add_argument('-verbose',default=False,action='store_true',help='talk a lot!')
    parser.add_argument('-useCA',default=False,action='store_true',help='use CA instead of backbone')
    parser.add_argument('-skip_check',default=False,action='store_true',help='skip initial check fo speed up on two chain examples')
    parser.add_argument('-perm1',default=False,action='store_true',help='use all chain1 permutations to find maximum DockQ (number of comparisons is n! = 24, if combined with -perm2 there will be n!*m! combinations')
    parser.add_argument('-perm2',default=False,action='store_true',help='use all chain2 permutations to find maximum DockQ (number of comparisons is n! = 24, if combined with -perm1 there will be n!*m! combinations')
#    parser.add_argument('-comb',default=False,action='store_true',help='use all cyclicchain permutations to find maximum DockQ (number of comparisons is n!*m! = 24*24 = 576 for two tetramers interacting')
    parser.add_argument('-chain1',metavar='chain1', type=str,nargs='+', help='pdb chain order to group together partner 1')
    parser.add_argument('-chain2',metavar='chain2', type=str,nargs='+', help='pdb chain order to group together partner 2 (complement to partner 1 if undef)')
    parser.add_argument('-native_chain1',metavar='native_chain1', type=str,nargs='+', help='pdb chain order to group together from native partner 1')
    parser.add_argument('-native_chain2',metavar='native_chain2', type=str,nargs='+', help='pdb chain order to group together from native partner 2 (complement to partner 1 if undef)')


    args = parser.parse_args()
    #bio_ver=1.64
    bio_ver=1.61
    if(float(Bio.__version__) < bio_ver):
        print "Biopython version (%s) is too old need at least >=%.2f" % (Bio.__version__,bio_ver)
        sys.exit()

#    if(len(sys.argv)!=3):
#        print "Usage: ./Dock.py <model> <native>"
#        sys.exit()

#    print args
#    print args.model[0]
#    sys.exit()
#    model=sys.argv[1]
#    native=sys.argv[2]

    exec_path=os.path.dirname(os.path.abspath(sys.argv[0]))    
    model=args.model[0]
    model_in=model
    native=args.native[0]
    native_in=native
    use_CA_only=args.useCA

    model_chains=[]
    native_chains=[]
    if(not args.skip_check):
        model_chains=get_pdb_chains(model)
        native_chains=get_pdb_chains(native)
    files_to_clean=[]

#    print native_chains
    if((len(model_chains) > 2 or len(native_chains) > 2) and
       (args.chain1 == None and args.native_chain1 == None)):
        print "Multi-chain model need sets of chains to group\nuse -chain1 and -native_chain1"
        sys.exit()
    if not args.skip_check and (len(model_chains) < 2 or len(native_chains)< 2):
        print "Need at least two chains in the two inputs\n";
        sys.exit()
        
    if len(model_chains) > 2 or len(native_chains)> 2:
        group1=model_chains[0]
        group2=model_chains[1]
        nat_group1=native_chains[0]
        nat_group2=native_chains[1]
        if(args.chain1 != None):
            group1=args.chain1
            nat_group1=group1
            if(args.chain2 != None):
                group2=args.chain2
            else:
                #will use the complement from group1
                group2=[]
                for c in model_chains:
                    if c not in group1:
                        group2.append(c)
            nat_group1=group1
            nat_group2=group2
            

        if(args.native_chain1 != None):
            nat_group1=args.native_chain1
            if(args.native_chain2 != None):
                nat_group2=args.native_chain2
            else:
                #will use the complement from group1
                nat_group2=[]
                for c in native_chains:
                    if c not in nat_group1:
                        nat_group2.append(c)
                        
        if(args.chain1 == None):
            group1=nat_group1
            group2=nat_group2

        #print group1
        #print group2

        #print "native"
        #print nat_group1
        #print nat_group2
        print str(group1) + ' -> ' + str(nat_group1)
        print str(group2) + ' -> ' + str(nat_group2)
        native=make_two_chain_pdb_perm(native,nat_group1,nat_group2)
        files_to_clean.append(native)
        pe=0
        if args.perm1 or args.perm2:
            best_DockQ=-1;
            best_g1=[]
            best_g2=[]
            
            iter_perm1=itertools.combinations(group1,len(group1))
            iter_perm2=itertools.combinations(group2,len(group2))
            if args.perm1:
                iter_perm1=itertools.permutations(group1)
            if args.perm2:
                iter_perm2=itertools.permutations(group2)
               
            combos1=[];
            combos2=[];
            for g1 in iter_perm1:#_temp:
                combos1.append(g1)
            for g2 in iter_perm2:
                combos2.append(g2)

            for g1 in combos1:
                for g2 in combos2:
                    pe=pe+1
                   # print str(g1)+' '+str(g2)
#            print pe
#            print group1
#            print group2
            pe_tot=pe
            pe=1
            #sys.exit()
            for g1 in combos1:
                for g2 in combos2:
                #g2=group2    
                    model=make_two_chain_pdb_perm(model_in,g1,g2)
                    test_dict=calc_DockQ(model,native,use_CA_only)
                    os.remove(model)
#                    if args.verbose:
                    print str(pe)+'/'+str(pe_tot) + ' ' + str(g1) + ' ' + str(g2) + ' ' + str(test_dict['DockQ'])
                    if(test_dict['DockQ'] > best_DockQ):
                        best_DockQ=test_dict['DockQ'];
                        dict=test_dict
                        best_g1=g1
                        best_g2=g2
 #                       if args.verbose:
                        print "Current best " + str(best_DockQ)
                    pe=pe+1
            print 'Best score ( ' + str(best_DockQ) +' ) found for ' + str(best_g1) + ' ' + str(best_g2)
        else:
            model=make_two_chain_pdb_perm(model,group1,group2)
            dict=calc_DockQ(model,native,use_CA_only)
        
            os.system('cp ' + native + ' native_multichain.pdb')
            os.system('cp ' + model + ' model_multichain.pdb')

            files_to_clean.append(model)
   #    sys.exit()
     
 #   print native
 #   print model
    else:
        dict=calc_DockQ(model,native,use_CA_only)
    
    irms=dict['irms']
    Lrms=dict['Lrms']
    fnat=dict['fnat']
    DockQ=dict['DockQ']
    
    
    if(args.short):
        print DockQ
    else:
        print '***********************************************************'
        print '*                       DockQ                             *'
        print '*   Scoring function for protein-protein docking models   *'
        print '*   Statistics on CAPRI data:                             *'
        print '*    0.00 <= DockQ <  0.23 - Incorrect                    *'
        print '*    0.23 <= DockQ <  0.49 - Acceptable quality           *'
        print '*    0.49 <= DockQ <  0.80 - Medium quality               *'
        print '*            DockQ >= 0.80 - High quality                 *'  
        print '*   Reference: Sankar Basu and Bjorn Wallner, DockQ:...   *'
        print '*   For comments, please email: bjornw@ifm.liu.se         *'
        print '***********************************************************\n'
        print 'Number of equivalent residues in chain ' + dict['chain1'] + ' ' + str(dict['len1']) + ' (' + dict['class1'] + ')'
        print 'Number of equivalent residues in chain ' + dict['chain2'] + ' ' + str(dict['len2']) + ' (' + dict['class2'] + ')'
        print("Fnat %.3f %d correct of %d native contacts" % (dict['fnat'],dict['nat_correct'],dict['nat_total']))
        print("Fnonnat %.3f %d non-native of %d model contacts" % (dict['fnonnat'],dict['nonnat_count'],dict['model_total']))
        print("iRMS %.3f" % irms)
        print("LRMS %.3f" % Lrms)
        print 'CAPRI ' + capri_class(fnat,irms,Lrms)
        print 'DockQ_CAPRI ' + capri_class_DockQ(DockQ)
        print("DockQ %.3f" % DockQ)

    for f in files_to_clean:
        os.remove(f)

if __name__ == '__main__':
  main()    
