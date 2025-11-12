from itertools import combinations


from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import IPythonConsole
from rdkit.DataStructs import FingerprintSimilarity

# insert SMILES
cocaine = ('CN1[C@H]2CC[C@@H]1[C@@H](C(OC)=O)[C@@H](OC(C3=CC=CC=C3)=O)C2') 
heroin = ('CC(OC1=C(O[C@@H]2[C@]34CCN(C)[C@@H]([C@@H]4C=C[C@@H]2OC(C)=O)C5)C3=C5C=C1)=O') 
fentanyl = ('O=C(CC)N(C1CCN(CC1)CCc2ccccc2)c3ccccc3')
meth = ('CNC(C)Cc1ccccc1')
THC= ('CCCCCc1cc(c2c(c1)OC([C@H]3[C@H]2C=C(CC3)C)(C)C)O')

smiles_list=[cocaine,heroin,fentanyl, meth, THC]

#added names so that it prints names instead of 0 1 2 from range(len())
names = ["cocaine", "heroin", "fentanyl", "meth", "THC"]

def tanimoto(smiles_list):
    #turn smiles into mols
    mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

    #define fingerprint generator
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)

    fp_list = [mfpgen.GetFingerprint(mol) for mol in mols]

    #something similar to I would write in Matlab
    #tanimotosim = [FingerprintSimilarity(i,i+1) for i-1 in fp_list]

    #something that would work in python but it still only compares adjacent smiles only
    #tanimotosim = [FingerprintSimilarity(fp_list[i],fp_list[i+1]) for i in range(len(fp_list)-1)]

    #a list of tuples i,j and tanimotosim
    tanimotomat=[]
    
    #new version with combinations
    for i,j in combinations(range(len(fp_list)),2):
        tanimotosim = DataStructs.FingerprintSimilarity(fp_list[i],fp_list[j])
        
        #generate the tuple for tuple unpacking below
        tanimotomat.append((i,j,tanimotosim))
    return tanimotomat


simlist=tanimoto(smiles_list)

#tuple unpacking to insert names
for i, j, tanimotosim in simlist:
    print(f"{names[i]} vs {names[j]}: {tanimotosim:.5f}")
