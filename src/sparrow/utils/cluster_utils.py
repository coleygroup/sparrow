from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

def compute_fps(smiles: list, radius=3, length=2048):
    """ Computes Morgan fingerprints for a list of smiles"""
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=length)
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    fps = [mfpgen.GetCountFingerprint(m) for m in mols]
    return fps 

def cluster_fps(fps, cutoff):
    """ Clusters fingerprints using Butina clustering algorithm from RDKit 
    Adapted from the RDKit Cookbook (https://rdkit.readthedocs.io/en/latest/Cookbook.html) """
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True)
    return cs

def cluster_smiles(smiles, radius=2, length=1024, cutoff=0.7) -> list:
    """ 
    Clusters smiles representations of molecules. 
    Returns a list of cluster ids corresponding to the order of entries in 'smiles' input 
    """
    fps = compute_fps(smiles, radius=radius, length=length)
    cs = cluster_fps(fps, cutoff=cutoff)
    clusters = [list(cluster) for cluster in cs]

    return clusters

if __name__=='__main__':
    smis_list = [
        'NC(=O)N1CCN(c2cccc(C(F)(F)F)c2)CC1',
        'O=C(Nc1cccc(F)c1)N1CCC(c2ccccc2)CC1',
        'O=C(CCl)Nc1ccccc1N1CCCCC1',
        'CCC(=O)N1CCN(c2ccc(Cl)cc2N)CC1',
        'Nc1cc(Cl)ccc1N1CCN(C(=O)C(F)(F)F)CC1',
        'CC(=O)N1CCN(c2ccc(Cl)cc2N)CC1',
        'NCCCN1CCN(c2ccccc2Cl)CC1',
        'CNC(=O)C(C(C)C)N(C(CCCC#N)c1ccc(C#N)cc1)C(c1ccc(N2CCNCC2)cc1)c1cccc(Cl)c1',
        'CC(Cl)C(=O)NCC(c1ccccc1)N1CCN(c2ccccc2F)CC1',
        'CC(Cl)C(=O)NCC(c1ccccc1)N1CCN(c2ccc(C(F)(F)F)cc2F)CC1',
        'CC(Cl)C(=O)NCC(c1ccccc1)N1CCN(c2cc(F)ccc2F)CC1',
        'CC(Cl)C(=O)NCC(c1ccccc1)N1CCN(c2cccc(F)c2F)CC1',
        'CC(Cl)C(=O)NCC(c1ccccc1)N1CCN(c2c(F)cc(F)cc2F)CC1',
        'CN(C)C(=O)c1cccc(CNc2cc(Cl)ccc2N2CCC(CO)CC2)c1',
    ]
    cluster_smiles(smis_list)