# Tag Reaction Center
 Converts an atom-mapped reaction into a Reaction SMILES where the product reaction center is marked

# Environment:

``` console
conda create -n rxnmarkcenter python=3.6 -y
conda activate rxnmarkcenter 
conda install -c rdkit rdkit=2020.03.3 -y
git clone https://github.com/DavKre/TagReactionCenter.git 
cd MarkReactionCenter
pip install -e .
```


# Example of use

``` python
from rxnmarkcenter import RXNMarkCenter
rxn_mark_center = RXNMarkCenter()

mapped_rxn = 'CN(C)C=O.F[c:5]1[n:6][cH:7][cH:8][cH:9][c:10]1[F:11].O=C([O-])[O-].[CH3:1][CH:2]([CH3:3])[SH:4].[K+].[K+]>>[CH3:1][CH:2]([CH3:3])[S:4][c:5]1[n:6][cH:7][cH:8][cH:9][c:10]1[F:11]'

result_1 = rxn_mark_center.TagMappedReactionCenter(mapped_rxn)
result_2 = rxn_mark_center.TagMappedReactionCenter(mapped_rxn, alternative_marking=True)
```


