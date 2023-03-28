from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdchem
import re
from rdkit.Chem import PeriodicTable
from rdkit.Chem import GetPeriodicTable
import random

from rdkit import RDLogger  
RDLogger.DisableLog('rdApp.*')  

class RXNMarkCenter:
    
    def __init__(self):
        self.mol_cache = {}
    
    def smi_tokenizer(self, smi):
        """
        Tokenize a SMILES molecule or reaction
        """
        pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\!|\$|\%[0-9]{2}|[0-9])"
        regex = re.compile(pattern)
        tokens = [token for token in regex.findall(smi)]
        assert smi == ''.join(tokens)
        return ' '.join(tokens)
    
    def GetAtomEnvironment(self, query_atom):
    
        atom_env = {}
        try:    atom_MapNum = query_atom.GetAtomMapNum()
        except: return atom_env #return empty env if the atom number does not exist in that molecule
        
        #For each bond:
        for bond in query_atom.GetBonds():
            atom_MapNum_1 = bond.GetBeginAtom().GetAtomMapNum()
            atom_MapNum_2 = bond.GetEndAtom().GetAtomMapNum()
            atom_connected_MapNum = 0
            if atom_MapNum == atom_MapNum_1:
                atom_connected_MapNum = atom_MapNum_2
            elif atom_MapNum == atom_MapNum_2:
                atom_connected_MapNum = atom_MapNum_1
        
            atom_env[atom_connected_MapNum] = bond.GetBondTypeAsDouble()
        
        return atom_env

    def GetAtomFromMapNum(self, molecules, query_atom):
        ''' Returns the atom  '''

        for molecule in molecules:
            for atom in molecule.GetAtoms():
                if query_atom == atom.GetAtomMapNum():
                    return atom
        return ''

    def ListMapNum(self, product_mol):
        ''' Returns the list of atom-mapped indices '''
        
        theList = []

        for atom in product_mol.GetAtoms():
            theList.append(atom.GetAtomMapNum())
        theList = list(set(theList))
        
        try: theList.remove(0)
        except: pass
        
        return theList
    
    def ListMapNum_count(self, product_mol):
        ''' Returns the count of atom-mapped atoms '''
        
        theList = []

        for atom in product_mol.GetAtoms():
            theList.append(atom.GetAtomMapNum())
        
        return theList.count(1)
    
    def mark_reaction_center_alt(self, rxn_SMILES_in, method=1):
        '''
        Takes as an input the "atom-mapped" traditionnal strategy marked reaction SMILES. 
        Outputs the alternative version of marking the atoms involved in the reaction.
        Having the advantage  of preserving the use of traditional atomic tokens.
        '''
        
        # Load:
        try:
            rxn_mapped = rdChemReactions.ReactionFromSmarts(rxn_SMILES_in, useSmiles=True)
            rxn_unmapped = rdChemReactions.ReactionFromSmarts(rxn_SMILES_in, useSmiles=True)
            rdChemReactions.SanitizeRxn(rxn_mapped)
            rdChemReactions.SanitizeRxn(rxn_unmapped)
        except:
            return ''
        
        # Get Marked Atoms Count:
        count_mapped_atoms = self.ListMapNum_count(rxn_mapped.GetProducts()[0])
        
        # Swap Num Radical Electrons:
        for product in rxn_mapped.GetProducts():
            for atom in product.GetAtoms():
                if atom.GetAtomMapNum() == 1:
                    if method == 1:     atom.SetAtomicNum(atom.GetAtomicNum()-1)
                    elif method == 2:   atom.SetNumRadicalElectrons(1)
                    elif method == 3:   pass
                    elif method == 4:   atom.SetAtomicNum(atom.GetAtomicNum()+1)
                    elif method == 5:
                        atom.SetAtomicNum(atom.GetAtomicNum()-1)
                        atom.SetNumRadicalElectrons(1)
                    elif method == 6:
                        atom.SetAtomicNum(atom.GetAtomicNum()+1)
                        atom.SetNumRadicalElectrons(1)
                    elif method == 7:
                        atom.InvertChirality()
                        atom.SetNumRadicalElectrons(1)
                    elif method == 8:
                        atom.InvertChirality()
                        atom.SetAtomicNum(atom.GetAtomicNum()+1)
                    elif method == 9:
                        atom.InvertChirality()
                        atom.SetAtomicNum(atom.GetAtomicNum()-1)
                    elif method == 10:
                        atom.InvertChirality()
                        atom.SetAtomicNum(atom.GetAtomicNum()-1)
                    elif method == 11:  atom.SetFormalCharge(0)
                    elif method == 12:
                        if random.randint(0, 1) == 1:   atom.SetIsAromatic(1)
                        if random.randint(0, 1) == 1:   atom.SetAtomicNum(atom.GetAtomicNum()-1)
                        if random.randint(0, 1) == 1:   atom.SetAtomicNum(atom.GetAtomicNum()+1)
                        if random.randint(0, 1) == 1:   atom.SetNumRadicalElectrons(1)
                        if random.randint(0, 1) == 1:   atom.InvertChirality()
                        if random.randint(0, 1) == 1:   atom.SetFormalCharge(0)
                        if random.randint(0, 1) == 1:   atom.SetFormalCharge(1)
                        if random.randint(0, 1) == 1:   atom.SetFormalCharge(-1)
                        
                    atom.SetAtomMapNum(0)
        
        try:
            # Get Smiles
            rdChemReactions.SanitizeRxn(rxn_mapped)
            smi_rxn_mapped = rdChemReactions.ReactionToSmiles(rxn_mapped, canonical=True)
            
            # Unmapped reference:
            rdChemReactions.RemoveMappingNumbersFromReactions(rxn_unmapped)
            rdChemReactions.SanitizeRxn(rxn_unmapped)
            smi_rxn_unmapped = rdChemReactions.ReactionToSmiles(rxn_unmapped, canonical=True)
        except:
            return ''
        
        # Special method 3:
        if method == 3:
            smi_rxn_mapped = rxn_SMILES_in

        # To list
        try:
            smi_rxn_mapped_list = self.smi_tokenizer(smi_rxn_mapped).split(' ')
            smi_rxn_unmapped_list = self.smi_tokenizer(smi_rxn_unmapped).split(' ')
        except: return ''
        
        #Check both lenghts
        if len(smi_rxn_mapped_list) != len(smi_rxn_unmapped_list): return ''
        
        # New Marking:
        MarkedSMILES_alt = ''
        count_mapped_atoms_changes = 0
        
        #Reconstruct the SMILES:
        for element in range (0, len(smi_rxn_mapped_list)):
            if smi_rxn_mapped_list[element] == smi_rxn_unmapped_list[element]:
                MarkedSMILES_alt += smi_rxn_unmapped_list[element] 
            else:
                MarkedSMILES_alt += smi_rxn_unmapped_list[element] + '!'
                count_mapped_atoms_changes += 1
        
        # Verify if the marked SMILES is the same as the "Unmapped" SMILES:
        if not MarkedSMILES_alt.replace('!', '') == smi_rxn_unmapped:
            return ''
        
        # Check all changes occurred:
        if count_mapped_atoms_changes != count_mapped_atoms:
            return ''
        
        return MarkedSMILES_alt
    
    def mark_reaction_center_alt_all(self, rxn_SMILES_in):
        
        # Try method 1 (least computationnaly hungry)
        MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=3)
        # For empty returns, try method 1
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=1)
        # For empty returns, try method 2
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=2)
        # For empty returns, try method 4
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=4)
        # For empty returns, try method 5
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=5)
        # For empty returns, try method 6
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=6)
        # For empty returns, try method 7
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=7)
        # For empty returns, try method 8
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=8)
        # For empty returns, try method 9
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=9)
        # For empty returns, try method 10
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=10)
        # For empty returns, try method 11
        if MarkedSMILES_alt == '':  MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=11)
        # For empty returns, try method 12
        if MarkedSMILES_alt == '':
            count = 0
            while MarkedSMILES_alt == '' and count < 5000:
                MarkedSMILES_alt = self.mark_reaction_center_alt(rxn_SMILES_in, method=12)
                count += 1
        
        return MarkedSMILES_alt
    
    def convert_smarts_to_smiles_alt_tagging(self, product_smarts_mapped_tag):
        '''
            Converts a product tagged by the atom-mapping SMARTS into the alternative tagging method:
        '''
        
        fake_rxn = self.mark_reaction_center_alt_all(rxn_SMILES_in='CCC>>' + product_smarts_mapped_tag)
        
        if '>>' in fake_rxn:
            return fake_rxn.split('>>')[1]
        else: return ''

    def TagMappedReactionCenter(self, MappedReaction, alternative_marking=False):
        
        try:    
            rxn = rdChemReactions.ReactionFromSmarts(MappedReaction, useSmiles=True)
            rdChemReactions.SanitizeRxn(rxn)
        except: return ''
        
        reactants = rxn.GetReactants()
        products = rxn.GetProducts()
        if len(products) > 1 : return '' 
        
        for molecule in reactants:
            try:    Chem.SanitizeMol(molecule)
            except: i=0

        for molecule in products:
            try:    Chem.SanitizeMol(molecule)
            except: i=0
            
        transformed_atoms = []
        
        for index in self.ListMapNum(products[0]):
            precursors_atom = self.GetAtomFromMapNum(reactants, index)
            products_atom = self.GetAtomFromMapNum(products, index)
        
            if self.GetAtomEnvironment(precursors_atom) != self.GetAtomEnvironment(products_atom):
                transformed_atoms.append(products_atom.GetIdx())
        
        # Clear atom mapping in Precursors:
        for molecule in reactants:
            for atom in molecule.GetAtoms():
                atom.SetAtomMapNum(0)
                
        # Clear atom mapping in Product (if unchanged):
        for molecule in products:
            for atom in molecule.GetAtoms():
                if not atom.GetIdx() in transformed_atoms:
                    atom.SetAtomMapNum(0)
                else:
                    atom.SetAtomMapNum(1)
        
        # Convert Mol to Smarts:
        rdChemReactions.SanitizeRxn(rxn)
        if not alternative_marking:
            return rdChemReactions.ReactionToSmiles(rxn, canonical=True)
        else:
            return self.mark_reaction_center_alt_all(rdChemReactions.ReactionToSmiles(rxn, canonical=True))
        
    def Mark_Random_Atoms(self, mol_SMILES, mark_count=2, neighbors=True, tokenized=False):
        '''
            Randomly mark an input SMILES single molecule.
            mark_count  ==> amount of atoms to be marked (from 1 to 3), default=2.
            neighbors   ==> whether or not to only mark multiple atoms that are next to each other, default=True.
            Returns a list of all possible random marking according to input parameters.
        '''
        
        New_SMILES = []
        
        # Load Molecule:
        try:    mol = Chem.MolFromSmarts(mol_SMILES)
        except: return []
        if mol is None: return []
        
        # PROCESS 
        if neighbors == False:
            if mark_count == 0: return [mol_SMILES]
            elif mark_count == 1:
                for first_atom in mol.GetAtoms():
                    first_atom.SetAtomMapNum(1)
                    New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                    first_atom.SetAtomMapNum(0)
            elif mark_count == 2:
                for first_atom in mol.GetAtoms():
                    first_atom.SetAtomMapNum(1)
                    for second_atom in mol.GetAtoms():
                        if second_atom.GetIdx() != first_atom.GetIdx():
                            second_atom.SetAtomMapNum(1)
                            New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                            second_atom.SetAtomMapNum(0)
                    first_atom.SetAtomMapNum(0)
            elif mark_count == 3:
                for first_atom in mol.GetAtoms():
                    first_atom.SetAtomMapNum(1)
                    for second_atom in mol.GetAtoms():
                        if second_atom.GetIdx() != first_atom.GetIdx():
                            second_atom.SetAtomMapNum(1)
                            for third_atom in mol.GetAtoms():
                                if third_atom.GetIdx() != second_atom.GetIdx() and third_atom.GetIdx() != first_atom.GetIdx():
                                    third_atom.SetAtomMapNum(1)
                                    New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                                    third_atom.SetAtomMapNum(0)
                            second_atom.SetAtomMapNum(0)
                    first_atom.SetAtomMapNum(0)
            else: return []
        
        elif neighbors == True:
            if mark_count == 0: return [mol_SMILES]
            elif mark_count == 1:
                for each_atom in mol.GetAtoms():
                    each_atom.SetAtomMapNum(1)
                    New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                    each_atom.SetAtomMapNum(0)
            elif mark_count == 2:
                for each_atom in mol.GetAtoms():
                    for each_bond in each_atom.GetBonds():
                        each_bond.GetBeginAtom().SetAtomMapNum(1)
                        each_bond.GetEndAtom().SetAtomMapNum(1)
                        New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                        each_bond.GetBeginAtom().SetAtomMapNum(0)
                        each_bond.GetEndAtom().SetAtomMapNum(0)
            elif mark_count == 3:
                for first_atom in mol.GetAtoms():
                    for each_bond in first_atom.GetBonds():
                        each_bond.GetBeginAtom().SetAtomMapNum(1)
                        each_bond.GetEndAtom().SetAtomMapNum(1)
                        
                        if first_atom.GetIdx() == each_bond.GetBeginAtom().GetIdx():
                            for each_bond2 in each_bond.GetEndAtom().GetBonds():
                                if first_atom.GetIdx() != each_bond2.GetBeginAtom().GetIdx() and first_atom.GetIdx() != each_bond2.GetEndAtom().GetIdx():
                                    each_bond2.GetBeginAtom().SetAtomMapNum(1)
                                    each_bond2.GetEndAtom().SetAtomMapNum(1)
                                    New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                                    each_bond2.GetBeginAtom().SetAtomMapNum(0)
                                    each_bond2.GetEndAtom().SetAtomMapNum(0)
                            
                        elif first_atom.GetIdx() == each_bond.GetEndAtom().GetIdx():
                            for each_bond2 in each_bond.GetBeginAtom().GetBonds():
                                if first_atom.GetIdx() != each_bond2.GetBeginAtom().GetIdx() and first_atom.GetIdx() != each_bond2.GetEndAtom().GetIdx():
                                    each_bond2.GetBeginAtom().SetAtomMapNum(1)
                                    each_bond2.GetEndAtom().SetAtomMapNum(1)
                                    New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                                    each_bond2.GetBeginAtom().SetAtomMapNum(0)
                                    each_bond2.GetEndAtom().SetAtomMapNum(0)
                            
                        each_bond.GetBeginAtom().SetAtomMapNum(0)
                        each_bond.GetEndAtom().SetAtomMapNum(0)
            else: return []
            
        New_SMILES = list(set(New_SMILES))
        
        # Convert Marking:
        New_SMILESs = []
        if tokenized:   
            for element in New_SMILES: New_SMILESs.append(self.smi_tokenizer(self.mark_reaction_center_alt_all(">>" + element).replace('>>', '')))
        else:           
            for element in New_SMILES: New_SMILESs.append(self.mark_reaction_center_alt_all(">>" + element).replace('>>', ''))
        
        return New_SMILESs
    
    def Mark_matching_substructures_OLD(self, mol_SMILES, list_substructures, tokenized=False):
        '''
            Really not the best code ever, we are getting a lot of extra possibilities, but that does not hurt much 
        '''
        
        try:        mol = Chem.MolFromSmarts(mol_SMILES)
        except:     return []
        
        New_SMILES = []

        for substructure in list_substructures:
            try:    
                sub = Chem.MolFromSmarts(substructure)
                matches = mol.GetSubstructMatches(sub)
            except: continue
            
            if len(matches)> 0:
                for match_set in matches:
                    # Mark all atoms:
                    for atom in mol.GetAtoms():
                        if atom.GetIdx() in match_set:
                            atom.SetAtomMapNum(1)
                    # Save SMILES:
                    New_SMILES.append(Chem.MolToSmiles(mol, canonical=True))
                    
                    # Unmark all atoms:
                    for atom in mol.GetAtoms():
                        if atom.GetIdx() in match_set:
                            atom.SetAtomMapNum(0)
                            
        New_SMILES = list(set(New_SMILES))
                            
        # Convert Marking:
        New_SMILESs = []
        if tokenized:   
            for element in New_SMILES: New_SMILESs.append(self.smi_tokenizer(self.mark_reaction_center_alt_all(">>" + element).replace('>>', '')))
        else:           
            for element in New_SMILES: New_SMILESs.append(self.mark_reaction_center_alt_all(">>" + element).replace('>>', ''))
            
        return New_SMILESs

    def Load_mol_in_cache(self, mol_SMILES):
        '''
            mol_SMILES: SMILES representation of the molecule we would like to load in the cache
            Load the molecule in the cache and return it
        '''
        mol = Chem.MolFromSmarts(mol_SMILES)
        Chem.SanitizeMol(mol)
        self.mol_cache[mol_SMILES] = mol
        return mol
    
    def Mark_matching_substructures(self, mol_SMILES, list_conditionnal_substructures_tags, tokenized=False):
        '''
            mol_SMILES: SMILES representation of the molecule we would like to tag
            list_conditionnal_substructures_tags: list of substructures and tags to be applied on the molecule if the substructure is present (list of list)
            tokenized: if True, the output will be tokenized
            
            Tags substructure[1] on the mol_SMILES under the condition that substructure[0] matches mol_SMILES. 
            This makes substructure[0] to me the conditionnal substructure to be present in order to tag the molecule. It limits the amount of tags. 
            You should provide substructure[0] with enough radious according to your needs.
            Really not the best code ever, we are getting a lot of extra no-matching tags , but that does not hurt much. 
        '''
        
        try:        mol1 = Chem.MolFromSmiles(mol_SMILES)
        except:     return []
        
        New_SMILES = []

        for substructure in list_conditionnal_substructures_tags:
            #Load both substructures and tags:
            try:
                conditionnal_substructure = substructure[0]
                substructures_tag =         substructure[1]
                
                # Load in cache if not already done:
                if conditionnal_substructure in self.mol_cache: mol2 = self.mol_cache[conditionnal_substructure]
                else: mol2 = self.Load_mol_in_cache(conditionnal_substructure)
                if substructures_tag in self.mol_cache: mol3 = self.mol_cache[substructures_tag]
                else: mol3 = self.Load_mol_in_cache(substructures_tag)
            except: continue

            try:    matches = mol1.GetSubstructMatches(mol2)
            except: continue
        
            if any(matches):
                try: matches2 = mol1.GetSubstructMatches(mol3)
                except: continue
                for match_set in matches:
                    for match_set2 in matches2:
                        # Only continue this set if all atoms from the tagging structure (matches2) are present in the matches:
                        if all(i in match_set for i in match_set2):
                            # Mark all atoms:
                            for atom in mol1.GetAtoms():
                                if atom.GetIdx() in match_set2:
                                    atom.SetAtomMapNum(1)
                            # Save SMILES:
                            New_SMILES.append(Chem.MolToSmiles(mol1, canonical=True))

                            # Unmark all atoms:
                            for atom in mol1.GetAtoms():
                                atom.SetAtomMapNum(0)
        
        New_SMILES = list(set(New_SMILES))
                            
        # Convert Marking:
        New_SMILESs = []
        if tokenized:   
            for element in New_SMILES: New_SMILESs.append(self.smi_tokenizer(self.mark_reaction_center_alt_all(">>" + element).replace('>>', '')))
        else:           
            for element in New_SMILES: New_SMILESs.append(self.mark_reaction_center_alt_all(">>" + element).replace('>>', ''))
            
        return New_SMILESs
        