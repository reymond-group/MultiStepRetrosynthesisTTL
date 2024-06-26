from __future__ import division, unicode_literals

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*') 
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import pandas as pd
import os
from random import sample
import re

import datetime
import subprocess
import onmt.bin.translate as trsl 
from typing import List, Tuple

from ttlretro.rxnmarkcenter import RXNMarkCenter
from ttlretro.utils.confidence_normalization_enzr_full import normalize_enzr_conf


class SingleStepRetrosynthesis:
    '''
        Class for Single Step Retrosynthesis.
        The input is passed over the Retrosynthesis Transformer model, followed by reagent prediction, then predictions are validated over the Forward Model.
    '''

    def __init__(
        self, 
        log_folder = '', 
        log_time_stamp = '', 
        list_substructures = [['', '']], 
        list_substructures_ENZR = [['', '']], 
        USPTO_AutoTag_path = '',
        USPTO_T1_path = '', 
        USPTO_T2_path = '',
        USPTO_T3_path = '',
        ENZR_AutoTag_path = '',
        ENZR_T1_path = '', 
        ENZR_T2_path = '',
        ENZR_T3_path = '',
        ENZR_confidence_threshold = 0.0, 
        tmp_file_path = 'tmp/',
        ):
        
        self.USPTO_AutoTag_path = USPTO_AutoTag_path
        self.USPTO_T1_path = USPTO_T1_path
        self.USPTO_T2_path = USPTO_T2_path
        self.USPTO_T3_path = USPTO_T3_path
        self.ENZR_AutoTag_path = ENZR_AutoTag_path
        self.ENZR_T1_path = ENZR_T1_path
        self.ENZR_T2_path = ENZR_T2_path
        self.ENZR_T3_path = ENZR_T3_path
        self.ENZR_confidence_threshold = ENZR_confidence_threshold
        
        #      Custom Model:
        self.Custom_Model = 'Custom_Model.pt'  # for round-trip accuracy testing
        self.tmp_file_path = tmp_file_path
        
        os.environ['MKL_THREADING_LAYER'] = 'GNU'
        self.log_folder = log_folder
        self.log_time_stamp = log_time_stamp

        self.list_substructures = list_substructures
        self.list_substructures_ENZR = list_substructures_ENZR
        self.rxn_mark_center = RXNMarkCenter()
        self.enzr_norm_func = normalize_enzr_conf()

    def write_logs(self, logs: str) -> None:
        
        if not os.path.exists('output/' + self.log_folder):
            os.makedirs('output/' + self.log_folder)
            
        file_object = open('output/' + self.log_folder + '/' + self.log_time_stamp + '__logs.txt', 'a')
        file_object.write('[' + str(datetime.datetime.now()) + '] \t' + logs + "\n")
        file_object.close()
        
    def smi_tokenizer(self, smi: str) -> str:
        """
        Tokenize a SMILES molecule or reaction. Modified for the special tagging character "!".
        """
        pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\!|\$|\%[0-9]{2}|[0-9])"
        regex = re.compile(pattern)
        tokens = [token for token in regex.findall(smi)]
        assert smi == ''.join(tokens)
        return ' '.join(tokens)

    def canonicalize_smiles(self, smiles: str) -> str:
        '''
        Molecule canonicalization that does not change the SMILES order of molecules in case of multiple molecules.
        Also neutralizes any charge of the molecules.
        
        Args:
            smiles (str): SMILES string of the molecule(s).
        
        Returns:
            str: Canonicalized SMILES string of the molecule(s).
        '''
        returned = []
        any_error = False
        for molecule in smiles.split('.'):
            molecule = self.neutralize_smi(molecule)
            mol = Chem.MolFromSmiles(molecule)
            if mol is not None:
                returned.append(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True))
            else: 
                any_error = True
        if not any_error:
            return '.'.join(returned)
        else:
            return ''
    
    def neutralize_smi(self, smiles: str) -> str:        # from: https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
        if '-' in smiles or '+' in smiles:
            try:
                mol = Chem.MolFromSmiles(smiles)
                pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
                at_matches = mol.GetSubstructMatches(pattern)
                at_matches_list = [y[0] for y in at_matches]
                if len(at_matches_list) > 0:
                    for at_idx in at_matches_list:
                        atom = mol.GetAtomWithIdx(at_idx)
                        chg = atom.GetFormalCharge()
                        hcount = atom.GetTotalNumHs()
                        atom.SetFormalCharge(0)
                        atom.SetNumExplicitHs(hcount - chg)
                        atom.UpdatePropertyCache()
                return Chem.MolToSmiles(mol)
            except:
                return smiles
        else:
            return smiles

    def Execute_Prediction(
        self, 
        SMILES_list: list, 
        Model_path: str, 
        beam_size: int = 3, 
        batch_size: int = 64, 
        untokenize_output: bool = True
        ) -> Tuple[List, List]:
        '''
        This function executes the prediction of a given model for a given list of SMILES.
        
        Inputs:
            SMILES_list: list of SMILES to be predicted
            Model_path: path to the model to be used for prediction
            beam_size: beam size for the prediction output
            batch_size: batch size for the prediction
            untokenize_output: whether the output should be untokenized (True) or not (False)
        '''

        # Check Input:
        if not isinstance(SMILES_list, list):
            print('Input should be a list of SMILES')
            return []

        # Prepare Input File names for OpenNMT:
        timestamp = str(datetime.datetime.now()).replace(' ', '__').replace('-', '_').replace(':', '#')
        input_file = '{}input__{}.txt'.format(self.tmp_file_path, timestamp)
        output_file = '{}output__{}.txt'.format(self.tmp_file_path, timestamp)

        # Generate INPUT:
        textfile = open(input_file, "w")
        for element in SMILES_list: textfile.write(element + "\n")
        textfile.close()

        # Prepare arguments for OpenNMT:
        args = [
            "-beam_size",   str(beam_size), 
            "-n_best",      str(beam_size), 
            "-model",       str(Model_path), 
            "-src",         str(input_file), 
            "-output",      str(output_file), 
            "-batch_size",  str(batch_size), 
            "-max_length",  "1000", 
            "-log_probs",
            "-replace_unk"
        ]
        
        # Execute Predictions:
        parser = trsl._get_parser()
        opt = parser.parse_args(args)
        trsl.translate(opt)
        
        # Read Predictions:
        predictions = [[] for i in range(beam_size)]
        with open(output_file, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if untokenize_output:
                    predictions[i % beam_size].append(''.join(line.strip().split(' ')))
                else:
                    predictions[i % beam_size].append(line.strip())
        probs = [[] for i in range(beam_size)]
        with open(output_file + '_log_probs', 'r') as f:
            for i, line in enumerate(f.readlines()):
                probs[i % beam_size].append(10**float(''.join(line.strip().split(' '))))

        try:
            os.remove(input_file)
            os.remove(output_file)
            os.remove(output_file + '_log_probs')
        except: pass

        return predictions, probs

    def get_reagent_prediction_USPTO(self, df_prediction_Forw, beam_size_output=3):

        # Prepare Input for Reagent prediction model:
        Retro = [self.smi_tokenizer(el) for el in df_prediction_Forw['Retro']]
        Target = [self.smi_tokenizer(el) for el in df_prediction_Forw['Target']]
        Reag_Pred_from_Reaction = [Retro[el] + ' > > ' + Target[el] for el in range(0, len(Retro))]

        # Execute Reagent Prediction:
        predictions, probs = self.Execute_Prediction(
            SMILES_list = Reag_Pred_from_Reaction,
            Model_path = self.USPTO_T2_path,
            beam_size = beam_size_output
        )

        pred = []
        for i in range(0, beam_size_output):
            pred.append(predictions[i])

        return pred
    
    def get_reagent_prediction_ENZR(self, df_prediction_Forw, beam_size_output=3):
        
        # Prepare Input for Reagent prediction model:
        Retro = [self.smi_tokenizer(el) for el in df_prediction_Forw['Retro']]
        Target = [self.smi_tokenizer(el) for el in df_prediction_Forw['Target']]
        Reag_Pred_from_Reaction = ['ENZYME ' + Retro[el] + ' > > ' + Target[el] + ' ENZYME' for el in range(0, len(Retro))]
        
        # Execute Reagent Prediction:
        predictions, probs = self.Execute_Prediction(
            SMILES_list = Reag_Pred_from_Reaction,
            Model_path = self.ENZR_T2_path,
            beam_size = beam_size_output, 
            untokenize_output = False
        )
        
        pred = []
        for i in range(0, beam_size_output):
            pred.append(predictions[i])
            
        return pred
    
    def get_AutoTags(self, target: str, ini_smiles: str, AutoTagging_Beam_Size: int, Model: str = '') -> list:
        '''
            Returns a list of AutoTags for a given target SMILES
        '''
        
        list_retro_auto_tag, _ = self.Execute_Prediction([self.smi_tokenizer(target)], Model_path=Model, beam_size=AutoTagging_Beam_Size)
            
        #Check if those generated tags are valid and still representing the same molecule:
        list_retro_auto_tag_curated = []
        
        # For every generated tag:
        for mol_tagged in list_retro_auto_tag:
            
            # Remove the tags, the canonicalize and compare to the original SMILES:
            untagged = ''
            try:
                rxn_mapped = rdChemReactions.ReactionFromSmarts('CCC>>' + mol_tagged[0], useSmiles=True)
                rdChemReactions.SanitizeRxn(rxn_mapped)
                
                for product in rxn_mapped.GetProducts():
                    for atom in product.GetAtoms():
                        if atom.GetAtomMapNum() == 1:
                            atom.SetAtomMapNum(0)
                
                # Get Smiles
                rdChemReactions.SanitizeRxn(rxn_mapped)
                untagged = rdChemReactions.ReactionToSmiles(rxn_mapped, canonical=True).split('>>')[1]
                untagged = self.canonicalize_smiles(self.neutralize_smi(untagged))
                
            except:
                # case if not able to load the molecule
                untagged = ''
            
            if untagged == ini_smiles:
                list_retro_auto_tag_curated.append(mol_tagged[0])
        
        list_retro_auto_tag_curated = list(set(list_retro_auto_tag_curated))
        
        return [self.smi_tokenizer(self.rxn_mark_center.convert_smarts_to_smiles_alt_tagging(el)) for el in list_retro_auto_tag_curated]
    
    def Make_Forward_Prediction(
        self, 
        df_prediction_Forw, 
        Fwd_ENZ_Reag_Pred=False, 
        Fwd_USPTO_Reag_Pred=True, 
        USPTO_Reag_Beam_Size=3, 
        log=False
        ):

        # Need entire reformating to cover cases when we want individual models but also multiple
        forw_df = []

        # Set of predicted Enzyme Names:        [[NEW]]
        if Fwd_ENZ_Reag_Pred:
            # Selecting ENZR subset:
            df_prediction_Forw_sub = df_prediction_Forw[df_prediction_Forw['T1_Model'] == self.ENZR_T1_path].copy()
            
            if log: self.write_logs('Predicting reagents for ENZR fwd Pred, {} predictions...'.format(len(df_prediction_Forw_sub)))
            reagent_set = self.get_reagent_prediction_ENZR(df_prediction_Forw_sub, beam_size_output=USPTO_Reag_Beam_Size)
            
            for i in range(0, len(reagent_set)):
                df_prediction_Forw_sub['Forward_Model'] = self.ENZR_T3_path
                df_prediction_Forw_sub['Reagents'] = reagent_set[i]
                forw_df.append(df_prediction_Forw_sub.copy())
            if log: self.write_logs('ENZR Reagents set predicted')
            
        # Set of predicted USPTO Reagents:
        if Fwd_USPTO_Reag_Pred:
            # Selecting USPTO subset:
            df_prediction_Forw_sub = df_prediction_Forw[df_prediction_Forw['T1_Model'] == self.USPTO_T1_path].copy()
            
            if log: self.write_logs('Predicting reagents for USPTO fwd Pred, {} predictions...'.format(len(df_prediction_Forw_sub)))
            reagent_set = self.get_reagent_prediction_USPTO(df_prediction_Forw_sub, beam_size_output=USPTO_Reag_Beam_Size)
            
            for i in range(0, len(reagent_set)):
                df_prediction_Forw_sub['Forward_Model'] = self.USPTO_T3_path
                df_prediction_Forw_sub['Reagents'] = reagent_set[i]
                forw_df.append(df_prediction_Forw_sub.copy())
            if log: self.write_logs('USPTO Reagents set predicted')

        # Concatenate it all into DF:
        forw_df = pd.concat(forw_df).sort_values(by=['index', 'Reagents']).reset_index(drop=True)

        forw_df['Forward_Model_Input'] = ''
        forw_df['Forward_Prediction'] = ''
        forw_df['Prob_Forward_Prediction_1'] = 0.0
        forw_df['Prob_Forward_Prediction_2'] = 0.0          # 2nd was useful for scoring based on condifence score difference between 1st and 2nd prediction

        # Tokenization:
        for element in range(0, len(forw_df)):
            if forw_df.at[element, 'Forward_Model'] == self.ENZR_T3_path:
                forw_df.at[element, 'Forward_Model_Input'] = 'ENZYME ' + self.smi_tokenizer(forw_df.at[element, 'Retro']) + ' > ' + forw_df.at[element, 'Reagents'] + ' ENZYME'
            elif  forw_df.at[element, 'Forward_Model'] == self.USPTO_T3_path:
                forw_df.at[element, 'Forward_Model_Input'] = self.smi_tokenizer(forw_df.at[element, 'Retro']) + ' > ' + self.smi_tokenizer(forw_df.at[element, 'Reagents'])

        # T3 ENZR Forward Prediction:
        if Fwd_ENZ_Reag_Pred:
            if log: self.write_logs('ENZR_Reag_Pred Forward prediction...')
            predictions, probs = self.Execute_Prediction(list(forw_df[forw_df['Forward_Model'] == self.ENZR_T3_path]['Forward_Model_Input']), Model_path=self.ENZR_T3_path, beam_size=3)
            
            if len(forw_df.loc[forw_df['Forward_Model'] == self.ENZR_T3_path, 'Forward_Prediction']) == len(predictions[0]):
                forw_df.loc[forw_df['Forward_Model'] == self.ENZR_T3_path, 'Forward_Prediction'] = predictions[0]
                forw_df.loc[forw_df['Forward_Model'] == self.ENZR_T3_path, 'Prob_Forward_Prediction_1'] = self.enzr_norm_func(probs[0])
                forw_df.loc[forw_df['Forward_Model'] == self.ENZR_T3_path, 'Prob_Forward_Prediction_2'] = self.enzr_norm_func(probs[1])
            else:
                if log: self.write_logs('Lenghts predictions mismatch: {}, vs {}'.format(len(forw_df.loc[forw_df['Forward_Model'] == self.ENZR_T3_path, 'Forward_Prediction']), len(predictions[0])))
        
        # T3 USPTO Forward Prediction:
        if Fwd_USPTO_Reag_Pred:
            if log: self.write_logs('USPTO_T3 Forward prediction...')
            predictions, probs = self.Execute_Prediction(list(forw_df[forw_df['Forward_Model'] == self.USPTO_T3_path]['Forward_Model_Input']), Model_path=self.USPTO_T3_path, beam_size=3)
            
            if len(forw_df.loc[forw_df['Forward_Model'] == self.USPTO_T3_path, 'Forward_Prediction']) == len(predictions[0]):
                forw_df.loc[forw_df['Forward_Model'] == self.USPTO_T3_path, 'Forward_Prediction'] = predictions[0]
                forw_df.loc[forw_df['Forward_Model'] == self.USPTO_T3_path, 'Prob_Forward_Prediction_1'] = probs[0]
                forw_df.loc[forw_df['Forward_Model'] == self.USPTO_T3_path, 'Prob_Forward_Prediction_2'] = probs[1]
            else:
                if log: self.write_logs('Lenghts predictions mismatch: {}, vs {}'.format(len(forw_df.loc[forw_df['Forward_Model'] == self.USPTO_T3_path, 'Forward_Prediction']), len(predictions[0])))

        if log: self.write_logs('All Forward predictions done, start canonicalization...')

        # Canonicalize_smiles:
        for element in range(0, len(forw_df)):
            forw_df.at[element, 'Forward_Prediction'] = self.canonicalize_smiles(forw_df.at[element, 'Forward_Prediction'])

        if log: self.write_logs('Canonicalization done.')

        return forw_df

    def get_highest_similarity_to_target(self, target_smiles, predictions):
        target_mol = Chem.MolFromSmarts(target_smiles)
        target_fp = Chem.RDKFingerprint(target_mol)
        highest_similarity = 0.0

        for element in predictions.split('.'):
            retro_mol = Chem.MolFromSmarts(element)
            retro_fp = Chem.RDKFingerprint(retro_mol)

            similarity = DataStructs.DiceSimilarity(target_fp, retro_fp)

            if similarity > highest_similarity:
                highest_similarity = similarity

        return highest_similarity

    def _remove_predictions_similar_to_target(self, df: pd.DataFrame) -> pd.DataFrame:
        '''
        Remove predictions that are too similar to the target compound.
        
        Args:
            df (pd.DataFrame): DataFrame containing the predictions.
            
        Returns:
            pd.DataFrame: DataFrame containing the predictions.
        '''
        
        # Compute Similarities:
        df['Highest_Similarity'] = [self.get_highest_similarity_to_target(df.at[x, 'Target'], df.at[x, 'Retro']) for x in range(0, len(df))]
        # Filter out what is too similar:
        del df['index']
        df_filtered = df[df['Highest_Similarity'] <= 0.95].reset_index(drop=True).reset_index().copy()
        
        return df_filtered
    
    def _get_list_tags(self, Random_Tagging, AutoTagging, AutoTagModel, Substructure_Tagging, mark_count, neighbors, list_substructures, target, SMILES, AutoTagging_Beam_Size, mark_locations_filter, log):
        '''
        The function tags atoms in the target molecule.
        
        Args:
            Random_Tagging (bool): If True, the function will mark random atoms in the target molecule.
            AutoTagging (bool): If True, the function will mark the atoms in the target molecule using the AutoTagging model.
            AutoTagModel (str): The path to the AutoTagging model.
            Substructure_Tagging (bool): If True, the function will mark the atoms in the target molecule using the Substructure_Tagging model.
            mark_count (int): The number of atoms to mark.
            neighbors (int): The number of neighbors to consider when marking atoms.
            target (str): The target molecule.
            SMILES (str): The SMILES of the target molecule.
            AutoTagging_Beam_Size (int): The beam size to use when predicting with the AutoTagging model.
            mark_locations_filter (list): The list of locations to filter out.
            log (bool): If True, the function will write logs.
            
        Returns:
            list_retro (list): The list of tagged atoms.
        '''
        list_retro = []
        
        if Random_Tagging:
            for i in range(mark_count):
                list_retro += self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=i+1, neighbors=neighbors, tokenized=True)
        
        if AutoTagging:
            list_retro += self.get_AutoTags(target=target, ini_smiles=SMILES, AutoTagging_Beam_Size=AutoTagging_Beam_Size, Model=AutoTagModel)
        
        if Substructure_Tagging:
            list_retro += self.rxn_mark_center.Mark_matching_substructures(mol_SMILES=target, list_conditionnal_substructures_tags=list_substructures, tokenized=True)
        
        # Remove duplicate tags:
        list_retro = list(set(list_retro))

        # Random tag filtering:
        if mark_locations_filter < 1:   list_retro = sample(list_retro, int(round(len(list_retro)*mark_locations_filter, 0)))
        
        # Check if list_retro is empty:
        if len(list_retro) == 0 and log: self.write_logs('WARNING: length of list_retro = 0')
        
        return list_retro
    
    @staticmethod
    def _filter_best_forward_by_id(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Filters the best forward predictions by unique ID among the multiple T2 predictions.
        
        Args:
            df (pd.DataFrame): The dataframe containing the forward predictions.
        
        Returns:
            pd.DataFrame: The filtered dataframe.
        '''
        assert isinstance(df, pd.DataFrame), 'df must be a pandas DataFrame.'
        
        df['Best_Forw'] = False
        df['Unique_ID'] = df['T1_Model'].astype(str) + '.' + df['ID_Tag'].astype(str) + '.' + df['ID_beam'].astype(str)

        for unique_id in set(df['Unique_ID']):
            index_ = df[df['Unique_ID'] == unique_id].sort_values(by='Prob_Forward_Prediction_1', ascending=False).index[0]
            df.at[index_, 'Best_Forw'] = True

        del df['index']
        df = df[df['Best_Forw'] == True].reset_index(drop=True).reset_index()
        del df['Best_Forw']
        del df['Unique_ID']
        
        return df
    
    @staticmethod
    def _remove_prediction_bad_confidence_score_ratio(df: pd.DataFrame) -> pd.DataFrame:
        '''
        Remove predictions where the confidence score ratio does not meet the requirements.
        
        Args:
            df (pd.DataFrame): DataFrame containing the predictions.
        
        Returns:
            pd.DataFrame: DataFrame containing the predictions after filtering.
        '''
        
        # Compute ratio:
        df['confidence_filter'] = [True if (df.at[x, 'Prob_Forward_Prediction_1'] > 0.6 or df.at[x, 'Prob_Forward_Prediction_1'] > df.at[x, 'Prob_Forward_Prediction_2'] + 0.2) else False for x in range(0, len(df))]
        
        del df['index']
        df_filtered2 = df[df['confidence_filter'] == True].reset_index(drop=True).reset_index().copy()
        del df_filtered2['confidence_filter']
        
        return df_filtered2

    def Execute_Retro_Prediction(
        self, 
        SMILES, 
        mark_count = 3, 
        neighbors = True, 
        Random_Tagging = True, 
        AutoTagging = True, 
        AutoTagging_Beam_Size = 50, 
        Substructure_Tagging = True, 
        Retro_ENZR = True, 
        Retro_USPTO = True, 
        Fwd_ENZ_Reag_Pred = True, 
        Fwd_USPTO_Reag_Pred = True, 
        USPTO_Reag_Beam_Size = 3, 
        similarity_filter = False, 
        confidence_filter = False, 
        Retro_beam_size = 5, 
        mark_locations_filter = 1, 
        log = False
        ):

        target = self.canonicalize_smiles(SMILES)
        if log: self.write_logs('Retro prediction on mol ' + str(target))
        
        if Retro_USPTO:
            list_retro_USPTO = self._get_list_tags(
                Random_Tagging = Random_Tagging, 
                AutoTagging = AutoTagging, 
                AutoTagModel = self.USPTO_AutoTag_path, 
                Substructure_Tagging = Substructure_Tagging, 
                mark_count = mark_count, 
                neighbors = neighbors, 
                list_substructures = self.list_substructures, 
                target = target, 
                SMILES = SMILES, 
                AutoTagging_Beam_Size = AutoTagging_Beam_Size, 
                mark_locations_filter = mark_locations_filter, 
                log = log, 
            )
        else: list_retro_USPTO = []

        if Retro_ENZR:
            list_retro_ENZR = self._get_list_tags(
                Random_Tagging = False,     # Random tagging is not used for ENZR
                AutoTagging = AutoTagging, 
                AutoTagModel = self.ENZR_AutoTag_path, 
                Substructure_Tagging = Substructure_Tagging, 
                mark_count = mark_count, 
                neighbors = neighbors, 
                list_substructures = self.list_substructures_ENZR, 
                target = target, 
                SMILES = SMILES, 
                AutoTagging_Beam_Size = AutoTagging_Beam_Size, 
                mark_locations_filter = mark_locations_filter, 
                log = log, 
            )
        else: list_retro_ENZR = []
        
        concat_models = []
        
        if len(list_retro_USPTO) + len(list_retro_ENZR) == 0:
            curr_model = pd.DataFrame(columns=['T1_Model', 'ID', 'ID_Tag', 'ID_beam', 'Target', 'Tagged_Target', 'Retro', 'Retro_Conf'])
            concat_models.append(curr_model.copy())
        
        # Execute Retro prediction on USPTO model:
        if Retro_ENZR and len(list_retro_ENZR) > 0:
            if log: self.write_logs('Retro prediction on ENZR model ' + str(len(list_retro_ENZR)) + ' marking examples...')
            current_model = self.ENZR_T1_path
            predictions, probs = self.Execute_Prediction(SMILES_list=['ENZYME ' + el + ' ENZYME' for el in list_retro_ENZR], Model_path=current_model, beam_size=Retro_beam_size)
            
            # Make DataFrame out of the predictions:
            curr_model = pd.DataFrame(['' for element in range(0, len(predictions)*len(predictions[0]))])
            curr_model['T1_Model'] =     str(current_model)
            curr_model['ID'] =           ['R_{}.{}'.format(str(each_TAG+1), str(each_beam+1)) for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['ID_Tag'] =       [each_TAG+1 for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['ID_beam'] =      [each_beam+1 for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Target'] =       target
            curr_model['Tagged_Target'] =   [list_retro_ENZR[each_TAG].replace(' ', '') for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Retro'] =        [self.canonicalize_smiles(predictions[each_beam][each_TAG]) for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Retro_Conf'] =   [probs[each_beam][each_TAG] for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            del curr_model[0]
            concat_models.append(curr_model.copy())
        
        if Retro_USPTO and len(list_retro_USPTO) > 0:
            if log: self.write_logs('Retro prediction on USPTO model ' + str(len(list_retro_USPTO)) + ' marking examples...')
            current_model = self.USPTO_T1_path
            predictions, probs = self.Execute_Prediction(SMILES_list=list_retro_USPTO, Model_path=current_model, beam_size=Retro_beam_size)

            # Make DataFrame out of the predictions:
            curr_model = pd.DataFrame(['' for element in range(0, len(predictions)*len(predictions[0]))])
            curr_model['T1_Model'] =     str(current_model)
            curr_model['ID'] =           ['R_{}.{}'.format(str(each_TAG+1), str(each_beam+1)) for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['ID_Tag'] =       [each_TAG+1 for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['ID_beam'] =      [each_beam+1 for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Target'] =       target
            curr_model['Tagged_Target'] =   [list_retro_USPTO[each_TAG].replace(' ', '') for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Retro'] =        [self.canonicalize_smiles(predictions[each_beam][each_TAG]) for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Retro_Conf'] =   [probs[each_beam][each_TAG] for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            del curr_model[0]
            concat_models.append(curr_model.copy())

        # Concatenate all Retro Predictions Across All Models, remove duplicates, remove too long sequences:
        df_prediction = pd.concat(concat_models)
        df_prediction = df_prediction.drop_duplicates(subset=['T1_Model', 'Retro'], keep='first').copy()
        df_prediction = df_prediction[df_prediction['Retro'] != ''].copy()
        df_prediction = df_prediction[df_prediction['Retro'].str.len() < 999]
        if log: self.write_logs('Resulted in {} unique T1 retro predictions to be forward validated.'.format(len(df_prediction)))

        df_prediction_Forw = df_prediction.copy().reset_index(drop=True).reset_index()     # Do not drop indexes
        df_prediction_Forw_2 = self.Make_Forward_Prediction(
            df_prediction_Forw = df_prediction_Forw, 
            Fwd_ENZ_Reag_Pred = Fwd_ENZ_Reag_Pred, 
            Fwd_USPTO_Reag_Pred = Fwd_USPTO_Reag_Pred, 
            USPTO_Reag_Beam_Size = USPTO_Reag_Beam_Size, 
            log = log
            )
        del df_prediction_Forw_2['Forward_Model_Input']

        # Save unvalidated routes for debugging or custom predictions: #DEBUG
        df_prediction_ALL_Debug_Backup = df_prediction_Forw_2.copy()
        
        # Keep predictions where T3 predicts the correct target, and target is not in the retro prediction:
        df_prediction_Forw_val = df_prediction_Forw_2[df_prediction_Forw_2['Target'] == df_prediction_Forw_2['Forward_Prediction']]
        df_prediction_Forw_val = df_prediction_Forw_val[df_prediction_Forw_val['Target'] != df_prediction_Forw_val['Retro']]
        
        # Remove low confident ENZR predictions:
        if Retro_ENZR:
            items = df_prediction_Forw_val[(df_prediction_Forw_val['Forward_Model'] == self.ENZR_T3_path) & (df_prediction_Forw_val['Prob_Forward_Prediction_1'] < self.ENZR_confidence_threshold)].index
            df_prediction_Forw_val.drop(items, inplace=True)
        
        del df_prediction_Forw_val['index']
        df_prediction_Forw_val = df_prediction_Forw_val.reset_index(drop=True).reset_index()

        # Keep The Best Forward Validated Prediction for each Tag of each model
        df_prediction_Forw_val = self._filter_best_forward_by_id(df_prediction_Forw_val)
        
        if log: self.write_logs('After forward validation, ' + str(len(df_prediction_Forw_val)) + ' unique predictions remain.')

        # Filter out precursors that are too similar to the target compound:
        if similarity_filter:   df_filtered = self._remove_predictions_similar_to_target(df_prediction_Forw_val)
        else:                   df_filtered = df_prediction_Forw_val.copy()
        if confidence_filter:   df_filtered2 = self._remove_prediction_bad_confidence_score_ratio(df_filtered)
        else:                   df_filtered2 = df_filtered.copy()

        return df_filtered2, df_prediction_ALL_Debug_Backup












