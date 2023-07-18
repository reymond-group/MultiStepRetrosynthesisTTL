from __future__ import division, unicode_literals

from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdChemReactions
from rdkit import RDLogger 
RDLogger.DisableLog('rdApp.*') 

import pandas as pd
import os
from random import sample
import re

import datetime
import subprocess

from rxnmarkcenter import RXNMarkCenter


class AugmentedSingleStepRetro:
    '''
        Class for Single Step Retrosynthesis.
        The input is passed over the Retrosynthesis Transformer model, followed by reagent prediction, then predictions are validated over the Forward Model.
    '''

    def __init__(self, log_folder = '', log_time_stamp = '', list_substructures = [['', '']]):
        
        self.Translate_py =  'makesingleretropredictions/onmt/translate.py'
        self.checkpoints =   'makesingleretropredictions/onmt/models/'

        ##      RETROSYNTHESIS:
        self.Retro_USPTO_No_Reag_Model = 'path_to_USPTO_retrosynthesis_model.pt'

        ##      REAGENT PREDICTION MODEL (from a given full disconnection / reaction SMILES):
        self.Reag_Prediction_USPTO_Model = 'path_to_reagent_prediction_model.pt'

        #       FORWARD MODELS:
        self.Frw_USPTO_with_Reag_Model = 'path_to_forward_USPTO_model_trained_including_reagents.pt'
        self.Frw_Std_USPTO_no_Reag_Model = 'path_to_forward_USPTO_model_trained_without_reagents.pt'
        
        #       AutoTag Model: (Thakkar)
        self.AutoTag_USPTO_Model = 'path_to_autotag_model.pt'

        os.environ['MKL_THREADING_LAYER'] = 'GNU'
        self.log_folder = log_folder
        self.log_time_stamp = log_time_stamp

        self.list_substructures = list_substructures
        self.rxn_mark_center = RXNMarkCenter()

    def write_logs(self, logs):
        
        if not os.path.exists('output/' + self.log_folder):
            os.makedirs('output/' + self.log_folder)
            
        file_object = open('output/' + self.log_folder + '/' + self.log_time_stamp + '__logs.txt', 'a')
        file_object.write('[' + str(datetime.datetime.now()) + '] \t' + logs + "\n")
        file_object.close()
        
    def smi_tokenizer(self, smi):
        """
        Tokenize a SMILES molecule or reaction. Modified for the special tagging character "!".
        """
        pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\!|\$|\%[0-9]{2}|[0-9])"
        regex = re.compile(pattern)
        tokens = [token for token in regex.findall(smi)]
        assert smi == ''.join(tokens)
        return ' '.join(tokens)

    def canonicalize_smiles(self, smiles):
        returned = []
        any_error = False
        for molecule in smiles.split('.'):
            molecule = self.neutralize_atoms(molecule)
            mol = Chem.MolFromSmiles(molecule)
            if mol is not None:
                returned.append(Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True))
            else: 
                any_error = True
        if not any_error:
            return '.'.join(returned)
        else:
            return ''
    
    def neutralize_atoms(self, smiles):        # from: https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
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

    def Execute_Prediction(self, SMILES_list, MODEL, beam_size = 3):
        '''
        This function executes the prediction of a given model for a given list of SMILES.
        
        '''

        # Check Input:
        if not isinstance(SMILES_list, list):
            print('Input should be a list of SMILES')
            return []

        # Prepare Input File names for OpenNMT:
        timestamp = str(datetime.datetime.now()).replace(' ', '__').replace('-', '_').replace(':', '#')
        input_file = 'makesingleretropredictions/onmt/predictions/input__' + timestamp + '.txt'
        output_file = 'makesingleretropredictions/onmt/predictions/output__' + timestamp + '.txt'

        # Generate INPUT:
        textfile = open(input_file, "w")
        for element in SMILES_list: textfile.write(element + "\n")
        textfile.close()

        batch_size = '64'
        bashCommand = ''

        # Execute Predictions:
        if MODEL == 'Retro_USPTO':
            bashCommand = 'python ' + self.Translate_py + ' -beam_size ' + str(beam_size) + ' -n_best ' + str(beam_size) + ' -model ' + self.checkpoints + self.Retro_USPTO_No_Reag_Model + ' -src ' + input_file + ' -output ' + output_file + ' -batch_size ' + batch_size + ' -replace_unk -max_length 1000 -log_probs'
        elif MODEL == 'Forw_USPTO_No_Reag ':         # FORWARD MOL TRANSFORMER ==> NO REAG
            bashCommand = 'python ' + self.Translate_py + ' -beam_size ' + str(beam_size) + ' -n_best ' + str(beam_size) + ' -model ' + self.checkpoints + self.Frw_Std_USPTO_no_Reag_Model + ' -src ' + input_file + ' -output ' + output_file + ' -batch_size ' + batch_size + ' -replace_unk -max_length 1000 -log_probs'
        elif MODEL == 'Forw_USPTO_Reag':            # FORWARD MOL TRANSFORMER ==> WITH REAG
            bashCommand = 'python ' + self.Translate_py + ' -beam_size ' + str(beam_size) + ' -n_best ' + str(beam_size) + ' -model ' + self.checkpoints + self.Frw_USPTO_with_Reag_Model + ' -src ' + input_file + ' -output ' + output_file + ' -batch_size ' + batch_size + ' -replace_unk -max_length 1000 -log_probs'        
        elif MODEL == 'Reagent_Pred_USPTO':         # PREDICT REAGENTS
            bashCommand = 'python ' + self.Translate_py + ' -beam_size ' + str(beam_size) + ' -n_best ' + str(beam_size) + ' -model ' + self.checkpoints + self.Reag_Prediction_USPTO_Model + ' -src ' + input_file + ' -output ' + output_file + ' -batch_size ' + batch_size + ' -replace_unk -max_length 1000 -log_probs'
        elif MODEL == 'AutoTag':                    # TAGGING TARGETS
            bashCommand = 'python ' + self.Translate_py + ' -beam_size ' + str(beam_size) + ' -n_best ' + str(beam_size) + ' -model ' + self.checkpoints + self.AutoTag_USPTO_Model + ' -src ' + input_file + ' -output ' + output_file + ' -batch_size ' + batch_size + ' -replace_unk -max_length 1000 -log_probs'

        # Execute prediction:
        subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE).communicate()
        
        # Read Predictions:
        predictions = [[] for i in range(beam_size)]
        with open(output_file, 'r') as f:
            for i, line in enumerate(f.readlines()):
                predictions[i % beam_size].append(''.join(line.strip().split(' ')))
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

    def get_reagent_prediction_USPTO(self, df_prediction_Forw, beam_size_output=1):

        # Prepare Input for Reagent prediction model:
        Retro = [self.smi_tokenizer(el) for el in df_prediction_Forw['Retro']]
        Target = [self.smi_tokenizer(el) for el in df_prediction_Forw['Target']]
        Reag_Pred_from_Reaction = [Retro[el] + ' > > ' + Target[el] for el in range(0, len(Retro))]

        # Execute Reagent Prediction:
        predictions, probs = self.Execute_Prediction(
            SMILES_list = Reag_Pred_from_Reaction,
            MODEL = 'Reagent_Pred_USPTO',
            beam_size = beam_size_output
        )

        pred = []
        for i in range(0, beam_size_output):
            pred.append(predictions[i])

        pred.append('')
        return pred
    
    def get_AutoTags(self, target, ini_smiles, AutoTagging_Beam_Size):
        '''
            Returns a list of AutoTags for a given target SMILES
        '''
        
        list_retro_auto_tag, _ = self.Execute_Prediction([self.smi_tokenizer(target)], MODEL='AutoTag', beam_size=AutoTagging_Beam_Size)
            
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
                untagged = self.canonicalize_smiles(self.neutralize_atoms(untagged))
                
            except:
                # case if not able to load the molecule
                untagged = ''
            
            if untagged == ini_smiles:
                list_retro_auto_tag_curated.append(mol_tagged[0])
        
        list_retro_auto_tag_curated = list(set(list_retro_auto_tag_curated))
        
        return [self.smi_tokenizer(self.rxn_mark_center.convert_smarts_to_smiles_alt_tagging(el)) for el in list_retro_auto_tag_curated]
    
    def Make_Forward_Prediction(self, df_prediction_Forw, Std_Fwd_USPTO=True, Fwd_USPTO_Reag_Pred=False, USPTO_Reag_Beam_Size=1, log=False):

        # Need entire reformating to cover cases when we want individual models but also multiple
        forw_df = []

        # Default USPTO without reagent:
        if Std_Fwd_USPTO:
            df_prediction_Forw['Forward_Model'] = 'Forw_USPTO_No_Reag'
            df_prediction_Forw['Reagents'] = ''
            forw_df.append(df_prediction_Forw.copy())
        
        # Set of predicted USPTO Reagents:
        if Fwd_USPTO_Reag_Pred:
            if log: self.write_logs('Predicting reagents for USPTO fwd Pred, ' + str(len(df_prediction_Forw)) +  ' predictions...')
            reagent_set = self.get_reagent_prediction_USPTO(df_prediction_Forw, beam_size_output=USPTO_Reag_Beam_Size)

            for i in range(0, len(reagent_set)):
                df_prediction_Forw['Forward_Model'] = 'Forw_USPTO_Reag'
                df_prediction_Forw['Reagents'] = reagent_set[i]
                forw_df.append(df_prediction_Forw.copy())
            if log: self.write_logs('Reagents set predicted')

        # Concatenate it all into DF:
        forw_df = pd.concat(forw_df).sort_values(by=['index', 'Reagents']).reset_index(drop=True)

        forw_df['Forward_Model_Input'] = ''
        forw_df['Forward_Prediction'] = ''
        forw_df['Prob_Forward_Prediction_1'] = 0.0
        forw_df['Prob_Forward_Prediction_2'] = 0.0

        # Tokenization:
        for element in range(0, len(forw_df)):
            if forw_df.at[element, 'Forward_Model'] == 'Forw_USPTO_No_Reag':
                forw_df.at[element, 'Forward_Model_Input'] = self.smi_tokenizer(forw_df['Retro'][element])
            elif  forw_df['Forward_Model'][element] == 'Forw_USPTO_Reag':
                forw_df.at[element, 'Forward_Model_Input'] = self.smi_tokenizer(forw_df['Retro'][element]) + ' > ' + self.smi_tokenizer(forw_df['Reagents'][element])

        # Std NO REAGENT Molecular Transformer:
        if Std_Fwd_USPTO:
            if log: self.write_logs('Forw_USPTO_No_Reag Forward prediction...')
            predictions, probs = self.Execute_Prediction(list(forw_df[forw_df['Forward_Model'] == 'Forw_USPTO_No_Reag']['Forward_Model_Input']), MODEL='Forw_USPTO_No_Reag', beam_size=3)
            
            if len(forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_No_Reag', 'Forward_Prediction']) == len(predictions[0]):
                forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_No_Reag', 'Forward_Prediction'] = predictions[0]
                forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_No_Reag', 'Prob_Forward_Prediction_1'] = probs[0]
                forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_No_Reag', 'Prob_Forward_Prediction_2'] = probs[1]
            else:
                if log: self.write_logs('Lenghts predictions mismatch: ' + str(len(forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_No_Reag', 'Forward_Prediction'])) + ', vs  ' + str(len(predictions[0])))
        
        # Molecular Transformer INCLUDING reagents (predicted):
        if Fwd_USPTO_Reag_Pred:
            if log: self.write_logs('Forw_USPTO_Reag Forward prediction...')
            predictions, probs = self.Execute_Prediction(list(forw_df[forw_df['Forward_Model'] == 'Forw_USPTO_Reag']['Forward_Model_Input']), MODEL='Forw_USPTO_Reag', beam_size=3)
            
            if len(forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_Reag', 'Forward_Prediction']) == len(predictions[0]):
                forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_Reag', 'Forward_Prediction'] = predictions[0]
                forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_Reag', 'Prob_Forward_Prediction_1'] = probs[0]
                forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_Reag', 'Prob_Forward_Prediction_2'] = probs[1]
            else:
                if log: self.write_logs('Lenghts predictions mismatch: ' + str(len(forw_df.loc[forw_df['Forward_Model'] == 'Forw_USPTO_Reag', 'Forward_Prediction'])) + ', vs  ' + str(len(predictions[0])))

        if log: self.write_logs('All Forward prediction done, start canonicalization...')

        # Canonicalize_smiles:
        for element in range(0, len(forw_df)):
            forw_df.at[element, 'Forward_Prediction'] = self.canonicalize_smiles(forw_df.at[element, 'Forward_Prediction'])

        if log: self.write_logs('Canonicalization done.')

        return forw_df

    def Execute_Retro_Prediction(self, SMILES, mark_count=2, neighbors=True, Random_Tagging=True, AutoTagging=False, AutoTagging_Beam_Size=100, Substructure_Tagging=True, Retro_USPTO=True, Std_Fwd_USPTO=True, Fwd_USPTO_Reag_Pred=False, USPTO_Reag_Beam_Size=1, confidence_filter=True, Retro_beam_size=15, mark_locations_filter=1, log=False):

        target = self.canonicalize_smiles(SMILES)
        if log: self.write_logs('Retro prediction on mol ' + str(target))
        list_retro = []
        
        if Random_Tagging:
            if mark_count == 1:
                list_retro +=   self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=mark_count, neighbors=neighbors, tokenized=True)
            elif mark_count == 2:
                list_retro +=   self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=mark_count-1, neighbors=neighbors, tokenized=True)
                list_retro +=   self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=mark_count, neighbors=neighbors, tokenized=True)
            elif mark_count == 3:
                list_retro +=   self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=mark_count-2, neighbors=neighbors, tokenized=True)
                list_retro +=   self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=mark_count-1, neighbors=neighbors, tokenized=True)
                list_retro +=   self.rxn_mark_center.Mark_Random_Atoms(target, mark_count=mark_count, neighbors=neighbors, tokenized=True)
                
        if AutoTagging:
            list_retro += self.get_AutoTags(target=target, ini_smiles=SMILES, AutoTagging_Beam_Size=AutoTagging_Beam_Size)
        
        if Substructure_Tagging:
            list_retro += self.rxn_mark_center.Mark_matching_substructures(mol_SMILES=target, list_conditionnal_substructures_tags=self.list_substructures, tokenized=True)
        
        # Remove duplicate tags:
        list_retro = list(set(list_retro))

        # Random tag filtering:
        if mark_locations_filter < 1:
            list_retro = sample(list_retro, int(round(len(list_retro)*mark_locations_filter, 0)))
        if len(list_retro) == 0:
            if log: self.write_logs('WARNING:  length of list_retro = 0')

        concat_models = []
        
        if len(list_retro) == 0:
            curr_model = pd.DataFrame(columns=['ID', 'ID_Tag', 'ID_beam', 'Target', 'Tag_Target', 'Retro', 'Retro_Conf'])
            concat_models.append(curr_model.copy())

        if Retro_USPTO and len(list_retro) > 0:
            if log: self.write_logs('Retro prediction on USPTO model ' + str(len(list_retro)) + ' marking examples...')
            model = 'Retro_USPTO'
            predictions, probs = self.Execute_Prediction(SMILES_list=list_retro, MODEL='Retro_USPTO', beam_size=Retro_beam_size)

            # Make DataFrame out of the predictions:
            curr_model = pd.DataFrame(['' for element in range(0, len(predictions)*len(predictions[0]))])
            curr_model['ID'] =           ['R_' + model + '_' + str(each_TAG+1) + '.' + str(each_beam+1) for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['ID_Tag'] =       [each_TAG+1 for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['ID_beam'] =      [each_beam+1 for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Target'] =       target
            curr_model['Tag_Target'] =   [list_retro[each_TAG].replace(' ', '') for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Retro'] =        [self.canonicalize_smiles(predictions[each_beam][each_TAG]) for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            curr_model['Retro_Conf'] =   [probs[each_beam][each_TAG] for each_TAG in range(0, len(predictions[0])) for each_beam in range(0, len(predictions))]
            del curr_model[0]
            concat_models.append(curr_model.copy())

        # Concatenate all Retro Predictions Across All Models, remove duplicates, remove too long sequences:
        df_prediction = pd.concat(concat_models)
        df_prediction = df_prediction.drop_duplicates(subset=['Retro'], keep='first').copy()
        df_prediction = df_prediction[df_prediction['Retro'] != ''].copy()
        df_prediction = df_prediction[df_prediction['Retro'].str.len() < 999]

        if log: self.write_logs('Resulted in ' + str(len(df_prediction)) + ' unique retro predictions to be forward validated.')

        df_prediction_Forw = df_prediction.copy().reset_index(drop=True).reset_index()     # Do not drop indexes

        # Save unvalidated routes for debugging:
        df_prediction_ALL_Debug_Backup = df_prediction_Forw.copy()

        df_prediction_Forw_2 = self.Make_Forward_Prediction(
            df_prediction_Forw=df_prediction_Forw, 
            Std_Fwd_USPTO=Std_Fwd_USPTO,
            Fwd_USPTO_Reag_Pred=Fwd_USPTO_Reag_Pred, 
            USPTO_Reag_Beam_Size=USPTO_Reag_Beam_Size,
            log=log)
        del df_prediction_Forw_2['Forward_Model_Input']

        # Keep The Best Forward Validated Prediction for each Tag
        df_prediction_Forw_val = df_prediction_Forw_2[df_prediction_Forw_2['Target'] == df_prediction_Forw_2['Forward_Prediction']]
        df_prediction_Forw_val = df_prediction_Forw_val[df_prediction_Forw_val['Target'] != df_prediction_Forw_val['Retro']]
        
        del df_prediction_Forw_val['index']
        df_prediction_Forw_val = df_prediction_Forw_val.reset_index(drop=True).reset_index()

        validated = list(set(list(df_prediction_Forw_val['ID'])))
        df_prediction_Forw_val['Best_Forw'] = False

        for best in validated:
            index__ = list(df_prediction_Forw_val[df_prediction_Forw_val['ID'] == best].sort_values(by='Prob_Forward_Prediction_1', ascending=False)['index'])[0]
            df_prediction_Forw_val.at[index__, 'Best_Forw'] = True

        del df_prediction_Forw_val['index']
        df_prediction_Forw_val = df_prediction_Forw_val[df_prediction_Forw_val['Best_Forw'] == True].reset_index(drop=True).reset_index()
        del df_prediction_Forw_val['Best_Forw']
        
        if log: self.write_logs('After forward validation, ' + str(len(df_prediction_Forw_val)) + ' unique predictions remain.')

        df_filtered = df_prediction_Forw_val.copy()

        if confidence_filter:
            # Compute ratio:
            df_filtered['confidence_filter'] = [True if (df_filtered['Prob_Forward_Prediction_1'][x] > 0.6 or df_filtered['Prob_Forward_Prediction_1'][x] > df_filtered['Prob_Forward_Prediction_2'][x] + 0.2) else False for x in range(0, len(df_filtered))]
            
            del df_filtered['index']
            df_filtered2 = df_filtered[df_filtered['confidence_filter'] == True].reset_index(drop=True).reset_index().copy()
            del df_filtered2['confidence_filter']
        else:
            df_filtered2 = df_filtered.copy()


        return df_filtered2, df_prediction_ALL_Debug_Backup





















