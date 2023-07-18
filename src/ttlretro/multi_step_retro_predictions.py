from __future__ import division, unicode_literals

from tqdm import tqdm
import pandas as pd
import numpy as np
import os
import copy
import datetime
import csv

from ttlretro.single_step_retro import SingleStepRetrosynthesis

from scscore import SCScorer
scs_scorer = SCScorer()

class MultiStepGraphRetro:
    '''
        Class for Multiple Step Retrosynthesis by Heuristic Best-First Tree Search.
    
    '''

    def __init__(
        self, 
        mark_count = 3, 
        neighbors = True, 
        Random_Tagging = True, 
        AutoTagging = False,
        AutoTagging_Beam_Size = 50,
        Substructure_Tagging = True, 
        Retro_USPTO = True, 
        Fwd_USPTO_Reag_Pred = True, 
        USPTO_Reag_Beam_Size = 3, 
        similarity_filter = False, 
        confidence_filter = False, 
        Retro_beam_size = 5, 
        max_mol_per_iteration = 20, 
        mark_locations_filter = 1.0, 
        step_penalties_rate = 0.8, 
        tree_only_best = 0.50, 
        tree_min_best = 15000, 
        tree_max_best = 50000, 
        tree_max_steps = 30, 
        branching_max_expansion = 40, 
        project_name = '', 
        Commercial_exclusions = [], 
        list_substructures_path = '', 
        log = False, 
        USPTO_AutoTag_path = '',
        USPTO_T1_path = '', 
        USPTO_T2_path = '',
        USPTO_T3_path = '',
        tmp_file_path = 'tmp/', 
        commercial_file_path = ''
        ):
        
        with open(commercial_file_path) as f:
            self.Commercial = f.read().splitlines()
        self.Commercial_exclusions = Commercial_exclusions

        self.mark_count = mark_count
        self.neighbors = neighbors
        self.Random_Tagging = Random_Tagging
        self.AutoTagging = AutoTagging
        self.AutoTagging_Beam_Size = AutoTagging_Beam_Size
        self.Substructure_Tagging = Substructure_Tagging
        self.Retro_USPTO = Retro_USPTO
        self.Fwd_USPTO_Reag_Pred = Fwd_USPTO_Reag_Pred
        self.USPTO_Reag_Beam_Size = USPTO_Reag_Beam_Size
        self.similarity_filter = similarity_filter
        self.confidence_filter = confidence_filter
        self.Retro_beam_size = Retro_beam_size
        self.max_mol_per_iteration = max_mol_per_iteration
        self.mark_locations_filter = mark_locations_filter
        self.step_penalties_rate = step_penalties_rate
        self.tree_only_best = tree_only_best
        self.tree_min_best = tree_min_best
        self.tree_max_best = tree_max_best
        self.tree_max_steps = tree_max_steps
        self.branching_max_expansion = branching_max_expansion

        self.project_name = project_name
        self.log = log
        
        
        self.USPTO_AutoTag_path = USPTO_AutoTag_path
        self.USPTO_T1_path = USPTO_T1_path
        self.USPTO_T2_path = USPTO_T2_path
        self.USPTO_T3_path = USPTO_T3_path
        
        self.tmp_file_path = tmp_file_path

        self.log_time_stamp = str(datetime.datetime.now()).split('.')[0].replace(' ', '__').replace('-', '_').replace(':', '') 

        if self.project_name == '':
            self.project_name = 'No_Name__' + self.log_time_stamp

        if not os.path.exists('output/' + self.project_name):
            os.makedirs('output/' + self.project_name)

        if list_substructures_path != '':
            with open(list_substructures_path, 'r') as f:
                reader = csv.reader(f)
                self.list_substructures = list(reader)
        else:
            self.list_substructures = [['', '']] 

        self.make_single_retropredictions = SingleStepRetrosynthesis(
            log_folder = self.project_name,
            log_time_stamp = self.log_time_stamp, 
            list_substructures = self.list_substructures, 
            USPTO_AutoTag_path = self.USPTO_AutoTag_path, 
            USPTO_T1_path = self.USPTO_T1_path, 
            USPTO_T2_path = self.USPTO_T2_path, 
            USPTO_T3_path = self.USPTO_T3_path, 
            tmp_file_path = self.tmp_file_path
        )

    def _if_commercial(self, smiles):
        
        if smiles in self.Commercial_exclusions: return False
        if smiles in self.Commercial:   return True
        else:                           
            # Check for charges, neutralize it all:
            if '-' in smiles or '+' in smiles:
                smiles_2 = self.make_single_retropredictions.neutralize_smi(smiles=smiles)
                if smiles_2 in self.Commercial_exclusions:  return False
                if smiles_2 in self.Commercial:             return True
                else:                                       return False
            else:
                return False

    def solved_status_update_reaction_predictions(self, predictions, step):

        predictions['Solved'] = [[self._if_commercial(mol) for mol in retro] for retro in [mol.split('.') for mol in predictions['Retro']]]
        predictions['Solved_overall'] = [True if all(el) else False for el in predictions['Solved']]
        predictions['Next_Iterated'] = [['Comm' if self._if_commercial(mol) else False for mol in retro] for retro in [mol.split('.') for mol in predictions['Retro']]]
        predictions['Step'] = step
        #predictions['Score_old'] = [predictions.at[el, 'Prob_Forward_Prediction_1'] * (np.prod([1-((scs_scorer.get_score_from_smi(smi)[1]-1)/4) for smi in predictions.at[el, 'Retro'].split('.')])) / (1-((scs_scorer.get_score_from_smi(predictions.at[el, 'Target'])[1]-1)/4)) for el in range(0, len(predictions))]
        predictions['Score'] = [predictions.at[el, 'Prob_Forward_Prediction_1'] * (np.prod([1-((scs_scorer.get_score_from_smi(smi)[1]-1)/4) if not self._if_commercial(smi) else 1 for smi in predictions.at[el, 'Retro'].split('.')])) for el in range(0, len(predictions))]        

        return predictions

    def connect_linked_reaction_predictions(self, predictions, next_pred_mol):
        
        #   UPDATE CONNECTIONS WITH PREVIOUS REACTIONS: (avoids repeating the same predictions)
        # For each predictions rxn:
        for each_retro in range(0, len(predictions)):
            # For each molecule in the Retro column...
            for each_mol in range(0, len(predictions.at[each_retro, 'Next_Iterated'])):
                #...that were potentially further iterated...
                if predictions.at[each_retro, 'Next_Iterated'][each_mol] == False:
                    #...that were for real further iterated...
                    if predictions.at[each_retro, 'Retro'].split('.')[each_mol] in next_pred_mol:
                        # Get the molecule:
                        mol = predictions.at[each_retro, 'Retro'].split('.')[each_mol]
                        
                        # Get the corresponding predictions in which it is involved:
                        sub_current_preds = predictions[predictions['Target'] == mol].copy()
                        if len(list(sub_current_preds['index'])) != 0:
                            # if entires are found, update the Next_Iterated:
                            predictions.at[each_retro, 'Next_Iterated'][each_mol] = list(sub_current_preds['index'])
                        else:
                            # if no entries are found, the retrosynthesis must have failed:
                            predictions.at[each_retro, 'Next_Iterated'][each_mol] = 'Failed'
            
            # Then update the potential new retro molecules just generated in advance before it get further iterated to avoid repeating the same predictions:
            if predictions.at[each_retro, 'Target'] in next_pred_mol:
                # then, for each of the Retro molecules...
                for each_mol in range(0, len(predictions.at[each_retro, 'Next_Iterated'])):
                    #...that will be potentially further iterated next...
                    if predictions.at[each_retro, 'Next_Iterated'][each_mol] == False:
                        
                        
                        # Get the molecule:
                        mol = predictions.at[each_retro, 'Retro'].split('.')[each_mol]
                        
                        # Get the corresponding predictions in which it is involved:
                        sub_current_preds = predictions[predictions['Target'] == mol].copy()
                        if len(list(sub_current_preds['index'])) != 0:
                            # if entires are found, update the Next_Iterated:
                            predictions.at[each_retro, 'Next_Iterated'][each_mol] = list(sub_current_preds['index'])
                        else:
                            # if no entries are found, it does not necessarily mean that the retrosynthesis failed as the molecule was not iterated yet:
                            i=0
            
        return predictions
    
    def multiply(self, List_of_List):
        LoL = List_of_List.copy()
        while len(LoL) > 1:
            take_list_1, take_list_2, New = LoL[0].copy(), LoL[1].copy(), []
            for el1 in take_list_1:
                is_el1_list = isinstance(el1, list)
                for el2 in take_list_2:
                    if is_el1_list:   
                        if isinstance(el2, list):   New.append(el1 + el2)
                        else:                       New.append(el1 + [el2])
                    else:
                        if isinstance(el2, list):   New.append([el1] + el2) 
                        else:                       New.append([el1, el2])
            LoL[0] = New
            del LoL[1]
        return LoL[0]

    def is_tree_expansion_possible(self, list_Next_Iterated):
        no_false = True
        one_list = False
        for element in list_Next_Iterated:
            if isinstance(element, list):
                one_list = True
            if element == False:
                no_false = False
        if one_list == True and no_false == True: return True
        else: return False

    def something_in_Still_Expandable(self, a_list):
        for el in a_list:
            if len(el) > 0:
                return True
        return False

    def get_solved_status_from_set_of_rxn(self, list_of_rxn, predictions):
        '''
        Returns if the reactions path have all been explored and if resolved. It looks for "False" or failed in the Next_Iterated column.
        #todo: function uses half the tree calculation time, should be optimized.
        '''

        any_false = False
        any_failed = False
        any_list_not_connected = False
        list_of_rxn_forward = copy.deepcopy(list_of_rxn)

        for rxn in list_of_rxn:
            list_of_rxn_forward.remove(rxn)
            Next_Iterated = predictions.at[rxn, 'Next_Iterated']
            for mol in Next_Iterated:
                if isinstance(mol, list):
                    found = False
                    for sub_rxn in mol:
                        if sub_rxn in list_of_rxn_forward: found = True
                    if found == False: any_list_not_connected = True

            if 'Failed' in Next_Iterated:    any_failed = True
            if 'False' in Next_Iterated:     any_false = True
            if False in Next_Iterated:       any_false = True

        if any_failed:  
            Solved_status = 'Failed'
        else:
            if any_false or any_list_not_connected:         Solved_status = 'No'
            else:                                           Solved_status = 'Yes'
            
        return Solved_status

    def Get_Tree_From_Reaction_DF(self, predictions, target_cpd):

        ini_reactions = list(predictions[predictions['Target'] == target_cpd]['index'])
        ini_base = []
        if self.log: self.make_single_retropredictions.write_logs('Start calculating tree...')

        for el in ini_reactions:
            ini_base.append([
                [], [el], 
                self.get_solved_status_from_set_of_rxn(list_of_rxn=[el], predictions=predictions), 
                float(predictions.at[el, 'Score']) 
                ])

        # MOVE TO THE LEFT, WHAT IS UNEXPANDABLE (COMMERCIAL or NOT EXPLORED)
        for el in ini_base:
            for Reaction in el[1]:
                if not self.is_tree_expansion_possible(list(predictions.at[Reaction, 'Next_Iterated'])):
                    el[0].append(Reaction)
                    el[1].remove(Reaction)

        # Convert into DataFrame for being able to sort by Score values:
        a = pd.DataFrame(ini_base, columns=['Route', 'Still_Expandable', 'Solved', 'Score'])

        # ___________________________________________________________________________
        #
        #                      LOOP EXPANDING NODES OVER 
        # ___________________________________________________________________________
        step_count = 0
        while self.something_in_Still_Expandable(list(a[(a['Solved'] != 'Failed') & (a['Solved'] != 'Failed_Loop')]['Still_Expandable'])):
            
            #Check for max step count:
            if step_count >= self.tree_max_steps: 
                if self.log: self.make_single_retropredictions.write_logs('WARNING: tree calculation interupted after max step searches')
                return a            
            step_count += 1
            
            New = []
            for el in range(0, len(a)):
                # If expandable:
                if len(list(a.at[el, 'Still_Expandable'])) > 0 and a.at[el, 'Solved'] != 'Failed_Loop' and a.at[el, 'Solved'] != 'Failed':# and a.at[el, 'Solved'] != 'Yes':
                    over_multiplication = []

                    # For each Reaction present in the Expandables, Multiply the given item, then, multiply item's results to each other:
                    for Reaction in list(a.at[el, 'Still_Expandable']):
                        
                        raw_lists_expand = predictions.at[Reaction, 'Next_Iterated']
                        currated_lists_expand = []

                        # Remove 'Commercial', 'False' or 'Failed':
                        for item in raw_lists_expand:
                            if isinstance(item, list):
                                currated_lists_expand.append(item)

                        # Drop unrealistic reactions before multiplying, especially when the amount of option are too large:
                        filtered_currated_lists_expand = []
                        for each_branch_list in currated_lists_expand:
                            if len(each_branch_list) > self.branching_max_expansion:
                                # Get top branch list:
                                New_branches = []
                                for rxn in each_branch_list:
                                    New_branches.append([rxn, predictions.at[rxn, 'Score']])
                                d = pd.DataFrame(New_branches, columns=['rxn', 'Score']).sort_values('Score', ascending=False).reset_index(drop=True)
                                filtered_currated_lists_expand.append(copy.deepcopy(list(d['rxn'])[0:self.branching_max_expansion]))
                            else:
                                filtered_currated_lists_expand.append(each_branch_list)

                        over_multiplication.append(self.multiply(filtered_currated_lists_expand))
                    final_possibilities_for_reaction_el = self.multiply(over_multiplication)

                    # Append each newly multiplied possibility into the new "New" second column (Still_Expandable), move the previous elements, calculate score:
                    for each_combo in final_possibilities_for_reaction_el:
                        if not isinstance(each_combo, list):
                            each_combo = [each_combo]

                        # Create New Routes, and append previous "Still_Expandable":
                        Route = a.at[el, 'Route'].copy()
                        for each in a.at[el, 'Still_Expandable']:#.copy():
                            Route.append(each)

                        # New solved status, only if all of the new "each_combo" items are all solved
                        New_solved = self.get_solved_status_from_set_of_rxn(list_of_rxn=Route+each_combo, predictions=predictions)
                        
                        # New Score is the previous one multiplied by each element from the new "":
                        New_Score = a.at[el, 'Score'].copy()
                        for exp in each_combo:
                            New_Score *= float(predictions.at[exp, 'Score'])
                            New_Score *= self.step_penalties_rate

                        New.append([Route, each_combo, New_solved, New_Score])

                # If not expandable, append the current row:
                else:
                    New.append([
                        a.at[el, 'Route'],
                        a.at[el, 'Still_Expandable'],
                        a.at[el, 'Solved'],
                        a.at[el, 'Score']
                    ])

            # MOVE TO THE LEFT, WHAT IS UNEXPANDABLE (COMMERCIAL or NOT EXPLORED)
            moved = True
            while moved:
                moved = False
                for el in New:
                    for Reaction in el[1]:
                        if not self.is_tree_expansion_possible(list(predictions.at[Reaction, 'Next_Iterated'])) or el[2] == 'Failed' or el[2] == 'Failed_Loop':
                            moved = True
                            el[0].append(Reaction)
                            el[1].remove(Reaction)

            a = pd.DataFrame(New, columns=['Route', 'Still_Expandable', 'Solved', 'Score'])
            a['Steps'] = [len(x) for x in a['Route']]

            # Drop low scoring fraction:
            top_to_keep = self.tree_min_best
            if top_to_keep < self.tree_only_best*len(a):
                top_to_keep = int(self.tree_only_best*len(a))
            if self.tree_max_best < self.tree_only_best*len(a):
                top_to_keep = self.tree_max_best
                
            solved = a[a['Solved'] == 'Yes'].copy()
            #print(len(solved), 'route solved')  #DEBUG
            unsolved = a[a['Solved'] == 'No'].sort_values('Score', ascending=False).copy()[0:top_to_keep].reset_index(drop=True)
            a = pd.concat([solved, unsolved]).sort_values('Score', ascending=False).copy().reset_index(drop=True)

            # Look for intra-tree loops, allow few re-usages of the same reaction in a given branch:
            for each_tree in range(0, len(a)):
                targets = []
                for each_rxn in a.at[each_tree, 'Route']:
                    targets.append(predictions.at[each_rxn, 'Target'])
                for each_rxn in a.at[each_tree, 'Still_Expandable']:
                    targets.append(predictions.at[each_rxn, 'Target'])
                if len(targets) - len(list(set(targets))) >= 1:
                    a.at[each_tree, 'Solved'] = 'Failed_Loop'

        if self.log: self.make_single_retropredictions.write_logs('Tree calculation done.')
        return a

    def is_further_explorable(self, reaction, predictions):
        ''' Return if a reaction has precursor molecules that are further explorable
        '''
        if 'Failed' in predictions.at[reaction, 'Next_Iterated']:
            return False
        if False in predictions.at[reaction, 'Next_Iterated']:
            return True
        else: return False

    def any_failed(self, reaction_list, predictions):
        for rxn in reaction_list:
            if 'Failed' in predictions.at[rxn, 'Next_Iterated']:
                return True
        return False

    def Get_Best_Scoring_Reactions_From_Tree(self, tree, predictions, proportion=0.3, max_mol=20):
        ''' Get as input the Possibility Tree, 
            Output the unexplored reactions of the top 'proportion' % of the possibilities
        '''
        lenght = len(tree)
        stop_len = int(round(lenght*proportion, 0))
        list_of_list_of_routes = list(tree[tree['Solved'] == 'No'].sort_values('Score', ascending=False)['Route'])[0:stop_len]
        
        conc_list = []
        for route_set in list_of_list_of_routes:
            # Append only if limit not reached, doesn't mean we exactly get to the exact max_mol value since multiple molecules could be appended for one given Route.
            if not self.any_failed(reaction_list=route_set, predictions=predictions):
                if len(conc_list) < max_mol:
                    for reaction in route_set:
                        if self.is_further_explorable(reaction=reaction, predictions=predictions):
                            conc_list.append(reaction)
                            conc_list = copy.deepcopy(list(set(conc_list)))

        return conc_list

    def Get_Best_Precursors_From_Reactions_List(self, Reaction_list, predictions):
        precursors = []
        for reaction in Reaction_list:
            for el in range(0, len(predictions.at[reaction, 'Retro'].split('.'))):
                # Only select molecules that are not already explored:
                if not isinstance(predictions.at[reaction, 'Next_Iterated'][el], list):
                    # Just to be sure it's not a failed retro molecule:
                    if predictions.at[reaction, 'Next_Iterated'][el] == False:
                        single_pre = predictions.at[reaction, 'Retro'].split('.')[el]
                        if single_pre != '' and not self._if_commercial(single_pre): 
                            precursors.append(single_pre)
        return list(set(precursors))

    def init_predictions_dataframe(self, target_cpd):

        predictions, df_prediction_backup = self.make_single_retropredictions.Execute_Retro_Prediction(
            SMILES = target_cpd, 
            mark_count = self.mark_count,
            neighbors = self.neighbors,
            Random_Tagging = self.Random_Tagging,
            AutoTagging = self.AutoTagging, 
            AutoTagging_Beam_Size = self.AutoTagging_Beam_Size, 
            Substructure_Tagging = self.Substructure_Tagging,
            Retro_USPTO = self.Retro_USPTO, 
            Fwd_USPTO_Reag_Pred = self.Fwd_USPTO_Reag_Pred, 
            USPTO_Reag_Beam_Size = self.USPTO_Reag_Beam_Size, 
            similarity_filter = self.similarity_filter,
            confidence_filter = self.confidence_filter,
            Retro_beam_size = self.Retro_beam_size,
            mark_locations_filter = self.mark_locations_filter,
            log = self.log
        )

        predictions = self.solved_status_update_reaction_predictions(predictions=predictions, step=0)
        predictions = self.connect_linked_reaction_predictions(predictions=predictions, next_pred_mol=[target_cpd])

        return predictions, df_prediction_backup

    def _save_checkpoint(self, predictions, tree, final=False):

        if not final:   tmp = 'tmp_'
        else:           tmp = ''

        predictions.to_pickle(  'output/' + self.project_name + '/' + self.log_time_stamp + '__' + tmp + 'prediction.pkl')
        tree.to_pickle(         'output/' + self.project_name + '/' + self.log_time_stamp + '__' + tmp + 'tree.pkl')

    def multistep_search(self, target_cpd, min_solved_routes=20, max_iteration=4, predictions_OLD=''):
        ''' 
            Main function to run the multistep retrosynthesis search
            Input: 
                - target_cpd: SMILES string of the target compound
                - min_solved_routes: minimum number of solved routes to stop the search
                - max_iteration: maximum number of iterations to run the search
                - predictions_OLD: if not empty, the function will continue the search from the previous predictions (Pandas DataFrame)
        '''
        
        if self.log: self.make_single_retropredictions.write_logs('Starting multistep retrosynthesis on: ' + str(target_cpd))

        target_cpd = self.make_single_retropredictions.canonicalize_smiles(target_cpd)
        if self._if_commercial(target_cpd): 
            if self.log: self.make_single_retropredictions.write_logs('Target cpd is already commercial.')

        if self.log: self.make_single_retropredictions.write_logs('Canonical version: ' + str(target_cpd))

        self.Commercial_exclusions.append(target_cpd)

        shift = 0
        if not isinstance(predictions_OLD, pd.DataFrame):
            # Initiate the First Step:
            predictions, _ = self.init_predictions_dataframe(target_cpd=target_cpd)
            tree = self.Get_Tree_From_Reaction_DF(predictions, target_cpd)
            if len(tree[tree['Solved'] == 'Yes']) >= min_solved_routes: 
                if self.log: self.make_single_retropredictions.write_logs('reached min_solved_routes')
                self._save_checkpoint(predictions, tree, final=True)
                return predictions, tree
        else:
            predictions = predictions_OLD.copy()
            tree = self.Get_Tree_From_Reaction_DF(predictions, target_cpd)
            if len(tree[tree['Solved'] == 'Yes']) >= min_solved_routes: 
                if self.log: self.make_single_retropredictions.write_logs('reached min_solved_routes')
                self._save_checkpoint(predictions, tree, final=True)
                return predictions, tree
            shift = max(list(predictions['Step']))

        self._save_checkpoint(predictions, tree, final=False)

        # Iterate over:
        for iteration in range(0+shift, max_iteration+shift):

            if self.log: self.make_single_retropredictions.write_logs('Iteration ' + str(iteration))
            next_pred_rxn = self.Get_Best_Scoring_Reactions_From_Tree(tree, predictions=predictions, proportion=0.35, max_mol=self.max_mol_per_iteration)
            next_pred_mol = self.Get_Best_Precursors_From_Reactions_List(Reaction_list=next_pred_rxn, predictions=predictions)

            if self.log: self.make_single_retropredictions.write_logs('devp rxns: ' + str(next_pred_rxn))
            if self.log: self.make_single_retropredictions.write_logs('devp mols: ' + str(next_pred_mol))

            # Prepare overall current round prediction DataFrame:
            current_preds = predictions[predictions['index'] == 999999999].copy()
            
            # For each selected precursor:
            for mol in tqdm(next_pred_mol):
                # Make prediction
                node_pred, _ = self.make_single_retropredictions.Execute_Retro_Prediction(
                    SMILES = mol, 
                    mark_count = self.mark_count,
                    neighbors = self.neighbors,
                    Random_Tagging = self.Random_Tagging,
                    AutoTagging = self.AutoTagging, 
                    AutoTagging_Beam_Size = self.AutoTagging_Beam_Size, 
                    Substructure_Tagging = self.Substructure_Tagging,
                    Retro_USPTO = self.Retro_USPTO, 
                    Fwd_USPTO_Reag_Pred = self.Fwd_USPTO_Reag_Pred, 
                    USPTO_Reag_Beam_Size = self.USPTO_Reag_Beam_Size, 
                    similarity_filter = self.similarity_filter,
                    confidence_filter = self.confidence_filter,
                    Retro_beam_size = self.Retro_beam_size,
                    mark_locations_filter = self.mark_locations_filter,
                    log = self.log
                    )

                # Calculate stuff:
                node_pred = self.solved_status_update_reaction_predictions(predictions=node_pred, step=iteration+1)
                node_pred['index'] = [i for i in range(len(predictions)+len(current_preds), len(predictions)+len(current_preds)+len(node_pred))]
                node_pred.set_index(node_pred['index'], inplace=True)
                
                # Append node_pred to current_preds:
                current_preds = pd.concat([current_preds, node_pred])

            predictions = pd.concat([predictions, current_preds])
            predictions = self.connect_linked_reaction_predictions(predictions=predictions, next_pred_mol=next_pred_mol)

            self._save_checkpoint(predictions, tree, final=False)
            tree = self.Get_Tree_From_Reaction_DF(predictions, target_cpd)
            if self.log: self.make_single_retropredictions.write_logs('Solved routes: ' + str(len(tree[tree['Solved'] == 'Yes'])))
            
            if len(tree[tree['Solved'] == 'Yes']) >= min_solved_routes: 
                if self.log: self.make_single_retropredictions.write_logs('Reached min_solved_routes')
                self._save_checkpoint(predictions, tree, final=True)
                return predictions, tree

            self._save_checkpoint(predictions, tree, final=False)

        if self.log: self.make_single_retropredictions.write_logs('normal termination')
        self._save_checkpoint(predictions, tree, final=True)
        return predictions, tree
