

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
import pandas as pd

def display_branch(
    tree: pd.DataFrame, 
    predictions: pd.DataFrame, 
    branch_tree_index_or_list_rxn_id, 
    forwarddirection: bool = False, 
    printsmiles:bool = False
) -> None:
    '''
    Display a branch of the tree. 
    If branch_tree_index_or_list_rxn_id is an int, it will display the branch of the tree with that index.
    If branch_tree_index_or_list_rxn_id is a list, it will display the branch of the tree with the rxn_id in the list.
    
    Args:
        tree: The tree dataframe.
        predictions: The predictions dataframe.
        branch_tree_index_or_list_rxn_id: The index of the branch or a list of rxn_id's.
        forwarddirection: If True, it will display the forward direction of the reaction, otherwise it will display the retro direction.
        printsmiles: If True, it will print the smiles of the reactants and products.
    '''
    
    if type(branch_tree_index_or_list_rxn_id) == int:
        route = list(tree.at[branch_tree_index_or_list_rxn_id, 'Route']).copy()
        print('Overall forward confidence score =', round(tree.at[branch_tree_index_or_list_rxn_id, 'Fwd_Conf_Score'], 4))
        print('Overall Guiding RPScore =', round(tree.at[branch_tree_index_or_list_rxn_id, 'Score'], 4))
        print('Overall Penalties =', round(tree.at[branch_tree_index_or_list_rxn_id, 'Simplicity_Score'], 4))
        print('Number of steps =', tree.at[branch_tree_index_or_list_rxn_id, 'Steps'])
        if forwarddirection: print('Forward direction')
        if not forwarddirection: print('WARNING: Retro direction, displayed arrows are not correct! Consider arrows as retrosynthesis arrows.')
        print('_________________________________')
    elif type(branch_tree_index_or_list_rxn_id) == list:
        route = branch_tree_index_or_list_rxn_id
    
    if forwarddirection:
        route.reverse()
    
    for el in route:
        if 'ENZR' in predictions.at[el, 'Forward_Model']:   
            reagents = ''
            enz_rgt = predictions.at[el, 'Reagents'].replace(' ', '').replace('_', ' ')#.replace('_', '')
        else:   
            reagents = predictions.at[el, 'Reagents'] 
            enz_rgt = ''
            
        if forwarddirection: rxn = '{}>{}>{}'.format(predictions.at[el, 'Retro'], reagents, predictions.at[el, 'Target'])
        if not forwarddirection: rxn = '{}>{}>{}'.format(predictions.at[el, 'Target'], reagents, predictions.at[el, 'Retro'])
            
        print('Reaction #{}'.format(round(el, 0)))
        print('Confidence score =', round(predictions.at[el, 'Prob_Forward_Prediction_1'], 3))
        if enz_rgt != '': print('Enzyme(s) = {}'.format(enz_rgt))
        if printsmiles: print(rxn)
        # TODO: Modify the RDkit ReactionToImage to display double retrosynthesis arrows.
        display(Chem.Draw.ReactionToImage(AllChem.ReactionFromSmarts(rxn, useSmiles=True), subImgSize=(400, 210)))
        
    return None

def get_best_first_branches(
    tree: pd.DataFrame,
    predictions: pd.DataFrame,
    num_branches: int = 10,
    score_metric: str = 'Score'
) -> pd.DataFrame:
    '''
    Get the best scoring first branches of the tree. 
    
    Args:
        tree: The tree dataframe.
        predictions: The predictions dataframe.
        num_branches: The max number of branches to return.
        score_metric: The metric to sort the branches by, must be a column of tree (Fwd_Conf_Score, Score, Simplicity_Score).
        
    Returns:
        A dataframe of the best scoring first branches
    '''
    # Check if score_metric is a column of tree:
    assert score_metric in tree.columns, 'score_metric must be a column of tree (Score, Fwd_Conf_Score, Simplicity_Score)'

    tree_copy = tree.copy()
    tree_copy['first_best'] = False

    for each_first in predictions[predictions['Step'] == 0]['index']:
        for el in tree_copy[tree_copy['Solved'] == 'Yes'].sort_values(score_metric, ascending=False).index:
            if each_first in tree_copy.at[el, 'Route']:
                tree_copy.at[el, 'first_best'] = True
                break

    return tree_copy[tree_copy['first_best']].sort_values(score_metric, ascending=False).head(num_branches).drop(columns=['first_best']).reset_index(drop=True)
    
def get_advanced_scores(
    tree: pd.DataFrame,
    predictions: pd.DataFrame
) -> pd.DataFrame:
    '''
    Computes advanced scores for the tree. 
    
    Args:
        tree: The tree dataframe.
        predictions: The predictions dataframe.
    
    Returns:
        - The tree dataframe with the advanced scores added:
            - Fwd_Conf_Score: Overall confidence score of the branch, does not consider molecular simplicity 
              or penalty score. Best sorting metric after the search is complete to get the most reasonable branches.
            - Simplicity_Score: Score without considering the confidence score or penalty score. Useful to get efficient routes that are not 
    '''
    tree['Fwd_Conf_Score'] = 0.0
    for branch in range(0, len(tree['Route'])):
        prod = 1
        for rxn in tree.at[branch, 'Route']:
            prod *= predictions.at[rxn, 'Prob_Forward_Prediction_1']
        tree.at[branch, 'Fwd_Conf_Score'] = prod
    
    tree['Simplicity_Score'] = 0.0
    for branch in range(0, len(tree['Route'])):
        prod = 1
        for rxn in tree.at[branch, 'Route']:
            prod *= (predictions.at[rxn, 'Score'] / predictions.at[rxn, 'Prob_Forward_Prediction_1'])

        tree.at[branch, 'Simplicity_Score'] = prod
    
    return tree
    
