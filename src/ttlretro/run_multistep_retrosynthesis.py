import sys
from os.path import exists
import pandas as pd
import argparse
from ttlretro.config_yaml import Config
from ttlretro import MultiStepGraphRetro

# Usage: python run_multistep_retrosynthesis.py -c <configfile>



# Dictionnary of default values:
default_values = {
    'project_name': '',
    'min_solved_routes': 100,
    'max_iteration': 10,
    'mark_count': 3,
    'neighbors': True,
    'Random_Tagging': True,
    'AutoTagging': False,
    'AutoTagging_Beam_Size': 50,
    'Substructure_Tagging': True,
    'Retro_USPTO': True,
    'Fwd_USPTO_Reag_Pred': True,
    'USPTO_Reag_Beam_Size': 3,
    'similarity_filter': False,
    'confidence_filter': False,
    'Retro_beam_size': 3,
    'max_mol_per_iteration': 30,
    'mark_locations_filter': 1,
    'step_penalties_rate': 0.8,
    'tree_only_best': 0.50,
    'tree_min_best': 500,
    'tree_max_best': 5000,
    'branching_max_expansion': 20,
    'Commercial_exclusions': [],
    'list_substructures_path': '',
    'log': True,
    'Resume_predictions_from_path': '', 
    'Model_Folder': '',
    'USPTO_AutoTag_path': '',
    'USPTO_T1_path': '', 
    'USPTO_T2_path': '',
    'USPTO_T3_path': '',
    'tmp_file_path': 'tmp/',
    'commercial_file_path': 'stocks/Commercial_canonical.smi', 
    }

def _load_check_config(configfile):
    
    # Input file check:
    if exists(configfile):
        conf_dict = Config(configfile) 
    elif exists(configfile + '.yaml'):
        conf_dict = Config(configfile + '.yaml')
    else:
        print('Cannot find input config file: {}'.format(configfile))
        sys.exit(2)
        
    if not exists(conf_dict.commercial_file_path):
        print('Cannot find commercially available SMILES file: {}'.format(conf_dict.commercial_file_path))
        sys.exit(2)
            
    # Config Checks:
    if not hasattr(conf_dict, 'target_cpd') or conf_dict.target_cpd == '':
        print('No target selected in {}, please assign a SMILES.'.format(configfile))
        sys.exit(2)
        
    return conf_dict
    
def _load_values(conf_dict, configfile):
    
    # Assign default values if not specified in the input file:
    for key, value in default_values.items():
        if not hasattr(conf_dict, key):
            setattr(conf_dict, key, value)
    
    # Assign project output folder name if not specified in the input file as the name of the config file:
    if conf_dict.project_name == '':
        conf_dict.project_name = configfile.replace('.yaml', '').split('/')[-1]
        
    return conf_dict

def _load_previous_preds(conf_dict):
    
    if conf_dict.Resume_predictions_from_path != '':
        try:
            conf_dict.predictions_pickle = pd.read_pickle(conf_dict.Resume_predictions_from_path)
        except:
            print('Invalid Resume_predictions_from_path, cannot open:', conf_dict.Resume_predictions_from_path)
            sys.exit(2)
            
        if not isinstance(conf_dict.predictions_pickle, pd.DataFrame):
            print('Invalid Resume_predictions_from_path, not a Pandas DataFrame:', conf_dict.Resume_predictions_from_path)
            sys.exit(2)
    else:
        conf_dict.predictions_pickle = ''
        
    return conf_dict

def _check_retro_config(conf_dict):
    if not conf_dict.Random_Tagging and not conf_dict.Substructure_Tagging:
        print('Invalid configuration, all tagging modes disabled.')
        sys.exit(2)
    if conf_dict.Substructure_Tagging and conf_dict.list_substructures_path == '':
        print('Invalid configuration, no path to substructures.')
        sys.exit(2)
    
    return conf_dict

def _load_multistep_graph_retro(conf_dict):
    
    return MultiStepGraphRetro(
        mark_count = conf_dict.mark_count, 
        neighbors = conf_dict.neighbors, 
        Random_Tagging = conf_dict.Random_Tagging, 
        AutoTagging = conf_dict.AutoTagging, 
        AutoTagging_Beam_Size = conf_dict.AutoTagging_Beam_Size, 
        Substructure_Tagging = conf_dict.Substructure_Tagging, 
        Retro_USPTO = conf_dict.Retro_USPTO, 
        Fwd_USPTO_Reag_Pred = conf_dict.Fwd_USPTO_Reag_Pred, 
        USPTO_Reag_Beam_Size = conf_dict.USPTO_Reag_Beam_Size, 
        similarity_filter = conf_dict.similarity_filter, 
        confidence_filter = conf_dict.confidence_filter, 
        Retro_beam_size = conf_dict.Retro_beam_size, 
        max_mol_per_iteration = conf_dict.max_mol_per_iteration, 
        mark_locations_filter = conf_dict.mark_locations_filter, 
        step_penalties_rate = conf_dict.step_penalties_rate, 
        tree_only_best = conf_dict.tree_only_best, 
        tree_min_best = conf_dict.tree_min_best, 
        tree_max_best = conf_dict.tree_max_best, 
        branching_max_expansion = conf_dict.branching_max_expansion, 
        project_name = conf_dict.project_name, 
        Commercial_exclusions = conf_dict.Commercial_exclusions, 
        list_substructures_path = conf_dict.list_substructures_path, 
        log = conf_dict.log, 
        USPTO_AutoTag_path = conf_dict.Model_Folder + conf_dict.USPTO_AutoTag_path,
        USPTO_T1_path = conf_dict.Model_Folder + conf_dict.USPTO_T1_path, 
        USPTO_T2_path = conf_dict.Model_Folder + conf_dict.USPTO_T2_path,
        USPTO_T3_path = conf_dict.Model_Folder + conf_dict.USPTO_T3_path,
        tmp_file_path = conf_dict.tmp_file_path, 
        commercial_file_path = conf_dict.commercial_file_path,
        )

def _write_logs_before_start(multi_step_retro_predictions, conf_dict, configfile):
    
    if multi_step_retro_predictions.log:
        if conf_dict.Resume_predictions_from_path != '':    multi_step_retro_predictions.make_single_retropredictions.write_logs('Resuming process {} from {} using following settings:'.format(configfile, conf_dict.Resume_predictions_from_path))
        else:                                               multi_step_retro_predictions.make_single_retropredictions.write_logs('Running {} using following settings:'.format(configfile))
        
        for key in dir(conf_dict):
            if not key.startswith('__') and not key.endswith('__'):
                multi_step_retro_predictions.make_single_retropredictions.write_logs('{} = {}'.format(key, getattr(conf_dict, key)))
    
    return None
    
def run_retrosynthesis(configfile):
    """
    Run multistep retrosynthesis from a configuration file.
    """
    
    conf_dict = _load_check_config(configfile)
    conf_dict = _load_values(conf_dict, configfile)
    conf_dict = _load_previous_preds(conf_dict)
    conf_dict = _check_retro_config(conf_dict)

    # Initiate Multistep Module:
    multi_step_retro_predictions = _load_multistep_graph_retro(conf_dict)

    _write_logs_before_start(multi_step_retro_predictions, conf_dict, configfile)
    
    # RUN PREDICTIONS:
    _, _ = multi_step_retro_predictions.multistep_search(
        target_cpd = conf_dict.target_cpd, 
        min_solved_routes = conf_dict.min_solved_routes, 
        max_iteration = conf_dict.max_iteration, 
        predictions_OLD = conf_dict.predictions_pickle
        )

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', help='Path to the configuration file')
    args = parser.parse_args()

    run_retrosynthesis(args.config)

if __name__ == '__main__':
    main()
