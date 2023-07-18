import sys, getopt
import yaml
import pandas as pd
from multistepretropredictions import Config


def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hc:",["c="])
    except getopt.GetoptError:
        print('run_multistep_retrosynthesis.py -c <configfile.yaml>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('run_multistep_retrosynthesis.py -c <configfile>')
            sys.exit(2)
        elif opt in ("-c", "--ifile"):
            inputfile = arg


    if inputfile != '':
        if not inputfile.endswith('.yaml'):
            inputfile += '.yaml'

        conf_dict = Config('configs/' + inputfile) 
    else:
        print('invalid input config file:', inputfile)
        sys.exit(2)

    # Config Checks:
    if not hasattr(conf_dict, 'target_cpd') or conf_dict.target_cpd == '':
        print('No target selected')
        sys.exit(2)

    if not hasattr(conf_dict, 'project_name'): conf_dict.project_name = ''
    if not hasattr(conf_dict, 'min_solved_routes'): conf_dict.min_solved_routes = 100
    if not hasattr(conf_dict, 'max_iteration'): conf_dict.max_iteration = 10
    if not hasattr(conf_dict, 'mark_count'): conf_dict.mark_count = 3
    if not hasattr(conf_dict, 'neighbors'): conf_dict.neighbors = True
    if not hasattr(conf_dict, 'Random_Tagging'): conf_dict.Random_Tagging = True
    if not hasattr(conf_dict, 'AutoTagging'): conf_dict.AutoTagging = False
    if not hasattr(conf_dict, 'AutoTagging_Beam_Size'): conf_dict.AutoTagging_Beam_Size = 50
    if not hasattr(conf_dict, 'Substructure_Tagging'): conf_dict.Substructure_Tagging = True
    if not hasattr(conf_dict, 'Retro_USPTO'): conf_dict.Retro_USPTO = True
    if not hasattr(conf_dict, 'Std_Fwd_USPTO'): conf_dict.Std_Fwd_USPTO = False
    if not hasattr(conf_dict, 'Fwd_USPTO_Reag_Pred'): conf_dict.Fwd_USPTO_Reag_Pred = True
    if not hasattr(conf_dict, 'USPTO_Reag_Beam_Size'): conf_dict.USPTO_Reag_Beam_Size = 3
    if not hasattr(conf_dict, 'confidence_filter'): conf_dict.confidence_filter = False
    if not hasattr(conf_dict, 'Retro_beam_size'): conf_dict.Retro_beam_size = 3
    if not hasattr(conf_dict, 'max_mol_per_iteration'): conf_dict.max_mol_per_iteration = 30
    if not hasattr(conf_dict, 'mark_locations_filter'): conf_dict.mark_locations_filter = 1
    if not hasattr(conf_dict, 'step_penalties_rate'): conf_dict.step_penalties_rate = 0.8
    if not hasattr(conf_dict, 'tree_only_best'): conf_dict.tree_only_best = 0.50
    if not hasattr(conf_dict, 'tree_min_best'): conf_dict.tree_min_best = 500
    if not hasattr(conf_dict, 'tree_max_best'): conf_dict.tree_max_best = 5000
    if not hasattr(conf_dict, 'branching_max_expansion'): conf_dict.branching_max_expansion = 20
    if not hasattr(conf_dict, 'Commercial_exclusions'): conf_dict.Commercial_exclusions = []
    if not hasattr(conf_dict, 'list_substructures_path'): conf_dict.list_substructures_path = ''
    if not hasattr(conf_dict, 'log'): conf_dict.log = True
    if not hasattr(conf_dict, 'Resume_predictions_from_path'): conf_dict.Resume_predictions_from_path = ''

    if conf_dict.project_name == '':
        conf_dict.project_name = inputfile.replace('.yaml', '')

    if conf_dict.Resume_predictions_from_path != '':
        try:
            conf_dict.predictions_pickle = pd.read_pickle(conf_dict.Resume_predictions_from_path)
        except:
            print('Invalid Resume_predictions_from_path', conf_dict.Resume_predictions_from_path)
            sys.exit(2)
        if not isinstance(conf_dict.predictions_pickle, pd.DataFrame):
            print('Invalid Resume_predictions_from_path', conf_dict.Resume_predictions_from_path)
            sys.exit(2)
    else:
        conf_dict.predictions_pickle = ''

    if not conf_dict.Random_Tagging and not conf_dict.Substructure_Tagging:
        print('Invalid configuration, all tagging modes disabled.')
        sys.exit(2)
    if conf_dict.Substructure_Tagging and conf_dict.list_substructures_path == '':
        print('Invalid configuration, no path to substructures.')
        sys.exit(2)

    # Import Multistep Module:
    from multistepretropredictions import AugmentedMultiStepRetro
    multi_step_retro_predictions = AugmentedMultiStepRetro(
        mark_count = conf_dict.mark_count, 
        neighbors = conf_dict.neighbors, 
        Random_Tagging = conf_dict.Random_Tagging, 
        AutoTagging = conf_dict.AutoTagging, 
        AutoTagging_Beam_Size = conf_dict.AutoTagging_Beam_Size, 
        Substructure_Tagging = conf_dict.Substructure_Tagging, 
        Retro_USPTO = conf_dict.Retro_USPTO, 
        Std_Fwd_USPTO = conf_dict.Std_Fwd_USPTO, 
        Fwd_USPTO_Reag_Pred = conf_dict.Fwd_USPTO_Reag_Pred, 
        USPTO_Reag_Beam_Size = conf_dict.USPTO_Reag_Beam_Size, 
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
        log = conf_dict.log
        )

    if conf_dict.Resume_predictions_from_path != '':
        if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('Resuming process ' + str(inputfile) + ' from ' + str(conf_dict.Resume_predictions_from_path) + ' using following settings:')
    else:
        if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('Running ' + str(inputfile) + ' using following settings:')
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.project_name = ' + str(conf_dict.project_name))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.min_solved_routes = ' + str(conf_dict.min_solved_routes))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.max_iteration = ' + str(conf_dict.max_iteration))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.mark_count = ' + str(conf_dict.mark_count))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.neighbors = ' + str(conf_dict.neighbors))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Random_Tagging = ' + str(conf_dict.Random_Tagging))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.AutoTagging = ' + str(conf_dict.AutoTagging))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.AutoTagging_Beam_Size = ' + str(conf_dict.AutoTagging_Beam_Size))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Substructure_Tagging = ' + str(conf_dict.Substructure_Tagging))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Retro_USPTO = ' + str(conf_dict.Retro_USPTO))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Std_Fwd_USPTO = ' + str(conf_dict.Std_Fwd_USPTO))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Fwd_USPTO_Reag_Pred = ' + str(conf_dict.Fwd_USPTO_Reag_Pred))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.USPTO_Reag_Beam_Size = ' + str(conf_dict.USPTO_Reag_Beam_Size))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.confidence_filter = ' + str(conf_dict.confidence_filter))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Retro_beam_size = ' + str(conf_dict.Retro_beam_size))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.max_mol_per_iteration = ' + str(conf_dict.max_mol_per_iteration))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.mark_locations_filter = ' + str(conf_dict.mark_locations_filter))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.step_penalties_rate = ' + str(conf_dict.step_penalties_rate))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.tree_only_best = ' + str(conf_dict.tree_only_best))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.tree_min_best = ' + str(conf_dict.tree_min_best))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.tree_max_best = ' + str(conf_dict.tree_max_best))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.branching_max_expansion = ' + str(conf_dict.branching_max_expansion))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Commercial_exclusions = ' + str(conf_dict.Commercial_exclusions))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.list_substructures_path = ' + str(conf_dict.list_substructures_path))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.log = ' + str(conf_dict.log))
    if multi_step_retro_predictions.log: multi_step_retro_predictions.make_single_retropredictions.write_logs('conf_dict.Resume_predictions_from_path = ' + str(conf_dict.Resume_predictions_from_path))

    # RUN PREDICTIONS:
    predictions, tree = multi_step_retro_predictions.multistep_search(
        target_cpd = conf_dict.target_cpd, 
        min_solved_routes = conf_dict.min_solved_routes,
        max_iteration = conf_dict.max_iteration, 
        predictions_OLD = conf_dict.predictions_pickle
        )


if __name__ == '__main__':
    main(sys.argv[1:])

