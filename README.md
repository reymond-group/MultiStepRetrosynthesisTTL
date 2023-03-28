# Multi-Step Retrosynthesis by a Disconnection Aware Triple Transformer Loop

The goal of this tool is to predict multi-step retrosynthesis routes using a tree search strategy and exploiting the Disconnection-Aware Retrosynthesis model. This single-step retrosynthesis model is augmented by combining a bias-free systematic tagging as well as a template-based-inspired tagging using a reaction center substructure identification from known reactions, and an AutoTag Transformer model.

## Environment:

``` console
conda env create --name graphretrosynthesis --file environment.yml
conda activate graphretrosynthesis
```

## Requirements

Disconnection-Aware Retrosynthesis, Reagent Prediction and Forward Validation Transformer models should be placed into `/makesingleretropredictions/onmt/models`. Names of those models should be changed accordingly in `/makesingleretropredictions/make_single_retro_prediction.py`.

Commercial compounds should be canonicalized using the same environment and located as one SMILES per line in `/stocks/Commercial_canonical.smi`.

## Usage for Multi-Step Prediction

Edit the target molecule SMILES such as shown in the default configuration file `/configs/config_example.yaml`, then run the multistep prediction as:

``` console
python run_multistep_retrosynthesis.py -c config_example.yaml
```

## Results

Results are written into `output/project_name/` as pickle files. Forward validated single step-reaction predictions are stored as `output/project_name/DayJob__prediction.pkl`, full predicted route paths are stored as `output/project_name/DayJob__tree.pkl` which refers to reaction indexes from `prediction.pkl`. Routes could be sorted by scores to get the best ones. Temporary files are constantly written after each iteration to monitor the progress of retrosynthesis, and also to resume a job starting from those files. If enabled, logs are written into `output/project_name/`.

```python
import pandas as pd

predictions = pd.read_pickle('output/project_name/Date__prediction.pkl')
tree = pd.read_pickle('output/project_name/Date__tree.pkl')

# Print best routes:
print(tree[tree['Solved'] == 'Yes'].sort_values('Score', ascending=False).head(5))

# Print all SMILES steps of a given route and solved status for each molecule:
route = 0
for rxn in tree.at[route, 'Route']:
	print(predictions.at[rxn, 'Retro'], predictions.at[rxn, 'Solved'])
```




