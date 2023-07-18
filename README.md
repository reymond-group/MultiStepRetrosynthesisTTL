# Multistep Retrosynthesis by a Disconnection Aware Triple Transformer Loop

This repo complements the ChemRxiv preprint ["Multistep retrosynthesis combining a disconnection aware triple transformer loop with a route penalty score guided tree search"](https://chemrxiv.org/engage/chemrxiv/article-details/6422d09a62fecd2a83937199).

The goal of this tool is to predict multistep retrosynthesis routes using a tree search strategy and exploiting the Disconnection-Aware Retrosynthesis model. This single-step retrosynthesis model is augmented by combining a bias-free systematic tagging as well as a template-based-inspired tagging using a reaction center substructure identification from known reactions.

## Setup Environment


``` bash
conda create -n MultiStepRetro python=3.8.16 -y
conda activate MultiStepRetro

git clone https://github.com/reymond-group/MultiStepRetrosynthesisTTL.git
cd MultiStepRetrosynthesisTTL
pip install -e .
```

## Download Models and Tagging Templates

Download models from Zenodo:

``` bash
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T0_AutoTag_260000.pt?download=1 -O models/USPTO_STEREO_separated_T0_AutoTag_260000.pt
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T1_Retro_255000.pt?download=1 -O models/USPTO_STEREO_separated_T1_Retro_255000.pt
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T2_Reagent_Pred_225000.pt?download=1 -O models/USPTO_STEREO_separated_T2_Reagent_Pred_225000.pt
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T3_Forward_255000.pt?download=1 -O models/USPTO_STEREO_separated_T3_Forward_255000.pt

```

## Commercial Building Blocks

The list of commercial compounds should be requested and downloaded from [MolPort](https://www.molport.com/shop/access-databases) and/or from [Enamine](https://enamine.net/building-blocks/building-blocks-catalog). SMILES should be canonicalized using the same environment and located as one SMILES per line and the path. The file should be referenced in the config file as "commercial_file_path".

## Usage for Multistep Prediction

Edit the target molecule SMILES as shown in the default configuration file `/configs/config_example.yaml`, as well as search parameters. Then, start the multistep prediction in a terminal:

``` bash
conda activate MultiStepRetro
retrosynthesis --config configs/config_example.yaml
```

## Visualizing Results

Results are written into `output/project_name/` as pickle files. Forward validated single step-reaction predictions are stored as `output/project_name/DayJob__prediction.pkl`, and full predicted route paths are stored as `output/project_name/DayJob__tree.pkl`, which refers to reaction indexes from `prediction.pkl`. Routes could be sorted by scores to get the best ones. Temporary checkpoints are constantly written in the `output/project_name/` folder after each iteration to monitor the progress of retrosynthesis, it also serves to resume a job starting from a checkpoint. If logs are enabled, those are written into `output/project_name/`.

To visualize predicted routes, check this notebook [`/notebooks/visualize_results.ipynb`](https://github.com/DavKre/TTLRetro/blob/main/notebook/visualize_results.ipynb) or the following example:

``` python
import pandas as pd
import ttlretro.view_routes as vr

project_name = 'config_example'
log_time_stamp = 'YYYY_MM_DD__HHMMSS'

predictions = pd.read_pickle('output/{}/{}__prediction.pkl'.format(project_name, log_time_stamp))
tree = pd.read_pickle('output/{}/{}__tree.pkl'.format(project_name, log_time_stamp))
tree = vr.get_advanced_scores(tree=tree, predictions=predictions)

bests = vr.get_best_first_branches(
    tree=tree, 
    predictions=predictions, 
    num_branches=10, 
    score_metric='Fwd_Conf_Score'
)

vr.display_branch(
    branch_tree_index_or_list_rxn_id=0, 
    tree=bests, 
    predictions=predictions, 
    forwarddirection=True
)

```


## Citations

This repository makes use of existing projects:

### OpenNMT-py

[OpenNMT: Neural Machine Translation Toolkit](https://arxiv.org/pdf/1805.11462.pdf)

[OpenNMT technical report](https://www.aclweb.org/anthology/P17-4012/)

``` bash
@inproceedings{klein-etal-2017-opennmt,
    title = "{O}pen{NMT}: Open-Source Toolkit for Neural Machine Translation",
    author = "Klein, Guillaume  and
      Kim, Yoon  and
      Deng, Yuntian  and
      Senellart, Jean  and
      Rush, Alexander",
    booktitle = "Proceedings of {ACL} 2017, System Demonstrations",
    month = jul,
    year = "2017",
    address = "Vancouver, Canada",
    publisher = "Association for Computational Linguistics",
    url = "https://www.aclweb.org/anthology/P17-4012",
    pages = "67--72",
}
```

### SCScore

Publication: [SCScore: Synthetic Complexity Learned from a Reaction Corpus](https://pubs.acs.org/doi/10.1021/acs.jcim.7b00622)

GitHub repository: [SCScore](https://github.com/connorcoley/scscore)

``` bash
@article{coley_scscore_2018,
	title = {{SCScore}: {Synthetic} {Complexity} {Learned} from a {Reaction} {Corpus}},
	author = {Coley, Connor W. and Rogers, Luke and Green, William H. and Jensen, Klavs F.},
	volume = {58},
	issn = {1549-9596},
	shorttitle = {{SCScore}},
	url = {https://doi.org/10.1021/acs.jcim.7b00622},
	doi = {10.1021/acs.jcim.7b00622},
	number = {2},
	urldate = {2022-09-02},
	journal = {Journal of Chemical Information and Modeling},
	month = feb,
	year = {2018},
	note = {Publisher: American Chemical Society},
	pages = {252--261},
}
```

### Cite this work

``` bash
@article{kreutter_multistep_2023,
	title = {Multistep retrosynthesis combining a disconnection aware triple transformer loop 
        with a route penalty score guided tree search},
	author = {Kreutter, David and Reymond, Jean-Louis},
	url = {https://chemrxiv.org/engage/chemrxiv/article-details/6422d09a62fecd2a83937199},
	doi = {10.26434/chemrxiv-2022-8khth-v2},
	publisher = {ChemRxiv},
	month = mar,
	year = {2023},
}
```
