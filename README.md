# Chemoenzymatic Multistep Retrosynthesis by Disconnection Aware Triple Transformer Loops


This repository complements the _Chemical Science_ article ["Multistep retrosynthesis combining a disconnection aware triple transformer loop with a route penalty score guided tree search"](https://pubs.rsc.org/en/content/articlelanding/2023/sc/d3sc01604h). The work has been extended to Enzymatic reactions as described in the ChemRxiv preprint ["Chemoenzymatic Multistep Retrosynthesis with Transformer Loops"](https://chemrxiv.org/engage/chemrxiv/article-details/6617b48391aefa6ce157c2b4).

The goal of this tool is to predict multistep retrosynthesis routes combining a tree search strategy and Transformer Disconnection-Aware Retrosynthesis models. The single-step retrosynthesis models are augmented by coupling a bias-free systematic tagging as well as a template-based tagging using a reaction center substructure identification from known reactions. The work was extended to biocatalysis by running two TTLs in parallel.

## Setup Environment

Create the conda environment, clone this repo and install:

``` bash
conda create -n MultiStepRetro_ENZ python=3.8.16 -y
conda activate MultiStepRetro_ENZ

git clone https://github.com/reymond-group/MultiStepRetrosynthesisTTL.git

cd MultiStepRetrosynthesisTTL

pip install .

```

## Download USPTO Models

Auto-Tag, Disconnection-Aware Retrosynthesis (T1), Reagent Prediction (T2), and Forward Validation Transformer (T3) models are required and should be placed into your `models` folder and referenced in the config file under parameters: `USPTO_AutoTag_path`, `USPTO_T1_path`, `USPTO_T2_path` and `USPTO_T3_path`.


The USPTO models can be downloaded from Zenodo and placed in the `models` folder with the following commands:

``` bash
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T0_AutoTag_260000.pt?download=1 -O models/USPTO_STEREO_separated_T0_AutoTag_260000.pt
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T1_Retro_255000.pt?download=1 -O models/USPTO_STEREO_separated_T1_Retro_255000.pt
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T2_Reagent_Pred_225000.pt?download=1 -O models/USPTO_STEREO_separated_T2_Reagent_Pred_225000.pt
wget https://zenodo.org/record/8160148/files/USPTO_STEREO_separated_T3_Forward_255000.pt?download=1 -O models/USPTO_STEREO_separated_T3_Forward_255000.pt
```

## ENZR Models

Models can not be made public as trained from data under Reaxys commercial subscription. If you have a Reaction Reaxys API licence, you can download the ENZR enzymatic reaction dataset using `notebook/API_reaxys_ENZR_dataset.py` with your API credentials. The list of reaction IDs used for this work can be found in `notebook/ENZR_Rxn_IDs.txt`. The raw data obtained needs to be converted to SMILES to train the models following the preprocessing described in the [OpenNMT](https://github.com/reymond-group/OpenNMT-py) repository. 


## Commercial Building Blocks

The list of commercial compounds should be requested and downloaded from [MolPort](https://www.molport.com/shop/access-databases) and/or from [Enamine](https://enamine.net/building-blocks/building-blocks-catalog). SMILES should be canonicalized using the same environment and located as one SMILES per line in the `stocks` folder. The file should be referenced in the config file as "commercial_file_path".


## Usage for Multistep Prediction

Edit the `/configs/config_example.yaml` configuration file, change the target compound (`target_cpd`) and tree search parameters as needed.
Start the multistep search as follow:

``` bash
conda activate MultiStepRetro_ENZ
retrosynthesis --config configs/config_example.yaml
```

## Visualizing Results

Results are written into `output/project_name/` as pickle files. Forward validated single step-reaction predictions are stored as `output/project_name/[DateJob]__prediction.pkl`, and full predicted route paths are stored as `output/project_name/[DateJob]__tree.pkl`, which refers to reaction indexes from `[DateJob]__prediction.pkl`. Routes could be sorted by scores to get the best ones. Temporary checkpoints are written in the `output/project_name/` folder after each iteration to monitor the progress of retrosynthesis, but also to resume a job starting from checkpoints. If logs are enabled, those are written into `output/project_name/`.

To visualize predicted routes, check this notebook [`/notebooks/visualize_results.ipynb`](https://github.com/reymond-group/MultiStepRetrosynthesisTTL/blob/main/notebook/visualize_results.ipynb) or the following example:

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

Publication: [OpenNMT: Neural Machine Translation Toolkit](https://arxiv.org/pdf/1805.11462.pdf)

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

## Cite this work

Publication of the original Triple Transformer Loop: [Multistep retrosynthesis combining a disconnection aware triple transformer loop with a route penalty score guided tree search](https://pubs.rsc.org/en/content/articlelanding/2023/sc/d3sc01604h).


``` bash
@article{kreutter_multistep_2023,
	title = {Multistep Retrosynthesis Combining a Disconnection Aware Triple Transformer Loop 
		with a Route Penalty Score Guided Tree Search},
	author = {Kreutter, David and Reymond, Jean-Louis},
	url = {https://pubs.rsc.org/en/content/articlelanding/2023/sc/d3sc01604h},
	doi = {10.1039/D3SC01604H},
	date = {2023-09-20},
	journaltitle = {Chemical Science},
	shortjournal = {Chem. Sci.},
	volume = {14},
	number = {36},
	pages = {9959--9969},
	publisher = {{The Royal Society of Chemistry}},
}
```

Preprint of the Enzymatic work: [Chemoenzymatic Multistep Retrosynthesis with Transformer Loops](https://chemrxiv.org/engage/chemrxiv/article-details/6617b48391aefa6ce157c2b4).

``` bash
@article{kreutterChemoenzymaticMultistepRetrosynthesis2024,
  title = {Chemoenzymatic {{Multistep Retrosynthesis}} with {{Transformer Loops}}},
  author = {Kreutter, David and Reymond, Jean-Louis},
  date = {2024-04-12},
  eprinttype = {ChemRxiv},
  doi = {10.26434/chemrxiv-2024-svr99},
  url = {https://chemrxiv.org/engage/chemrxiv/article-details/6617b48391aefa6ce157c2b4},
  pubstate = {preprint}
}
```