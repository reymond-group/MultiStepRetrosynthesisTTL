# Chemoenzymatic Multistep Retrosynthesis by Disconnection Aware Triple Transformer Loops

This repo complements the ChemRxiv preprint ["Chemoenzymatic Multistep Retrosynthesis with Transformer Loops"](https://chemrxiv.org/engage/chemrxiv/article-details/6617b48391aefa6ce157c2b4). This work is based on the _Chemical Science_ article ["Multistep retrosynthesis combining a disconnection aware triple transformer loop with a route penalty score guided tree search"](https://pubs.rsc.org/en/content/articlelanding/2023/sc/d3sc01604h).

The goal of this tool is to predict multistep retrosynthesis routes using a tree search strategy and exploiting the Disconnection-Aware Retrosynthesis model. This single-step retrosynthesis model is augmented by combining a bias-free systematic tagging as well as a template-based-inspired tagging using a reaction center substructure identification from known reactions. The work was extended to biocatalysis by running two TTL in parallel. 

## Setup Environment


``` bash
conda create -n MultiStepRetro_ENZ python=3.8.16 -y
conda activate MultiStepRetro_ENZ

git clone https://github.com/reymond-group/MultiStepRetrosynthesisTTL.git 
cd TTLRetro
pip install -e .
```

## Model Requirements

Auto-Tag, Disconnection-Aware Retrosynthesis (T1), Reagent Prediction (T2), and Forward Validation Transformer (T3) models should be placed into your `models` folder and referenced in the config file. 

The list of commercial compounds should be canonicalized using the same environment and located as one SMILES per line and the path referenced in the config file.

## Usage for Multistep Prediction

Edit the target molecule SMILES as shown in the default configuration file `/configs/config_example.yaml`, as well as search parameters. Then, start the multistep prediction in a terminal:

``` bash
conda activate MultiStepRetro
retrosynthesis --config configs/config_example.yaml
```

## Visualizing Results

Results are written into `output/project_name/` as pickle files. Forward validated single step-reaction predictions are stored as `output/project_name/DayJob__prediction.pkl`, and full predicted route paths are stored as `output/project_name/DayJob__tree.pkl`, which refers to reaction indexes from `prediction.pkl`. Routes could be sorted by scores to get the best ones. Temporary checkpoints are constantly written in the `output/` folder after each iteration to monitor the progress of retrosynthesis, and also to resume a job starting from those files. If logs are enabled, those are written into `output/project_name/`.

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

## Cite this work

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
