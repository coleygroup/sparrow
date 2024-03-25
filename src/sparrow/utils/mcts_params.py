"""Parameters for MCTS."""


MCTS_PARAMS = {
  # "smiles": "CN(C)CCOC(c1ccccc1)c1ccccc1",
  "smiles": "",
  "description": "",
  "tags": "",
  "expand_one_options": {
    "template_count": 100,
    "max_cum_template_prob": 0.995,
    "forbidden_molecules": [],
    "known_bad_reactions": [],
    "retro_backend_options": [
      {
        "retro_backend": "template_relevance",
        "retro_model_name": "reaxys",
        "max_num_templates": 2000,
        "max_cum_prob": 0.995,
        "attribute_filter": []
      }
    ],
    "use_fast_filter": True,
    "filter_threshold": 0.75,
    "retro_rerank_backend": "relevance_heuristic",
    "cluster_precursors": False,
    "cluster_setting": {
      "feature": "original",
      "cluster_method": "hdbscan",
      "fp_type": "morgan",
      "fp_length": 512,
      "fp_radius": 1,
      "classification_threshold": 0.2
    },
    "extract_template": False,
    "return_reacting_atoms": True,
    "selectivity_check": False
  },
  "build_tree_options": {
    "expansion_time": 60,
    "max_branching": 25,
    "max_depth": 5,
    "exploration_weight": 1,
    "return_first": False,
    "max_trees": 500,
    "buyable_logic": "and",
    "max_ppg_logic": "none",
    "max_ppg": 0,
    "max_scscore_logic": "none",
    "max_scscore": 0,
    "chemical_property_logic": "none",
    "max_chemprop_c": 0,
    "max_chemprop_n": 0,
    "max_chemprop_o": 0,
    "max_chemprop_h": 0,
    "chemical_popularity_logic": "none",
    "min_chempop_reactants": 5,
    "min_chempop_products": 5,
    "custom_buyables": []
  },
  "enumerate_paths_options": {
    "path_format": "json",
    "json_format": "nodelink",
    "sorting_metric": "plausibility",
    "validate_paths": True,
    "score_trees": False,
    "cluster_trees": False,
    "cluster_method": "hdbscan",
    "min_samples": 5,
    "min_cluster_size": 5,
    "paths_only": False,
    "max_paths": 200
  },
  "run_async": False,
  "result_id": "60407ac1-146d-4939-8cf2-04769c6a7cc5"
}
