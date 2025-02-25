from sparrow.cli.run import * 

path_to_config = ''


def run_from_config(path_to_config):
    params = vars(get_args("--config " + path_to_config))

    print('SPARROW will be run with the following parameters:')    
    for k, v in sorted(params.items()):
        if v is not None: 
            print(f"  {k}: {v}")
    print(flush=True)
    output_dir = Path(params['output_dir'])
    output_dir.mkdir(exist_ok=True, parents=True)
    
    save_args(params)

    target_dict, targets = get_target_dict(params['target_csv']) 
    if params['diversity_weight'] > 0: 
        print(f'Overwriting parameter "--cluster" to be similarity for diversity objective')
        params['cluster'] = 'similarity'
    clusters = get_clusters(
        cluster_type=params['cluster'], 
        filepath=params['target_csv'], 
        cutoff=params['cluster_cutoff'], 
        outdir=params['output_dir']
    )
    storage_path = get_path_storage(params, targets)
    selector = build_selector(params, target_dict, storage_path, clusters)
    selector.formulate_and_optimize(extract_routes=not params['no_routes'])

if __name__=='__main__':
    run_from_config(path_to_config)