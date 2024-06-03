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

    target_dict, targets, clusters = get_target_dict(params['target_csv'], clusters=params['custom_cluster']) 
    storage_path = get_path_storage(params, targets)
    selector = build_selector(params, target_dict, storage_path, clusters)
    selector = optimize(selector, params)

if __name__=='__main__':
    run_from_config(path_to_config)