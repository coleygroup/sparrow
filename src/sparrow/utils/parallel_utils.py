""" parallel_utils.py """
import logging

from tqdm import tqdm
from pathos import multiprocessing as mp

# import multiprocess.context as ctx
# ctx._force_start_method('spawn')

def simple_parallel(input_list, function, max_cpu=16,
                    timeout=4000, max_retries=3, desc=None):
    """ Simple parallelization.

    Use map async and retries in case we get odd stalling behavior.

    input_list: Input list to op on
    function: Fn to apply
    max_cpu: Num cpus
    timeout: Length of timeout
    max_retries: Num times to retry this

    """
    from multiprocess.context import TimeoutError

    cpus = min(mp.cpu_count(), max_cpu)
    pool = mp.Pool(processes=cpus)
    async_results = [pool.apply_async(function, args=(i, ))
                     for i in input_list]
    pool.close()

    retries = 0
    while True:
        try:
            list_outputs = []
            for async_result in tqdm(async_results, total=len(input_list), desc=desc):
                result = async_result.get(timeout)
                list_outputs.append(result)

            break
        except TimeoutError:
            retries += 1
            logging.info(f"Timeout Error (s > {timeout})")
            if retries <= max_retries:
                pool = mp.Pool(processes=cpus)
                async_results = [pool.apply_async(function, args=(i, ))
                                 for i in input_list]
                pool.close()
                logging.info(f"Retry attempt: {retries}")
            else:
                raise ValueError()

    return list_outputs


def chunked_parallel(input_list, function, chunks=100,
                     max_cpu=16, timeout=4000, max_retries=3, desc=None):
    """chunked_parallel.

    Args:
        input_list : list of objects to apply function
        function : Callable with 1 input and returning a single value
        chunks: number of hcunks
        max_cpu: Max num cpus
        timeout: Length of timeout
        max_retries: Num times to retry this
    """
    # Adding it here fixes somessetting disrupted elsewhere

    def batch_func(list_inputs):
        outputs = []
        for i in list_inputs: 
            outputs.append(function(i))
        return outputs

    list_len = len(input_list)
    num_chunks = min(list_len, chunks)
    step_size = len(input_list) // num_chunks

    chunked_list = [input_list[i:i+step_size]
                    for i in range(0, len(input_list), step_size)]

    list_outputs = simple_parallel(chunked_list, batch_func, max_cpu=max_cpu,
                                   timeout=timeout, max_retries=max_retries, desc=desc)
    # Unroll
    full_output = [j for i in list_outputs for j in i]

    return full_output