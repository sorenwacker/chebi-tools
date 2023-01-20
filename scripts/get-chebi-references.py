import logging

import pandas as pd
import numpy as np

from chebi_tools.ChEBIGraph import ChEBIGraph



def main():

    graph = ChEBIGraph()
    graph.remove_unnecessary_nodes()
    graph.update_df()

    all_nodes = graph.df.index.to_list()
    N_all_nodes = len(all_nodes)

    results = pd.DataFrame(index=all_nodes, columns=['ref_chebi', 'ref_name'])
    n = 0

    errors = []

    while len(all_nodes) > 0: 
        node = all_nodes[0]
        current_n = len(all_nodes)
        percent = np.round(100*(N_all_nodes-current_n)/N_all_nodes,2)
        
        print(f'{n}: {current_n} {percent}%', end="\r")

        group = graph.get_group(node)
        group_size = len(group)
        
        group_is_large_size = 10
        if group_size > group_is_large_size:
            logging.warning(f'Group {node} size {group_size} > {group_is_large_size}')
            errors.append(node)
            continue
            
        ref_id, ref_name = graph.get_reference_chebi_of_group(group)
        
        results.loc[group, ['ref_chebi', 'ref_name']] = ref_id, ref_name
        [all_nodes.remove(e) for e in group]
        n += 1
        if n%100 == 0:
            results.to_parquet('reference-chebis.parquet')
        
    results.to_parquet('reference-chebis.parquet')

    with open(r'errors.log', 'w') as fp:
        for e in errors:
            fp.write(f"{e}\n")
    
    print('Done')


if __name__ == '__main__':
    main()
