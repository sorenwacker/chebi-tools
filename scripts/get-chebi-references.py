import shutil
import logging

import pandas as pd
import numpy as np

from pathlib import Path as P

from chebi_tools.ChEBIGraph import ChEBIGraph



def main():

    path = P('data')
    path.mkdir(parents=True, exist_ok=True)
    fn_out = path/"reference-chebis.parquet"
    print(f'Output file: {fn_out}')

    chebi_id_to_int = lambda x: int(x.lower().replace('chebi:', ''))

    # Prepare Graph and dataframe
    graph = ChEBIGraph()
    graph.remove_unnecessary_nodes()
    graph.update_df()

    # List of all nodes to work through
    all_nodes = graph.df.index.to_list()
    N_all_nodes = len(all_nodes)

    # Dataframe to fill in results
    results = pd.DataFrame(index=[chebi_id_to_int(e) for e in all_nodes], columns=["CHEBI", "REF_CHEBI", "REF_NAME"])
    results.index.name = "CHEBI"

    n = 0
    errors = []
    while len(all_nodes) > 0:
        # Current note to work find group
        node = all_nodes[0]

        # Progress print
        current_n = len(all_nodes)
        percent = np.round(100 * (N_all_nodes - current_n) / N_all_nodes, 2)
        print(f"{n}: {current_n} {percent}%", end="\r")

        # Get the reference group and its size
        group = graph.get_group(node)
        group_size = len(group)

        # Log large groups for future reference
        group_is_large_size = 10
        if group_size > group_is_large_size:
            logging.warning(f"Group {node} size {group_size} > {group_is_large_size}")
            errors.append(node)
            with open("errors.log", "a") as f:
                f.write(f"{node}\n")

        # Get reference data from graph
        ref_id, ref_name = graph.get_reference_chebi_of_group(group)

        # Add data to results table
        ndx = [chebi_id_to_int(e) for e in group]
        results.loc[ndx, ["REF_CHEBI", "REF_NAME"]] = ref_id, ref_name
        results.loc[ndx, "CHEBI"] = group
        
        # Remove all nodes from the to-work-on list
        [all_nodes.remove(e) for e in group]

        # Every 100th step write temporary data
        n += 1
        if n % 100 == 0:
            results.to_parquet(path/"tmp_reference-chebis.parquet")
        
    # Format dataframe
    final_results = results.dropna()
    
    # Write to disk
    final_results.to_parquet(fn_out)
    print("Done")


if __name__ == "__main__":
    main()
