# chebi-tools

Tools to extract subgraphs from the ChEBI database and for standardization of molecular entities e.g. for applications in metabolomics.


## ChEBIStandardizer

    from chebi_tools import ChEBIStandardizer
    std.process('mevalonic acid')
    # >>> {'CHEBI': 'CHEBI:25350', 'REF_CHEBI': 'CHEBI:25351','REF_NAME': 'mevalonic acid', 'EXACT_MATCH': True, 'COMPOUND_ID': 25350, 'SMILES': 'CC(O)(CCO)CC([O-])=O', 'QUERY': 'mevalonate'}
    
    std.process('mevalonate')
    # >>> {'CHEBI': 'CHEBI:17710', 'REF_CHEBI': 'CHEBI:25351', 'REF_NAME': 'mevalonic acid', 'EXACT_MATCH': True, 'COMPOUND_ID': 17710, 'SMILES': 'C[C@@](O)(CCO)CC(O)=O', 'QUERY': 'mevalonic acid'}
   
   
## ChEBIGraph

    from chebi_tools import ChEBIStandardizer
    graph.remove_unnecessary_nodes()
    graph.update_df()

### Mevalonic acid

    a, fig = graph.get_subgraph('CHEBI:25351', show=True)
    fig

![image](https://user-images.githubusercontent.com/3391614/216475726-f89e211c-bc4e-4288-a670-5415852ed1ed.png)


### Aspartic acid

    a, fig = graph.get_subgraph('CHEBI:22660', show=True)
    fig

![image](https://user-images.githubusercontent.com/3391614/216479405-9824c30d-dcf7-4ae9-9973-daccf0744111.png)
