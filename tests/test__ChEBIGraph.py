from chebi_tools import ChEBIGraph


graph = ChEBIGraph()
graph.update_df()
graph.remove_unnecessary_nodes()
graph.update_df()



def test__creatine():
    group = graph.get_subgraph('CHEBI:16919', show=False)
    ref_id, ref_name = graph.get_reference_chebi_of_group(group)
    assert ref_name == 'creatine', ref_name


def test__mevalonate():
    group = graph.get_subgraph('CHEBI:25351', show=False)
    ref_id, ref_name = graph.get_reference_chebi_of_group(group)
    assert ref_name == 'mevalonic acid', ref_name


def test__arabinitol():
    group = graph.get_subgraph('CHEBI:18333', show=False)
    ref_id, ref_name = graph.get_reference_chebi_of_group(group)
    assert ref_name == 'arabinitol', ref_name


def test__aspartic_acid():
    group = graph.get_subgraph('CHEBI:22660', show=False)
    ref_id, ref_name = graph.get_reference_chebi_of_group(group)
    assert ref_name == 'aspartic acid', ref_name


def test__ornithine():
    group = graph.get_subgraph('CHEBI:15729', show=False)
    ref_id, ref_name = graph.get_reference_chebi_of_group(group)
    assert ref_name == 'ornithine', ref_name



