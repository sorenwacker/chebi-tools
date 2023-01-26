import logging
import obonet

import networkx as nx
import pandas as pd

from networkx.exception import NodeNotFound
from pyvis.network import Network
from tqdm import tqdm

from . import ChEBIDownloader


class ChEBIGraph:
    def __init__(self, download_dir=None, G=None):
        self.downloader = ChEBIDownloader(download_dir=download_dir)
        self.downloader.download_missing()
        self.download_dir = self.downloader.download_dir
        self.G = G
        if G is None:
            self.load_graph()

        # self.df = self.get_data_from_graph()

    def check_node(self, node):
        data = self.G.nodes[node]
        if "property_value" not in data.keys():
            return False
        smiles = self.get_property_from_node(node, "smiles")
        mass = self.get_property_from_node(node, "monoisotopicmass")
        if (smiles is None) or ("*" in smiles) or (smiles == ""):
            return False
        if (mass is None) or (mass < 50) or (mass > 1500):
            return False
        return True

    def load_graph(self):
        fn = self.downloader.get_path("obo")
        print(f"Loading topology from: {fn}")
        self.G = obonet.read_obo(fn)

    def remove_unnecessary_nodes(self):
        self.remove_deuterated_compounds()
        self.remove_oligopeptides()
        self.remove_compound_classes_nodes()

    def remove_compound_classes_nodes(self):
        # Remove nodes without propertie values
        print("Removing compound classes nodes")
        nodes_to_remove = [n for n in tqdm(self.nodes) if not self.check_node(n)]
        self.G.remove_nodes_from(nodes_to_remove)

        ## Remove edges which are not in the list
        edge_types_to_keep = [
            "is_conjugate_acid_of",
            "is_conjugate_base_of",
            "is_enantiomer_of",
            "is_tautomer_of",
            "is_a",
        ]
        edges_to_remove = [
            edge for edge in self.G.edges if edge[2] not in edge_types_to_keep
        ]
        self.G.remove_edges_from(edges_to_remove)

    def remove_deuterated_compounds(self):
        print("Removing deuterated compounds")
        self.remove_subgraph("CHEBI:76107")

    def remove_dipeptides(self):
        print("Removing dipeptides")
        self.remove_subgraph("CHEBI:46761")

    def remove_tripeptides(self):
        print("Removing tripeptides")
        self.remove_subgraph("CHEBI:47923")

    def remove_tetrapeptides(self):
        print("Removing tetrapeptides")
        self.remove_subgraph("CHEBI:48030")

    def remove_pentapeptides(self):
        print("Removing pentapeptides")
        self.remove_subgraph("CHEBI:48545")

    def remove_oligopeptides(self):
        print("Removing oligopeptides")
        self.remove_dipeptides()
        self.remove_tripeptides()
        self.remove_tetrapeptides()
        self.remove_pentapeptides()
        self.remove_subgraph("CHEBI:7755")

    #def remove_isotopically_modified_compounds(self):
    #    print('Removing isotopically modified compounds')
    #    self.remove_subgraph("CHEBI:139358", depth=2)

    def remove_subgraph(self, token, depth=1):
        try:
            a = self.get_subgraph(token, depth=depth)
            self.G.remove_nodes_from(a.nodes)
        except NodeNotFound as e:
            logging.warning(e)

    def get_subgraph(
        self,
        token="CHEBI:25350",
        name=None,
        depth=10,
        undirected=True,
        show=False,
        **kwargs,
    ):
        H = nx.ego_graph(self.G, token, depth, undirected=undirected)
        if show:
            if name is None:
                name = token
            return H, self.show_graph(H, name=name, **kwargs)
        return H

    def show_graph(
        self,
        G,
        name="graph",
        height="800px",
        width="1000px",
        notebook=True,
        directed=True,
    ):

        n_nodes = len(G.nodes)

        if n_nodes > 100:
            logging.warning(f"To many nodes to plot (n={n_nodes})")

        fn = f'{"".join([e for e in name if e.isalnum()])}.html'
        for n in G.nodes(data=True):
            n[1]["title"] = n[0]  # add hoovering to graph
            n[1]["label"] = '\n'.join([n[0], n[1]["name"]])  # add hoovering to graph

        nt = Network(height, width, notebook=notebook, directed=directed)
        nt.from_nx(G)
        return nt.show(fn)

    def update_df(self):
        self.df = self.get_data_from_graph()

    def get_data_from_graph(self, props_to_extract=None):
        results = []
        for node in self.G.nodes:
            results.append(
                self.get_data_from_node(node, props_to_extract=props_to_extract)
            )
        return pd.DataFrame.from_records(results).set_index("ChEBI")

    def get_data_from_node(self, token, props_to_extract=None):
        if props_to_extract is None:
            props_to_extract = [
                "monoisotopicmass",
                "charge",
                "formula",
                "inchikey",
                "smiles",
            ]

        data = self.G.nodes[token]

        result = dict(
            ChEBI=token,
            compound_id=int(token.lower().replace("chebi:", "")),
            name=data["name"],
        )

        for prop_name in props_to_extract:
            value = self.get_property_from_node(token, prop_name)
            result.update({prop_name: value})

        return result

    def get_property_from_node(self, token, prop="smiles"):
        properties = self.G.nodes[token]["property_value"]
        value = None
        for string in properties:
            if "/" + prop + " " in string:
                value = string.split('"')[1]
        if value is not None:
            try:
                value = float(value)
            except ValueError:
                pass
        return value

    def get_reference_chebi_of_group(self, group_ids):
        data = self.df.loc[group_ids].copy()
        data["abs_charge"] = data.charge.abs()
        data["name_alpha"] = data.name.str.isalpha()
        data["name_length"] = data.name.apply(len)
        data = data.sort_values("name_length")
        grps = data.groupby(["abs_charge", "name_alpha"])
        for ndx, grp in grps:
            return grp.index[0], grp.name[0]

    def get_group(self, token):
        H = self.get_subgraph(token)
        return list(H.nodes)

    @property
    def nodes(self):
        return self.G.nodes

    @property
    def edges(self):
        return self.G.edges
