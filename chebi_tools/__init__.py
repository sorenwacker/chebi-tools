from .ChEBIStandardizer import ChEBIStandardizer
from .ChEBIDownloader import ChEBIDownloader
from .ChEBIGraph import ChEBIGraph

from pathlib import Path as P

DATA_PATH = P(__file__).parent/'data'
