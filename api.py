from chebi_tools.ChEBIStandardizer import ChEBIStandardizer
from chebi_tools.ChEBIGraph import ChEBIGraph

from fastapi import FastAPI
from pydantic import BaseModel
from typing import Optional


class Item(BaseModel):
    name: str
    description: Optional[str] = None
    price: float
    tax: Optional[float] = None


STD = ChEBIStandardizer()
# graph = ChEBIGraph()

app = FastAPI()


@app.get("/")
async def root():
    return "ChEBI tools"


@app.get("/name/{name}/")
async def get_name(name: str):
    print(name)
    suggested_name = STD.suggest_name(name)
    return {"suggested_name": suggested_name}


"""
@app.get("/graph/{name}/")
async def get_subgraph(name:str):
    try:
        suggested_name = graph.get_subgraph(name, show=True)
        return {"suggested_name": suggested_name}
    except KeyError:
        return {"suggested_name": None}
"""
