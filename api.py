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

std = ChEBIStandardizer()
graph = ChEBIGraph()

app = FastAPI()

@app.get("/")
async def root():
    return {"message": "ChEBI tool"}

@app.get("/name/{name}/")
async def get_name(name:str):
    print(name)
    suggested_name = std.suggest_name(name)
    return {"suggested_name": suggested_name}


@app.get("/graph/{name}/")
async def get_subgraph(name:str):
    graph.get_subgraph(name, show=True)
    return {"suggested_name": suggested_name}
