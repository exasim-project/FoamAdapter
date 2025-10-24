
from .registry import get_model
from pydantic import BaseModel, Field
from typing import List, Dict

class RegionConfig(BaseModel):
    name: str

class SimulationConfig(BaseModel):
    regions: List[RegionConfig]
    solvers: Dict[str, str]

class RunTime:
    def loop(self):
        yield True

class Domain:
    def __init__(self, name, model):
        self.name = name
        self.model = model

class Simulation:
    def __init__(self, config):
        # Accept dict or SimulationConfig
        if isinstance(config, dict):
            config = SimulationConfig(**config)
        self.config = config
        self.runTime = RunTime()
        self.domains = {}
        for region in config.regions:
            name = region.name
            model_name = config.solvers[name]
            model_cls = get_model(model_name)
            self.domains[name] = Domain(name, model_cls())
