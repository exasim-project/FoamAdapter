_MODEL_REGISTRY = {}

def register_model(name, cls):
    _MODEL_REGISTRY[name] = cls

def get_model(name):
    return _MODEL_REGISTRY[name]

def list_models():
    return list(_MODEL_REGISTRY.keys())
