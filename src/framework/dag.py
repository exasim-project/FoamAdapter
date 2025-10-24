def build_global_dag(domains):
    """
    Build a global DAG from multiple domain models, supporting interdomain dependencies.
    Each node is named as 'domain.step'.
    """
    G = nx.DiGraph()
    for domain_name, domain in domains.items():
        for step in domain.__class__._steps:
            node_name = f"{domain_name}.{step.name}"
            G.add_node(node_name, meta=step)
            for dep in step.depends_on:
                G.add_edge(dep, node_name)
    return G
import networkx as nx

def build_step_dag(model):
    steps = getattr(model.__class__, "_steps", [])
    G = nx.DiGraph()
    # Add nodes
    for step in steps:
        G.add_node(step.name, meta=step)
    # Add edges by depends_on
    for step in steps:
        for dep in step.depends_on:
            G.add_edge(dep, step.name)
    # If no depends_on, use order for linear chain
    if not any(s.depends_on for s in steps):
        steps_sorted = sorted(steps, key=lambda s: s.order)
        for i in range(1, len(steps_sorted)):
            G.add_edge(steps_sorted[i-1].name, steps_sorted[i].name)
    return G