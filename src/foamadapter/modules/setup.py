#!/usr/bin/env python3
"""
Validate and visualize a dependency DAG using graphlib + networkx.

- Uses graphlib.TopologicalSorter for cycle-safe validation and ordering.
- Uses networkx.multipartite_layout for a clean layered diagram.
- Sets a per-node 'subset' attribute (topological level) required by multipartite_layout.
"""

from __future__ import annotations

from graphlib import TopologicalSorter, CycleError
import networkx as nx
import matplotlib.pyplot as plt
from typing import Dict, Set, Iterable, List, Tuple, Any, Optional
from pathlib import Path
from foamadapter.modules.fields import Field
from foamadapter.modules.models import Models

# ----------------------------
# Validation (graphlib)
# ----------------------------
def validate_with_graphlib(deps: Dict[str, Set[str]]) -> List[str]:
    """Validate the dependency dict and return a topological order. Raises CycleError if a cycle exists."""
    ts = TopologicalSorter(deps)
    order = list(ts.static_order())
    return order


# ----------------------------
# Build networkx graph
# ----------------------------
def build_nx_graph(deps: Dict[str, Set[str]]) -> nx.DiGraph:
    """Convert the dependency mapping into a networkx.DiGraph with edges dep -> node."""
    G = nx.DiGraph()
    for node, dep_set in deps.items():
        if node not in G:
            G.add_node(node)
        for dep in dep_set:
            if dep not in G:
                G.add_node(dep)
            G.add_edge(dep, node)
    return G


# ----------------------------
# Compute layered positions (TOP -> BOTTOM)
# ----------------------------
def layered_positions(
    G: nx.DiGraph,
    *,
    top_to_bottom: bool = True,
) -> Tuple[Dict[Any, Tuple[float, float]], Dict[Any, int], List[Set[str]]]:
    """
    Compute layered positions using networkx.multipartite_layout.
    - Assigns node attribute 'subset' = topological generation index.
    - Uses align='vertical' so subsets are stacked along Y (levels).
    - Ensures generation 0 is at the TOP if top_to_bottom=True.
    Returns: (pos, subset_by_node, layers)
    """
    if not nx.is_directed_acyclic_graph(G):
        raise nx.NetworkXUnfeasible("Graph is not a DAG (cycle detected).")

    layers = list(nx.topological_generations(G))  # [set(nodes), ...]
    subset_by_node = {node: i for i, layer in enumerate(layers) for node in layer}
    nx.set_node_attributes(G, subset_by_node, name="subset")

    # Stack along Y (vertical): nodes of the same subset share X? No—here they share Y (different Y per subset).
    pos = nx.multipartite_layout(G, subset_key="subset", align="horizontal")

    if top_to_bottom:
        # Make sure layer 0 is visually at the TOP (largest y). If not, flip Y.
        layer_y = {i: sum(pos[n][1] for n in layer) / max(1, len(layer)) for i, layer in enumerate(layers)}
        top_layer_index_by_y = max(layer_y, key=lambda k: layer_y[k])  # layer with highest y
        if top_layer_index_by_y != 0:
            for n, (x, y) in pos.items():
                pos[n] = (x, -y)

    return pos, subset_by_node, layers


# ----------------------------
# Level guide lines
# ----------------------------
def draw_level_guides(
    ax,
    pos: Dict[Any, Tuple[float, float]],
    subset_by_node: Dict[Any, int],
    *,
    label_prefix: str = "Level ",
):
    """Draw faint horizontal guide lines for each unique Y (level) with a left-side label."""
    ys_by_level: Dict[int, List[float]] = {}
    for n, lvl in subset_by_node.items():
        ys_by_level.setdefault(lvl, []).append(pos[n][1])

    level_y = {lvl: sum(ys) / len(ys) for lvl, ys in ys_by_level.items()}

    xs = [x for x, _ in pos.values()]
    xmin, xmax = min(xs), max(xs)
    xpad = 0.15 * (xmax - xmin if xmax > xmin else 1.0)
    xmin -= xpad
    xmax += xpad

    for lvl, y in sorted(level_y.items()):
        ax.hlines(y, xmin, xmax, linestyles="dotted", linewidth=1.0, alpha=0.6)
        ax.text(xmin, y, f"{label_prefix}{lvl}", va="center", ha="right", fontsize=9, alpha=0.8)


# ----------------------------
# Visualization
# ----------------------------
def visualize_dag(
    deps: Dict[str, Set[str]],
    title: str = "Dependency Graph",
    filename: Optional[str | Path] = None,
    show: bool = True,
    *,
    top_to_bottom: bool = True,
    level_lines: bool = True,
    special_node: Optional[str] = "turbulenceModel",
) -> nx.DiGraph:
    """
    Validate with graphlib, build a networkx graph, compute layered layout,
    and draw the DAG. Saves to filename if provided.
    - top_to_bottom: place source layer at the top.
    - level_lines: draw dotted guides per level.
    - special_node: name to draw with distinct styling (None to disable).
    """
    try:
        order = validate_with_graphlib(deps)
        print("✅ Valid DAG. Topological order:", order)
    except CycleError as e:
        print("❌ Cycle detected:", e.args)
        raise

    G = build_nx_graph(deps)
    pos, subset_by_node, layers = layered_positions(G, top_to_bottom=top_to_bottom)

    # Create figure/axes
    fig, ax = plt.subplots(figsize=(10, 8))

    # Optional level guides
    if level_lines:
        draw_level_guides(ax, pos, subset_by_node)

    # Prepare node partitions
    all_nodes = list(G.nodes)
    specials = [special_node] if (special_node and special_node in G) else []
    others = [n for n in all_nodes if n not in specials]

    # Draw other nodes
    nx.draw_networkx_nodes(
        G, pos, nodelist=others,
        node_size=3800, node_color="#EEEEEE", edgecolors="black", ax=ax
    )

    # Draw special node(s) with distinct representation
    if specials:
        nx.draw_networkx_nodes(
            G, pos, nodelist=specials,
            node_size=4600, node_color="#FFE8B3",
            edgecolors="black", linewidths=2.0,
            node_shape="s",  # square (try: "D","^","v","p","h","o")
            ax=ax
        )

    # Labels (bold for specials)
    nx.draw_networkx_labels(G, pos, labels={n: n for n in others}, font_size=9, ax=ax)
    if specials:
        for s in specials:
            nx.draw_networkx_labels(G, pos, labels={s: s}, font_size=9, font_weight="bold", ax=ax)

    # Edges: dashed incoming edges to special, solid for others
    if specials:
        incoming_to_special = [(u, v) for (u, v) in G.edges if v in specials]
        other_edges = [(u, v) for (u, v) in G.edges if v not in specials]
    else:
        incoming_to_special = []
        other_edges = list(G.edges)

    nx.draw_networkx_edges(
        G, pos, edgelist=other_edges,
        arrows=True, arrowstyle='-|>', arrowsize=18, width=1.2,
        edge_color="black", min_source_margin=0, min_target_margin=24, ax=ax
    )
    if incoming_to_special:
        nx.draw_networkx_edges(
            G, pos, edgelist=incoming_to_special,
            arrows=True, arrowstyle='-|>', arrowsize=20, width=1.6,
            style="dashed", edge_color="black",
            min_source_margin=0, min_target_margin=28, ax=ax
        )

    ax.set_title(title)
    ax.axis("off")
    fig.tight_layout()

    if filename:
        filename = Path(filename)
        filename.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(filename, dpi=160, bbox_inches="tight")
        print(f"💾 Saved DAG to: {filename.resolve()}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return G



def initialize(fields: Fields, models: Models) -> None:
    # """
    # Topologically initialize fields and models in-place.

    # fields.entries / models.entries can hold either:
    #   - ready-built objects
    #   - factory instances (functions or classes with `.dependencies`)

    # Each factory is called with an Env containing the current
    # fields/models state.
    # """
    # def key(kind: str, name: str) -> str:
    #     return f"{kind}:{name}"

    # # --- Step 1: Build the full dependency DAG ---
    # dag: Dict[str, set[str]] = {}
    # is_factory: Dict[str, bool] = {}

    # # Helper: resolve a dependency reference (short name → namespaced key)
    # def resolve_ref(owner: str, ref: str) -> str:
    #     if ":" in ref:
    #         kind, n = ref.split(":", 1)
    #         kind, n = kind.strip(), n.strip()
    #         if kind not in ("field", "model"):
    #             raise ValueError(f"{owner}: invalid prefix '{kind}' in '{ref}'")
    #         return key(kind, n)
    #     # short name: disambiguate between field/model
    #     in_f, in_m = ref in fields.entries, ref in models.entries
    #     if not in_f and not in_m:
    #         raise KeyError(f"{owner}: unknown dependency '{ref}'")
    #     if in_f and in_m:
    #         raise KeyError(f"{owner}: ambiguous '{ref}' (exists in fields and models)")
    #     return key("field" if in_f else "model", ref)

    # # Fields
    # for name, entry in fields.entries.items():
    #     k = key("field", name)
    #     if callable(entry) and hasattr(entry, "dependencies"):
    #         deps = getattr(entry, "dependencies")
    #         dag[k] = {resolve_ref(k, d) for d in deps}
    #         is_factory[k] = True
    #     else:
    #         dag[k] = set()
    #         is_factory[k] = False

    # # Models
    # for name, entry in models.entries.items():
    #     k = key("model", name)
    #     if callable(entry) and hasattr(entry, "dependencies"):
    #         deps = getattr(entry, "dependencies")
    #         dag[k] = {resolve_ref(k, d) for d in deps}
    #         is_factory[k] = True
    #     else:
    #         dag[k] = set()
    #         is_factory[k] = False

    # # --- Step 2: Topological sort ---
    # try:
    #     order = list(TopologicalSorter(dag).static_order())
    # except CycleError as e:
    #     raise RuntimeError(f"Dependency cycle detected: {e}") from e

    # # --- Step 3: Build in order ---
    # built_fields: Dict[str, Field] = {
    #     k: v for k, v in fields.entries.items()
    #     if not (callable(v) and hasattr(v, "dependencies"))
    # }
    # built_models: Dict[str, Model] = {
    #     k: v for k, v in models.entries.items()
    #     if not (callable(v) and hasattr(v, "dependencies"))
    # }

    # for ref in order:
    #     kind, name = ref.split(":", 1)
    #     if not is_factory[ref]:
    #         continue

    #     env = Env(built_fields, built_models)

    #     if kind == "field":
    #         factory = fields.entries[name]
    #         built = factory(env)
    #         fields.entries[name] = built
    #         built_fields[name] = built
    #     else:
    #         factory = models.entries[name]
    #         built = factory(env)
    #         models.entries[name] = built
    #         built_models[name] = built