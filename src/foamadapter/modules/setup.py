
from graphlib import TopologicalSorter, CycleError
import networkx as nx
import matplotlib.pyplot as plt
from typing import Dict, Set, Iterable, List, Tuple, Any, Optional
from pathlib import Path
from foamadapter.modules.fields import Field
from foamadapter.modules.models import Models, Model
from typing import Protocol, Dict, Set, Any

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
) -> Tuple[Dict[Any, Tuple[float, float]], Dict[Any, int], List[Set[str]]]:
    """
    Compute layered positions using networkx.multipartite_layout.
    - Assigns node attribute 'subset' = topological generation index.
    - Uses align='vertical' so subsets are stacked along Y (levels).
    - Ensures generation 0 is at the TOP (top_to_bottom=True).
    Returns: (pos, subset_by_node, layers)
    """
    if not nx.is_directed_acyclic_graph(G):
        raise nx.NetworkXUnfeasible("Graph is not a DAG (cycle detected).")

    layers = list(nx.topological_generations(G))  # [set(nodes), ...]
    subset_by_node = {node: i for i, layer in enumerate(layers) for node in layer}
    nx.set_node_attributes(G, subset_by_node, name="subset")

    # Stack along Y (vertical): nodes of the same subset share X? No—here they share Y (different Y per subset).
    pos = nx.multipartite_layout(G, subset_key="subset", align="horizontal")

    if layers:
        # Make sure layer 0 is visually at the TOP (largest y). If not, flip Y.
        layer_y = {i: sum(pos[n][1] for n in layer) / max(1, len(layer)) for i, layer in enumerate(layers)}
        if layer_y:  # Only process if there are layers
            top_layer_index_by_y = max(layer_y, key=lambda k: layer_y[k])  # layer with highest y
            if top_layer_index_by_y != 0:
                for n, (x, y) in pos.items():
                    pos[n] = (x, -y)

    return pos, subset_by_node, layers


# ----------------------------
# Color palette generation
# ----------------------------
def generate_color_palette(n_colors: int = 8) -> List[str]:
    """
    Generate a visually distinct color palette with light, pastel colors.
    Uses HSV color space to ensure good visual separation.
    """
    import colorsys
    
    colors = []
    for i in range(n_colors):
        # Distribute hues evenly around the color wheel
        hue = i / n_colors
        # Use high saturation and value for vibrant but light colors
        saturation = 0.3  # Low saturation for pastel effect
        value = 0.95      # High value for light colors
        
        # Convert HSV to RGB
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        # Convert to hex
        hex_color = f"#{int(rgb[0]*255):02x}{int(rgb[1]*255):02x}{int(rgb[2]*255):02x}"
        colors.append(hex_color)
    
    return colors


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
    if not pos:  # Handle empty graphs
        return
        
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
    level_lines: bool = True,
) -> nx.DiGraph:
    """
    Validate with graphlib, build a networkx graph, compute layered layout,
    and draw the DAG. Saves to filename if provided.
    
    Automatically detects container-based naming (e.g., "fields:temperature", "models:viscosity")
    and applies container-specific styling when present.
    
    - level_lines: draw dotted guides per level.
    """
    try:
        order = validate_with_graphlib(deps)
        print("✅ Valid DAG. Topological order:", order)
    except CycleError as e:
        print("❌ Cycle detected:", e.args)
        raise

    G = build_nx_graph(deps)
    pos, subset_by_node, layers = layered_positions(G)

    # Create figure/axes
    fig, ax = plt.subplots(figsize=(12, 10))

    # Optional level guides
    if level_lines:
        draw_level_guides(ax, pos, subset_by_node)

    all_nodes = list(G.nodes)
    
    # Auto-detect container-based naming (CONTAINER:NAME format)
    container_nodes = {}
    simple_nodes = []
    
    for node in all_nodes:
        if ":" in node:
            container_type = node.split(":", 1)[0]
            if container_type not in container_nodes:
                container_nodes[container_type] = []
            container_nodes[container_type].append(node)
        else:
            simple_nodes.append(node)
    
    # If container-based nodes detected, use container-specific styling
    if container_nodes:
        # Auto-generate color palette based on number of container types
        n_containers = len(container_nodes)
        color_palette = generate_color_palette(max(n_containers, 8))
        
        # Assign colors to container types based on order of appearance
        container_types = list(container_nodes.keys())
        container_colors = {
            container_type: color_palette[i % len(color_palette)]
            for i, container_type in enumerate(container_types)
        }
        
        # Draw nodes by container type
        for container_type, nodes in container_nodes.items():
            color = container_colors[container_type]
            nx.draw_networkx_nodes(
                G, pos, nodelist=nodes,
                node_size=3800, node_color=color, 
                edgecolors="black", node_shape="o", ax=ax
            )
        
        # Draw any non-container nodes with default styling
        if simple_nodes:
            nx.draw_networkx_nodes(
                G, pos, nodelist=simple_nodes,
                node_size=3800, node_color="#EEEEEE", edgecolors="black", ax=ax
            )
    else:
        # No container naming detected, use uniform styling
        nx.draw_networkx_nodes(
            G, pos, 
            node_size=3800, node_color="#EEEEEE", 
            edgecolors="black", ax=ax
        )

    # Labels for all nodes
    nx.draw_networkx_labels(G, pos, font_size=9, ax=ax)

    # Draw all edges uniformly
    nx.draw_networkx_edges(
        G, pos, arrows=True, arrowstyle='-|>', arrowsize=18, width=1.2,
        edge_color="black", min_source_margin=0, min_target_margin=24, ax=ax
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


def containers_to_deps(*containers: 'DependencyContainer') -> Dict[str, Set[str]]:
    """
    Convert multiple dependency containers into a single dependency dictionary
    with namespaced keys (CONTAINER:NAME format).
    """
    # Helper: namespaced key  
    def key(container, name):
        return f"{container.__class__.__name__.lower()}:{name}"

    # Merge all dependencies
    merged_deps: Dict[str, Set[str]] = {}
    
    for container in containers:
        container_deps = container.dependencies()
        for name, deps in container_deps.items():
            namespaced_name = key(container, name)
            # Convert dependency names to namespaced keys by searching all containers
            namespaced_deps = set()
            for dep in deps:
                for c in containers:
                    if dep in c.entries:
                        namespaced_deps.add(key(c, dep))
                        break
                else:
                    # If not found in any container, keep as-is (might be external)
                    namespaced_deps.add(dep)
            merged_deps[namespaced_name] = namespaced_deps
    
    return merged_deps


# --- DependencyContainer Protocol ---
class DependencyContainer(Protocol):
    entries: Dict[str, Any]
    def dependencies(self) -> Dict[str, Set[str]]:
        ...

# --- Generic Initialization Function ---
def initialize_containers(*containers: DependencyContainer) -> None:
    """
    Initialize all dependency containers (Fields, Models, etc.) by resolving factories in topological order.
    """
    from graphlib import TopologicalSorter

    # Helper: namespaced key
    def key(container, name):
        return f"{container.__class__.__name__.lower()}:{name}"

    # Helper: find a dependency across all containers
    def find_dependency(dep_name: str, containers_list) -> str:
        for container in containers_list:
            if dep_name in container.entries:
                return key(container, dep_name)
        raise KeyError(f"Dependency '{dep_name}' not found in any container")

    dag: Dict[str, Set[str]] = {}
    is_factory: Dict[str, bool] = {}
    entry_map: Dict[str, Any] = {}
    container_map: Dict[str, Any] = {}

    # Collect all entries and dependencies
    for container in containers:
        for name, entry in container.entries.items():
            k = key(container, name)
            deps = getattr(entry, "dependencies", [])
            # Resolve dependencies across all containers
            dag[k] = {find_dependency(d, containers) for d in deps}
            is_factory[k] = callable(entry) and hasattr(entry, "dependencies")
            entry_map[k] = entry
            container_map[k] = container

    # Topological sort
    order = list(TopologicalSorter(dag).static_order())

    # Build in order
    built: Dict[str, Any] = {}
    for ref in order:
        entry = entry_map[ref]
        container = container_map[ref]
        if not is_factory[ref]:
            built[ref] = entry
            continue
        # Build with resolved dependencies from built
        dep_names = getattr(entry, "dependencies", [])
        env = {dep: built[find_dependency(dep, containers)] for dep in dep_names}
        built_entry = entry(env) if callable(entry) else entry
        container.entries[ref.split(":", 1)[1]] = built_entry
        built[ref] = built_entry
