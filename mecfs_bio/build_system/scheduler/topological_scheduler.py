from copy import deepcopy
from pathlib import Path
from typing import Mapping, Sequence

import emoji
import networkx as nx
import structlog
from loguru import logger

logger = structlog.get_logger()

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.base_rebuilder import Rebuilder
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)
from mecfs_bio.build_system.scheduler.asset_retrieval import get_asset_if_exists
from mecfs_bio.build_system.tasks.base_tasks import Tasks
from mecfs_bio.build_system.wf.base_wf import WF


def get_dependency_graph_from_tasks(tasks: Tasks) -> nx.DiGraph:
    G: nx.DiGraph = nx.DiGraph()
    for asset_id, task in tasks.items():
        G.add_node(asset_id)
        for dep in task.deps:
            G.add_edge(dep.asset_id, asset_id)
    return G


def dependencies_of_targets_from_tasks(
    tasks: Tasks, targets: Sequence[AssetId]
) -> nx.DiGraph:
    """
    Given a set of tasks, and a list of targets, build a minimal dependency graph containing the targets
    and all their transitive dependencies.
    """
    G = get_dependency_graph_from_tasks(tasks)
    return dependencies_of_targets(G, targets)


def dependees_of_targets_from_tasks(
    tasks: Tasks, targets: Sequence[AssetId]
) -> nx.DiGraph:
    """
    Given a set of tasks, and a list of targets, build a minimal dependency graph containing the targets
    and all their transitive dependees.
    """
    G = get_dependency_graph_from_tasks(tasks)
    return dependees_of_targets(G, targets)


def dependencies_of_targets(G: nx.DiGraph, targets: Sequence[AssetId]) -> nx.DiGraph:
    """
    Given a dependency graph of assets, find the transitive dependencies of given targets
    """
    reachable = set(targets)
    for target in targets:
        reachable = nx.ancestors(G, target) | reachable
    subgraph = nx.DiGraph(G.subgraph(reachable))
    return subgraph


def dependees_of_targets(G: nx.DiGraph, targets: Sequence[AssetId]) -> nx.DiGraph:
    """
    Given a dependency graph of assets, find the transitive dependees of given targets
    """
    reachable = set(targets)
    for target in targets:
        reachable = nx.descendants(G, target) | reachable
    subgraph = nx.DiGraph(G.subgraph(reachable))
    return subgraph


def _get_initial_frontier(g: nx.DiGraph) -> list[AssetId]:
    """
    Frontier is defined as the set of assets without dependencies, or whose dependencies have been built
    """
    generations = nx.topological_generations(g)
    list_of_generations = [sorted(gen) for gen in generations]
    return list_of_generations[0]


def _update_frontier_and_G(
    G: nx.DiGraph, completed: AssetId, frontier: list[AssetId]
) -> tuple[list[AssetId], nx.DiGraph]:
    """
    After competing an asset, any nodes whose only predecessor was that asset can be added to the frontier
    """
    frontier = deepcopy(frontier)
    G = deepcopy(G)
    succs = list(G.successors(completed))
    G.remove_node(completed)
    frontier.remove(completed)
    for node in succs:
        if len(list(G.predecessors(node))) == 0:
            frontier.append(node)
    return frontier, G


def _get_progress_list(todo: Sequence[AssetId], done: set[AssetId]) -> str:
    s = "\nWork Progress:\n"
    for asset_id in todo:
        if asset_id in done:
            s += f"{emoji.emojize(':check_mark_button:')} {asset_id}\n"
        else:
            s += f"{asset_id}\n"
    return s


def topological[
    Info,
](
    rebuilder: Rebuilder[Info],
    tasks: Tasks,
    targets: Sequence[AssetId],
    wf: WF,
    info: Info,
    meta_to_path: MetaToPath,
    incremental_save_path: Path | None = None,
) -> tuple[Mapping[AssetId, Asset], Info]:
    """
    A scheduler that builds a dependency graph of tasks, and executes them in topological order.
    Based on
    Mokhov, Andrey, Neil Mitchell, and Simon Peyton Jones.
    "Build systems à la carte: Theory and practice." Journal of Functional Programming 30 (2020): e11.
    """
    G = dependencies_of_targets_from_tasks(tasks, targets)
    todo: list[AssetId] = list(nx.topological_sort(G))
    done: set[AssetId] = set()
    store: dict[AssetId, Asset] = {}
    frontier = _get_initial_frontier(G)
    while len(G) > 0:
        logger.info(_get_progress_list(todo=todo, done=done))
        node = frontier[0]
        task = tasks[node]
        maybe_asset = get_asset_if_exists(
            meta=task.meta,
            meta_to_path=meta_to_path,
        )

        class SimpleFetch(Fetch):
            def __call__(self, asset_id: AssetId) -> Asset:
                return store[asset_id]

        new_asset, info = rebuilder.rebuild(
            task=task,
            asset=maybe_asset,
            fetch=SimpleFetch(),
            wf=wf,
            info=info,
            meta_to_path=meta_to_path,
        )
        if incremental_save_path is not None:
            rebuilder.save_info(info, incremental_save_path)
        store[node] = new_asset
        frontier, G = _update_frontier_and_G(G=G, completed=node, frontier=frontier)
        done.add(node)
    logger.info(_get_progress_list(todo=todo, done=done))
    return store, info
