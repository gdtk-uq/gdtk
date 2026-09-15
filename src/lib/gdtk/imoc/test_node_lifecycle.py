# test_node_lifecycle.py
"""
Tests for the node lifecycle functions of the IMOC kernel.

These check the behaviour against the original moc C code that they were
ported from, and against its Tcl procedure ExtendStreamLineToGivenLine.

Usage:
  pytest test_node_lifecycle.py
"""

import math
import pytest

import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit
from gdtk.ideal_gas_flow import PM1


@pytest.fixture(autouse=True)
def fresh_kernel():
    """Every test starts with an empty mesh."""
    kernel.delete_all_nodes()
    kernel.axisymmetric = False
    kernel.g = 1.4
    yield
    kernel.delete_all_nodes()


def make_line(n, attr_up='cminus_up', attr_down='cminus_down'):
    """
    Make a chain of n nodes linked along one characteristic family,
    running up from index 0.  Returns the list of indices.
    """
    ids = [kernel.create_node() for _ in range(n)]
    for a, b in zip(ids[:-1], ids[1:]):
        setattr(kernel.nodes[a], attr_down, b)
        setattr(kernel.nodes[b], attr_up, a)
    return ids


def walk(start, attr):
    """Follow a chain of links from start, returning the indices visited."""
    out = [start]
    while True:
        nxt = getattr(kernel.nodes[out[-1]], attr)
        if nxt is None: break
        out.append(nxt)
    return out


# ---------------------------------------------------------------- validity

def test_valid_node():
    i = kernel.create_node()
    assert kernel.valid_node(i)
    assert kernel.valid_node(kernel.nodes[i])
    kernel.delete_node(i)
    assert not kernel.valid_node(i)
    # Out of range and nonsense arguments are simply not valid.
    assert not kernel.valid_node(-1)
    assert not kernel.valid_node(len(kernel.nodes))
    assert not kernel.valid_node(None)
    assert not kernel.valid_node("3")


def test_number_of_nodes_is_the_live_count():
    """
    GetNumberOfNodes in the C kernel counts live nodes, not slots.
    Scripts ported from Tcl rely on the difference.
    """
    ids = [kernel.create_node() for _ in range(5)]
    assert kernel.number_of_nodes() == 5
    kernel.delete_node(ids[2])
    assert kernel.number_of_nodes() == 4
    assert len(kernel.nodes) == 5      # the slot is retained
    kernel.delete_node(ids[2])         # deleting twice changes nothing
    assert kernel.number_of_nodes() == 4


def test_get_next_node_id_skips_holes():
    ids = [kernel.create_node() for _ in range(4)]
    kernel.delete_node(ids[1])
    kernel.delete_node(ids[2])
    seen = []
    i = kernel.get_next_node_id(-1)
    while i >= 0:
        seen.append(i)
        i = kernel.get_next_node_id(i)
    assert seen == [ids[0], ids[3]]


# ------------------------------------------------------------------ create

def test_create_node_reuses_lowest_free_slot():
    """
    CreateNode(-1) in the C kernel takes the first available space in the
    pointer array.  Ported scripts can depend on the indices that come back.
    """
    ids = [kernel.create_node() for _ in range(5)]
    assert ids == [0, 1, 2, 3, 4]
    kernel.delete_node(3)
    kernel.delete_node(1)
    assert kernel.create_node() == 1   # lowest free slot, not the newest
    assert kernel.create_node() == 3
    assert kernel.create_node() == 5   # none free, so append
    assert len(kernel.nodes) == 6


def test_create_node_at_a_given_slot_pads_the_list():
    """
    The original moc C code had a preallocated pointer array, so
    CreateNode(101) worked on an empty mesh.  Its demo scripts use that,
    e.g. the one for Anderson's Exercise 11.1.
    """
    i = kernel.create_node(101)
    assert i == 101
    assert kernel.valid_node(101)
    assert kernel.nodes[101].indx == 101
    assert kernel.number_of_nodes() == 1      # the padding is not live
    assert len(kernel.nodes) == 102
    # The padding is free for reuse, lowest first.
    assert kernel.create_node() == 0


def test_create_node_at_an_occupied_slot_replaces_it():
    i = kernel.create_node()
    kernel.nodes[i].x = 5.0
    j = kernel.create_node(i)
    assert j == i
    assert kernel.nodes[i].x == 0.0           # a fresh, empty node
    assert kernel.number_of_nodes() == 1


def test_node_beyond_the_end_is_refused():
    """A Node stored at a position that does not match its indx corrupts links."""
    kernel.create_node()
    with pytest.raises(RuntimeError):
        kernel.Node(indx=7)


# ------------------------------------------------------------------ delete

def test_delete_node_splices_the_characteristic():
    """
    DeleteNode in the C kernel retains the linkages between the remaining
    nodes, so a line running through the deleted node stays walkable.
    """
    ids = make_line(5)
    kernel.delete_node(ids[2])
    assert walk(ids[0], 'cminus_down') == [ids[0], ids[1], ids[3], ids[4]]
    assert walk(ids[4], 'cminus_up') == [ids[4], ids[3], ids[1], ids[0]]


def test_delete_node_at_the_end_severs_cleanly():
    ids = make_line(3)
    kernel.delete_node(ids[2])
    assert walk(ids[0], 'cminus_down') == [ids[0], ids[1]]
    assert kernel.nodes[ids[1]].cminus_down is None


def test_delete_node_leaves_no_link_into_the_dead_slot():
    ids = make_line(4, 'cplus_up', 'cplus_down')
    victim = ids[1]
    kernel.delete_node(victim)
    for i in range(len(kernel.nodes)):
        if not kernel.valid_node(i): continue
        n = kernel.nodes[i]
        for attr in ('cplus_up', 'cplus_down', 'cminus_up', 'cminus_down',
                     'czero_up', 'czero_down'):
            assert getattr(n, attr) != victim, f"nodes[{i}].{attr} still points at it"


def test_delete_node_splices_all_three_families_at_once():
    up_m, down_m = kernel.create_node(), kernel.create_node()
    up_p, down_p = kernel.create_node(), kernel.create_node()
    up_z, down_z = kernel.create_node(), kernel.create_node()
    i = kernel.create_node()
    n = kernel.nodes[i]
    n.cminus_up, n.cminus_down = up_m, down_m
    n.cplus_up, n.cplus_down = up_p, down_p
    n.czero_up, n.czero_down = up_z, down_z
    kernel.nodes[up_m].cminus_down = i; kernel.nodes[down_m].cminus_up = i
    kernel.nodes[up_p].cplus_down = i;  kernel.nodes[down_p].cplus_up = i
    kernel.nodes[up_z].czero_down = i;  kernel.nodes[down_z].czero_up = i
    kernel.delete_node(i)
    assert kernel.nodes[up_m].cminus_down == down_m
    assert kernel.nodes[down_m].cminus_up == up_m
    assert kernel.nodes[up_p].cplus_down == down_p
    assert kernel.nodes[down_p].cplus_up == up_p
    assert kernel.nodes[up_z].czero_down == down_z
    assert kernel.nodes[down_z].czero_up == up_z


def test_delete_node_cleans_the_registries():
    ids = make_line(3, 'czero_up', 'czero_down')
    kernel.register_node_in_mesh(ids[1])
    kernel.register_streamline_start(ids[0])
    kernel.delete_node(ids[1])
    assert ids[1] not in kernel.char_mesh
    # Deleting the registered start promotes the next node along,
    # so the streamline is not lost.
    kernel.delete_node(ids[0])
    assert kernel.streamlines == [ids[2]]


def test_top_down_tail_deletion_keeps_survivor_indices():
    """
    The rollback idiom of the Tcl nozzle scripts: delete from the top down,
    back to a marker node.  The survivors must keep their indices.
    """
    ids = make_line(10)
    reset = ids[4]
    for i in range(len(kernel.nodes)-1, reset, -1):
        kernel.delete_node(i)
    assert kernel.number_of_nodes() == 5
    for i in range(reset+1):
        assert kernel.valid_node(i) and kernel.nodes[i].indx == i
    assert walk(ids[0], 'cminus_down') == ids[0:5]
    assert kernel.nodes[reset].cminus_down is None


# ------------------------------------------------------------------ search

def test_find_nodes_near_ignores_deleted_slots():
    a = kernel.create_node(); kernel.nodes[a].x, kernel.nodes[a].y = 0.0, 0.0
    b = kernel.create_node(); kernel.nodes[b].x, kernel.nodes[b].y = 0.01, 0.0
    c = kernel.create_node(); kernel.nodes[c].x, kernel.nodes[c].y = 1.0, 0.0
    kernel.delete_node(b)
    assert list(kernel.find_nodes_near(0.011, 0.0, tol=0.0)) == [a]
    near = list(kernel.find_nodes_near(0.0, 0.0, tol=0.5))
    assert b not in near and a in near


def test_find_nodes_near_with_kdtree_ignores_deleted_slots():
    for k in range(6):
        i = kernel.create_node()
        kernel.nodes[i].x, kernel.nodes[i].y = 0.1*k, 0.0
    kernel.delete_node(1)
    kernel.delete_node(2)
    kdt = kernel.create_kd_tree()
    near = list(kernel.find_nodes_near(0.15, 0.0, tol=0.5, max_count=10, kdtree=kdt))
    assert 1 not in near and 2 not in near
    assert 0 in near and 3 in near
    assert list(kernel.find_nodes_near(0.1, 0.0, tol=0.0, kdtree=kdt)) == [0]


def test_search_is_unchanged_when_nothing_is_deleted():
    for k in range(6):
        i = kernel.create_node()
        kernel.nodes[i].x, kernel.nodes[i].y = 0.1*k, 0.0
    kdt = kernel.create_kd_tree()
    # Nodes sit at x = 0.0, 0.1, ... 0.5; from x = 0.25 a radius of 0.16
    # takes in the four middle ones and leaves out the two at the ends.
    slow = list(kernel.find_nodes_near(0.25, 0.0, tol=0.16, max_count=10))
    fast = list(kernel.find_nodes_near(0.25, 0.0, tol=0.16, max_count=10, kdtree=kdt))
    assert sorted(slow) == sorted(fast) == [1, 2, 3, 4]
    assert kernel.find_nodes_near(0.31, 0.0, tol=0.0)[0] == 3
    assert kernel.find_nodes_near(0.31, 0.0, tol=0.0, kdtree=kdt)[0] == 3


# --------------------------------------------- extend_streamline_to_given_line

def uniform_flow_node(x, y, mach=2.0, theta=0.0):
    i = kernel.create_node()
    n = kernel.nodes[i]
    n.x = x; n.y = y; n.mach = mach; n.theta = theta; n.nu = PM1(mach, kernel.g)
    return i


def test_extend_streamline_picks_the_crossed_segment():
    """
    A streamline running along y = 0.5 towards +x should cross the middle
    segment of a vertical line of three nodes at x = 1.
    """
    start = uniform_flow_node(0.0, 0.5)
    line = [uniform_flow_node(1.0, y) for y in (0.0, 0.4, 0.9)]
    for a, b in zip(line[:-1], line[1:]):
        kernel.nodes[a].cplus_down = b
        kernel.nodes[b].cplus_up = a
    new = unit.extend_streamline_to_given_line(start, line)
    assert new is not None
    assert math.isclose(kernel.nodes[new].x, 1.0, abs_tol=1.0e-9)
    assert 0.4 <= kernel.nodes[new].y <= 0.9
    assert kernel.nodes[new].czero_up == start
    assert kernel.nodes[start].czero_down == new


def test_extend_streamline_returns_none_when_no_segment_is_crossed():
    """The Tcl procedure returns -1 here; we report None."""
    start = uniform_flow_node(0.0, 5.0)
    line = [uniform_flow_node(1.0, y) for y in (0.0, 0.4, 0.9)]
    assert unit.extend_streamline_to_given_line(start, line) is None


def test_extend_streamline_needs_at_least_two_nodes():
    start = uniform_flow_node(0.0, 0.5)
    assert unit.extend_streamline_to_given_line(start, []) is None
    assert unit.extend_streamline_to_given_line(start, [start]) is None


def test_extend_streamline_takes_the_first_valid_intersection():
    """
    The original Tcl procedure searches the segments in the order given and
    breaks on the first valid intersection, so the order of the list matters.
    """
    start = uniform_flow_node(0.0, 0.5)
    line = [uniform_flow_node(1.0, y) for y in (0.0, 0.4, 0.9, 2.0)]
    for a, b in zip(line[:-1], line[1:]):
        kernel.nodes[a].cplus_down = b
        kernel.nodes[b].cplus_up = a
    new = unit.extend_streamline_to_given_line(start, line)
    # The 0.4-to-0.9 segment comes first and contains y = 0.5.
    assert 0.4 <= kernel.nodes[new].y <= 0.9


# ------------------------------------------------------- delete in a real mesh

def test_delete_an_interior_node_of_a_computed_mesh():
    """
    Build a small mesh with the unit processes, delete a node from the middle
    of a C- line, and check the line is still walkable end to end.
    """
    kernel.axisymmetric = False
    left = []
    for k in range(4):
        i = kernel.create_node()
        n = kernel.nodes[i]
        n.x = 0.0; n.y = 0.2*k + 0.1
        n.nu = math.radians(10.0 + 2.0*k); n.theta = math.radians(2.0*k)
        n.mach = unit.igf.PM2(n.nu, kernel.g)
        left.append(i)
    for a, b in zip(left[:-1], left[1:]):
        kernel.nodes[a].cminus_down = b
        kernel.nodes[b].cminus_up = a
    start = kernel.create_node()
    s = kernel.nodes[start]
    s.x = 0.05; s.y = 0.05
    s.nu = math.radians(10.0); s.theta = 0.0
    s.mach = unit.igf.PM2(s.nu, kernel.g)
    new_line = unit.march_along_cminus(left[0], start, 'down')
    assert len(new_line) == len(left) + 1
    # interior() orients each link by the x-ordering of the two nodes, so the
    # new line may be chained either way round; take whichever it used.
    attr = ('cminus_down' if kernel.nodes[new_line[0]].cminus_down == new_line[1]
            else 'cminus_up')
    assert walk(new_line[0], attr) == new_line
    victim = new_line[2]
    kernel.delete_node(victim)
    walked = walk(new_line[0], attr)
    assert victim not in walked
    assert walked == [i for i in new_line if i != victim]
