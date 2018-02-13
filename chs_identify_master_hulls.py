"""
Module providing the code to identify potential master hulls for a
stack, using a simple overlap test.

"""

import six

import numpy as np

import chs_utils as utils


def make_graph():
    """Create an empty undirected graph.

    Returns
    -------
    gr : dict
        Don't peek under the covers!
    """

    return {'edges': {}}


def add_edge(gr, l1, l2):
    """Add an edge between two identifiers.

    Parameters
    ----------
    gr : dict
        The graph (e.g. as returned by make_graph)
    l1, l2 : identifiers
        They must be hashable.
    """

    assert 'edges' in gr

    # Rather than try to keep track of nodes here, just add
    # in the links and leave the hard work to get_nodes.
    #
    def add(x, y):
        try:
            gr['edges'][x].add(y)
        except KeyError:
            gr['edges'][x] = set([y])

    add(l1, l2)
    add(l2, l1)


def get_labels(gr):
    """Return the labels in the graph.

    Parameters
    ----------
    gr : dict
        The graph

    Returns
    -------
    labels : set of labels
        All the labels in the graph.
    """

    return set(gr['edges'].keys())


def get_nodes(gr):
    """Return the nodes in the graph.

    Parameters
    ----------
    gr : dict
        The graph

    Returns
    -------
    nodes : list of set
        Each element is the set of identifiers identified with a node.
        There is no guarantee to the ordering of this list.
    """

    # Here's where the simplicity of add_edge doesn't pay off...
    labels = gr['edges'].keys()
    store = {}
    nextNode = 1
    seen = set([])
    for label in labels:

        # Do we already have a node for this set of edges?
        #
        related = [label] + list(gr['edges'][label])
        nodes = set([])
        for l in related:
            n = -1
            try:
                n = store[l]
            except KeyError:
                continue

            nodes.add(n)

        nmatch = len(nodes)
        if nmatch == 0:
            node = nextNode
            nextNode += 1
        elif nmatch == 1:
            node = nodes.pop()
        else:
            # Have to merge two nodes.
            #
            # For now create a new node and retire the old ones
            # (by overwriting their labels).
            #
            # What are the labels associated with the old nodes?
            toadd = set([])
            for k, v in six.iteritems(store):
                if v not in nodes:
                    continue
                toadd.add(k)

            node = nextNode
            nextNode += 1

            # do not forgot to add in the extra nodes!
            related.extend(list(toadd))

        for l in related:
            store[l] = node
            seen.add(l)

    assert seen == set(labels)

    nodes = {}
    for label, node in six.iteritems(store):
        # do not need to use a set here, since by construction there's
        # only one label -> node possible, but we want sets in the
        # output
        try:
            nodes[node].add(label)
        except KeyError:
            nodes[node] = set([label])

    return list(nodes.values())


def show_graph(gr):
    """Print the graph for debugging.

    Parameters
    ----------
    gr : dict
        The graph
    """

    nodes = get_nodes(gr)
    nnodes = len(nodes)
    if nnodes == 0:
        print("The graph is empty.")
        return
    elif nnodes == 1:
        print("There is one node.")
    else:
        print("There are {} nodes.".format(nnodes))

    for i, node in enumerate(nodes):
        print("Node {}".format(i + 1))
        labels = " ".join("{}".format(l) for l in node)
        print("  {}".format(labels))

    print("")


def find_overlap_graph(hulls):
    """Return the overlapping hulls for this ensemble as a graph

    Parameters
    ----------
    hulls : list of dict
        The hulls for the stacks in an ensemble. Each entry is a
        separate hull, with the fields:
        stack, component, transform, pos, eqpos, and infile.
        It is assumed that the pos and eqpos arrays only contain
        finite values, and form closed regions (not yet sure if
        the last constraint is necessary).

    Returns
    -------
    gr, singles : dict, list of (cohort, component)
        The graph of the overlapping hulls and the hulls that do
        not overlap anything.

    """

    if hulls is None or hulls == []:
        raise ValueError("hulls can not be unset/empty")

    # Use a "simple" graph to represent the overlapping
    # hulls. This is *very* limited and - as it only has to
    # deal with small numbers - is not designed with efficiency
    # in mind.
    #
    # Note: ideally the hulls from the same stack would not
    # overlap, so we do not have to correlate hulls from the same
    # cohort. This is unfortunately not true, so this complicates
    # the following. Actually, we could probably enforce this
    # restriction now, but leaving in for now.
    #
    gr = make_graph()

    # What are the labels (i.e. cohort, component pairs)
    # for all the hulls. This lets the code identify those hulls
    # which have no overlaps.
    #
    def get_hull_key(h):
        return h['stack'], h['component']

    hullkeys = set([get_hull_key(h)
                    for h in hulls])

    # Considered as a matrix - comparing (stack,cpt)_i vs
    # (stack_cpt)_j - the overlap check should be symmetric (modulo
    # floating point errors and possible issues using a different
    # coordinate system), so we only need to loop through half
    # the matrix (excluding i=j).
    #
    nhulls = len(hulls)
    for i in range(0, nhulls - 1):

        # project the remaining hulls to this coordinate system;
        # this does lead to repeated work when there are multiple
        # hulls from the same stack, but it should not be a
        # significant time cost.
        #
        hulli = hulls[i]
        base_transform = hulli['transform']
        base_stack = hulli['stack']
        base_cpt = hulli['component']

        base_pos = hulli['pos']
        base_reg = utils.make_region_string(base_pos)
        base_key = (base_stack, base_cpt)
        hullkeys.add(base_key)

        for j in range(i + 1, nhulls):

            hullj = hulls[j]
            other_stack = hullj['stack']
            other_cpt = hullj['component']
            other_key = (other_stack, other_cpt)
            hullkeys.add(other_key)

            if other_stack == base_stack:
                other_pos = hullj['pos']

            else:
                # switch from 2D to 3D since this makes the transform
                # much easier to work with.
                #
                other_eqpos = hullj['eqpos'][np.newaxis]
                other_pos = base_transform.invert(other_eqpos)
                other_pos = other_pos[0]

                assert other_pos.shape == hullj['pos'].shape, \
                    'shapes differ: coding error!'

            other_reg = utils.make_region_string(other_pos)

            if utils.check_overlap(base_reg, other_reg):
                add_edge(gr, base_key, other_key)

    # Find the hulls which have no overlap; i.e. the ones not in
    # the overlap graph?
    #
    missing = sorted(hullkeys.difference(get_labels(gr)))

    return gr, missing
