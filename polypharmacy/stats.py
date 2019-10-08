import config
from collections import defaultdict
import networkx as nx


def load_combo_se(fname=config.POLY_ADR_PATH):
    combo2stitch = {}
    combo2se = defaultdict(set)
    se2name = {}
    drugs = set()
    fin = open(fname)
    print 'Reading: %s' % fname
    fin.readline()
    for line in fin:
        stitch_id1, stitch_id2, se, se_name = line.strip().split(',')
        drugs.add(stitch_id1)
        drugs.add(stitch_id2)
        combo = stitch_id1 + '_' + stitch_id2
        combo2stitch[combo] = [stitch_id1, stitch_id2]
        combo2se[combo].add(se)
        se2name[se] = se_name
    fin.close()
    n_interactions = sum([len(v) for v in combo2se.values()])
    print 'Drug combinations: %d Side effects: %d' % (len(combo2stitch), len(se2name))
    print 'Drug-drug interactions: %d' % (n_interactions)
    print 'Num drug: %d' % (len(drugs))
    return combo2stitch, combo2se, se2name


def load_ppi(fname=config.PRO_PRO_PATH):
    fin = open(fname)
    print 'Reading: %s' % fname
    fin.readline()
    edges = []
    for line in fin:
        gene_id1, gene_id2= line.strip().split(',')
        edges += [[gene_id1,gene_id2]]
    nodes = set([u for e in edges for u in e])
    print 'Edges: %d' % len(edges)
    print 'Nodes: %d' % len(nodes)
    net = nx.Graph()
    net.add_edges_from(edges)
    net.remove_nodes_from(nx.isolates(net))
    net.remove_edges_from(net.selfloop_edges())
    node2idx = {node: i for i, node in enumerate(net.nodes())}
    return net, node2idx


def load_targets(fname = config.DRUG_PRO_PATH):
    stitch2proteins = defaultdict(set)
    fin = open(fname)
    print 'Reading: %s' % fname
    fin.readline()
    gens = set()
    for line in fin:
        stitch_id, gene = line.strip().split(',')
        stitch2proteins[stitch_id].add(gene)
        gens.add(gene)
    print len(gens)



    fin = open(config.PRO_PRO_PATH)
    print 'Reading: %s' % fname
    fin.readline()
    edges = []
    for line in fin:
        gene_id1, gene_id2= line.strip().split(',')
        edges += [[gene_id1,gene_id2]]
    nodes = set([u for e in edges for u in e])
    print 'Edges: %d' % len(edges)
    print 'Nodes: %d' % len(nodes)

    cc = 0
    for node in nodes:
        if node not in gens:
            cc += 1
    print cc, len(nodes), len(gens)

    return stitch2proteins


if __name__ == "__main__":
    load_targets()