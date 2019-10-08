import config, utils
from collections import defaultdict
import networkx as nx
import scipy.sparse as sp
import numpy as np


def filter_combo_se(fname=config.POLY_ADR_PATH):
    combo2stitch = {}
    combo2se = defaultdict(set)
    se2name = {}
    seCounter = dict()
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
        utils.add_dict_counter(seCounter, se)
        se2name[se] = se_name
    fin.close()
    n_interactions = sum([len(v) for v in combo2se.values()])

    print 'Before'
    print 'Drug combinations: %d Side effects: %d' % (len(combo2stitch), len(se2name))
    print 'Drug-drug interactions: %d' % (n_interactions)
    print 'Num drug: %d' % (len(drugs))

    seCounterSorted = utils.sort_dict(seCounter)
    validSe = set()
    se2Id = dict()
    for i in range(config.NUM_SE):
        se = seCounterSorted[i][0]
        validSe.add(seCounterSorted[i][0])
        utils.get_update_dict_index(se2Id, se)
    print validSe
    drug2Id = dict()
    se2Combos = dict()
    nInteractions = 0
    combosCounter = dict()
    for combo, ses in combo2se.iteritems():
        t1, t2 = combo2stitch[combo]
        id1 = utils.get_update_dict_index(drug2Id, t1)
        id2 = utils.get_update_dict_index(drug2Id, t2)
        for se in ses:
            seId = utils.get_dict(se2Id, se, -1)
            if seId != -1:
                combos = utils.get_insert_key_dict(se2Combos, seId, [])
                combos.append([id1, id2])
                nInteractions += 1
                utils.add_dict_counter(combosCounter, "%s_%s" % (id1, id2))

    nDrug = len(drug2Id)

    print 'After'
    print 'Drug combinations: %d Side effects: %d' % (len(combosCounter), config.NUM_SE)
    print 'Drug-drug interactions: %d' % nInteractions
    print 'Num drug: %d' % (nDrug)

    print 'Save to file...'

    drug_drug_adj_list = []
    for se, combos in se2Combos.iteritems():
        drugdrugMatrix = np.zeros((nDrug, nDrug))
        for d1, d2 in combos:
            drugdrugMatrix[d1, d2] = drugdrugMatrix[d2, d1] = 1
        drug_drug_adj_list.append(sp.csr_matrix(drugdrugMatrix))

    utils.save_obj(drug_drug_adj_list, config.PROCESSED_COMBO_ADR)
    utils.save_obj(drug2Id, config.DRUG_MAP)


def filter_targets(fname=config.DRUG_PRO_PATH):
    gen2Id = dict()
    drug2Id = utils.load_obj(config.DRUG_MAP)
    fin = open(fname)
    print 'Reading: %s' % fname
    fin.readline()
    drug2Protein = [[] for i in range(len(drug2Id))]

    for line in fin:
        stitch_id, gene = line.strip().split(',')
        genId = utils.get_update_dict_index(gen2Id, gene)
        drugId = utils.get_dict(drug2Id, stitch_id)
        proteins = drug2Protein[drugId]
        proteins.append(genId)

    print len(gen2Id)
    print gen2Id

    mat = np.zeros((len(drug2Id), len(gen2Id)))
    for drugId, proteins in enumerate(drug2Protein):
        for proteinId in proteins:
            mat[drugId][proteinId] = 1

    drugProteinAdj = sp.csr_matrix(mat)
    proteinDrugAdj = drugProteinAdj.transpose(copy=True)

    utils.save_obj(drugProteinAdj, config.PROCESSED_DRUGPROTEIN)
    utils.save_obj(proteinDrugAdj, config.PROCESSED_PROTEINDRUG)
    utils.save_obj(gen2Id, config.PROTEIN_MAP)

    # fin = open(config.PRO_PRO_PATH)
    # print 'Reading: %s' % fname
    # fin.readline()
    # edges = []
    # for line in fin:
    #     gene_id1, gene_id2 = line.strip().split(',')
    #     edges += [[gene_id1, gene_id2]]
    # nodes = set([u for e in edges for u in e])
    # print 'Edges: %d' % len(edges)
    # print 'Nodes: %d' % len(nodes)
    #
    # cc = 0
    # for node in nodes:
    #     if node not in gens:
    #         cc += 1
    # print cc, len(nodes), len(gens)
    #
    # return stitch2proteins


def filterPPI(fname=config.PRO_PRO_PATH):
    fin = open(fname)
    print 'Reading: %s' % fname

    gen2Id = utils.load_obj(config.PROTEIN_MAP)
    fin.readline()
    nGen = len(gen2Id)
    mat = np.zeros((nGen, nGen))
    for line in fin:
        gene_id1, gene_id2 = line.strip().split(',')
        id1 = utils.get_dict(gen2Id, gene_id1, -1)
        id2 = utils.get_dict(gen2Id, gene_id2, -1)
        if id1 != -1 and id2 != -1:
            mat[id1, id2] = mat[id2, id1] = 1
    mat = sp.csr_matrix(mat)
    utils.save_obj(mat, config.PROCESSED_PROTEINPROTEIN)

if __name__ == "__main__":
    # filter_combo_se()
    # filter_targets()
    filterPPI()
