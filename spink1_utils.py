import numpy as np
import pandas as pd
import numpy.linalg as la
from scipy import stats
import gseapy as gp

MYELOID = ['CD68', 'AIF1', 'ITGAM', 'CSF1R', 'TYROBP', 'FCER1G',
           'C1QA', 'C1QB', 'C1QC', 'MRC1', 'CD163']

CELL_PANELS = {
    'Myeloid / macrophage': MYELOID,
    'Neutrophil': ['FCGR3B', 'CXCR2', 'CSF3R', 'SELL', 'S100A12', 'MNDA', 'FPR1'],
    'T / NK cell': ['CD2', 'CD3D', 'CD3E', 'CD247', 'IL2RB', 'GZMA', 'PRF1'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'CLDN5', 'ESAM'],
    'Stromal / ECM': ['COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'DCN', 'LUM'],
    'Astrocyte': ['AQP4', 'GJA1', 'SLC1A2', 'ALDH1L1'],
    'Neuron': ['SYT1', 'SNAP25', 'RBFOX3', 'NEFL', 'SYN1'],
    'Oligodendrocyte': ['MBP', 'PLP1', 'MOG', 'MAG', 'CNP'],
    'OPC / proneural': ['PDGFRA', 'CSPG4', 'OLIG1', 'OLIG2', 'SOX10'],
    'Proliferation': ['MKI67', 'TOP2A', 'CCNB1', 'BUB1', 'CDK1', 'PCNA'],
}

SUBTYPE_PANELS = {
    'Proneural': ['PDGFRA', 'OLIG2', 'SOX2', 'DLL3', 'ASCL1', 'NKX2-2', 'ERBB3',
                  'CDKN1C', 'SOX10', 'MYT1', 'GRIA2', 'ATCAY', 'DCX', 'CSPG4'],
    'Classical': ['EGFR', 'NES', 'NOTCH3', 'JAG1', 'LFNG', 'SMO', 'GAS1', 'GLI2',
                  'AKT2', 'FGFR3', 'SOX9', 'MEOX2', 'PTPRA', 'ELOVL2'],
    'Mesenchymal': ['CHI3L1', 'MET', 'CD44', 'SERPINE1', 'TIMP1', 'LOX', 'PLAUR',
                    'CAST', 'COL5A1', 'ITGB1', 'TRADD', 'RELB', 'TNFRSF1A', 'TGFB1'],
}

MES_CLEAN = ['CHI3L1', 'MET', 'CD44', 'TIMP1', 'CAST', 'ITGB1',
             'TRADD', 'RELB', 'TNFRSF1A', 'TGFB1']

GROUP_ORDER = ['Undetectable', 'Low', 'High']


def hallmark_sets():
    return gp.get_library(name='MSigDB_Hallmark_2020')


def proliferation_genes(lib):
    return sorted(set(lib['E2F Targets']) | set(lib['G2-M Checkpoint']))


def panel_score(L, genes, drop=None):
    g = [x for x in genes if x in L.index and x != drop]
    Z = L.loc[g].T
    return ((Z - Z.mean()) / Z.std(ddof=0)).mean(axis=1)


def partial_spearman(x, y, Z):
    rank = lambda a: pd.Series(np.asarray(a, dtype=float)).rank().values
    x, y = rank(x), rank(y)
    M = np.column_stack([np.ones(len(x))] + [rank(Z[c]) for c in Z.columns])
    rx = x - M @ la.lstsq(M, x, rcond=None)[0]
    ry = y - M @ la.lstsq(M, y, rcond=None)[0]
    r = np.corrcoef(rx, ry)[0, 1]
    dof = len(x) - 2 - Z.shape[1]
    return r, 2 * stats.t.sf(abs(r * np.sqrt(dof / (1 - r ** 2))), dof)


def spink1_groups(s):
    cut = s[s > 0].median()
    return pd.Series(np.where(s == 0, 'Undetectable',
                              np.where(s <= cut, 'Low', 'High')), index=s.index)
