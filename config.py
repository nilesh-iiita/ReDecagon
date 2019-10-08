import os
CUR_DIR = os.path.dirname(os.path.realpath(__file__))

DECAGON_DIR = "/home/anhnd/DTI Project/Data/DECAGON"

POLY_ADR_PATH = "%s/bio-decagon-combo.csv" % DECAGON_DIR
PRO_PRO_PATH = "%s/bio-decagon-ppi.csv" % DECAGON_DIR
DRUG_PRO_PATH = "%s/bio-decagon-targets.csv" % DECAGON_DIR

NUM_SE = 300

PROCESSED_DATA_DIR = "%s/data" % CUR_DIR

PROCESSED_COMBO_ADR = "%s/combo.dat" % PROCESSED_DATA_DIR
PROCESSED_DRUGPROTEIN = "%s/drugProtein.dat" % PROCESSED_DATA_DIR
PROCESSED_PROTEINDRUG = "%s/proteinDrug.dat" % PROCESSED_DATA_DIR
PROCESSED_PROTEINPROTEIN = "%s/proteinProtein.dat" % PROCESSED_DATA_DIR
DRUG_MAP = "%s/drug2Id.dat" % PROCESSED_DATA_DIR
PROTEIN_MAP = "%s/protein2Id.dat" % PROCESSED_DATA_DIR

