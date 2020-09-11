from py4pm.pmfutilities import PMF

DATA_DIR="./tests/data/"

def test_instantiate_PMF():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)

def test_read_base_contributions():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_base_contributions()

def test_read_base_profiles():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_base_profiles()

def test_read_constrained_contributions():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_constrained_contributions()

def test_read_constrained_profiles():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_constrained_profiles()

def test_read_base_bootstrap():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_base_bootstrap()

def test_read_constrained_bootstrap():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_constrained_bootstrap()

def test_read_base_uncertainties_summary():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_base_uncertainties_summary()

def test_read_constrained_uncertainties_summary():
    pmf = PMF(site="GRE-cb", BDIR=DATA_DIR)
    pmf.read.read_constrained_uncertainties_summary()
