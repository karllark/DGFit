
from ..DGFit_Models import DGFit_MRN


def test_mrn_initialize():
    dgmod = DGFit_MRN()
    assert dgmod.type == 'MRN'
