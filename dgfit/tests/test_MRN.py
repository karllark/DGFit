from dgfit.dustmodel import MRNDustModel


def test_mrn_initialize():
    dmod = MRNDustModel()
    assert dmod.sizedisttype == "MRN"
