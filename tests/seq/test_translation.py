import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from cb.seq.translation import check_residue, translate_cds

# First 8 HBB codons (M V H L T P E E) plus a stop codon.
CDS = "ATGGTGCACCTGACTCCTGAGGAGTAA"
UTR = "GCC"  # 5' UTR, so the CDS does not start at position 0


def make_record(cds: str, with_cds_feature: bool = True) -> SeqRecord:
    seq = Seq(UTR + cds)
    features = [SeqFeature(FeatureLocation(len(UTR), len(seq)), type="CDS")] if with_cds_feature else []
    return SeqRecord(seq, id="TEST", features=features)


def test_translate_cds_stops_at_stop_codon():
    assert str(translate_cds(make_record(CDS))) == "MVHLTPEE"


def test_translate_cds_ignores_utr():
    assert translate_cds(make_record(CDS))[0] == "M"


def test_translate_cds_without_cds_feature_raises():
    with pytest.raises(ValueError, match="no CDS feature"):
        translate_cds(make_record(CDS, with_cds_feature=False))


def test_reference_residue_6_is_glu():
    protein = translate_cds(make_record(CDS))
    assert protein[6] == "E"
    assert check_residue(protein, 6, "E", "V") == "reference"


def test_glu6val_is_detected():
    mutated = CDS[:18] + "GTG" + CDS[21:]  # codon index 6: GAG -> GTG
    protein = translate_cds(make_record(mutated))
    assert protein[6] == "V"
    assert check_residue(protein, 6, "E", "V") == "variant"


@pytest.mark.parametrize("protein", ["MVHLTPEE", Seq("MVHLTPEE")])
def test_check_residue_accepts_str_and_seq(protein):
    assert check_residue(protein, 6, "E", "V") == "reference"


def test_wrong_offset_fails_loudly():
    protein = translate_cds(make_record(CDS))
    with pytest.raises(AssertionError, match="index 5"):
        check_residue(protein, 5, "E", "V")


def test_index_out_of_range_raises():
    with pytest.raises(IndexError):
        check_residue("MV", 6, "E", "V")
