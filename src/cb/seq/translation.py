"""CDS extraction, translation and single-residue variant checks.

Numbering: indices are 0-based on the translated protein, which still has its
initiator Met at index 0. Legacy hemoglobin numbering (Glu6Val) skips the Met,
so residue 6 is index 6. HGVS counts the Met, so the same variant is p.Glu7Val.

Example:
    >>> protein = translate_cds(record)  # GenBank record, e.g. NM_000518
    >>> check_residue(protein, 6, ref="E", alt="V")
    'reference'
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def translate_cds(record: SeqRecord) -> Seq:
    """Translate the CDS feature of a GenBank record, stopping at the first stop codon.

    Uses the record's own CDS location (no hard-coded slice).
    Raises ValueError if the record has no CDS feature.
    """
    cds = next((f for f in record.features if f.type == "CDS"), None)
    if cds is None:
        raise ValueError(f"{record.id} has no CDS feature")
    return cds.extract(record.seq).translate(to_stop=True)


def check_residue(protein: Seq | str, index: int, ref: str, alt: str) -> str:
    """Return "reference" or "variant" for the residue at `index` (0-based, Met = 0).

    Raises AssertionError if the residue is neither `ref` nor `alt`, so a wrong
    numbering offset fails loudly instead of reporting a healthy result.
    """
    residue = str(protein[index])
    if residue == ref:
        return "reference"
    if residue == alt:
        return "variant"
    raise AssertionError(f"index {index} is {residue!r}, expected {ref!r} or {alt!r}")
