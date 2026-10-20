# 02 HBB variant validation

Notebook: [02_hbb_variant_validation.ipynb](../notebooks/02_hbb_variant_validation.ipynb)

## Question

Can I go from the HBB mRNA record to the protein and test programmatically for the sickle cell
substitution (Glu6Val)?

## Method

1. Fetch `NM_000518` from NCBI nucleotide (GenBank format) with Biopython `Entrez.efetch`.
2. Take the CDS location from the record's `CDS` feature instead of a hard-coded slice, extract it
   and translate with `to_stop=True`.
3. Compare the residue at index 6 with Val (`V`).

## Data

- Record used: `NM_000518.5`, *Homo sapiens* hemoglobin subunit beta (HBB), mRNA, the messenger RNA 
  transcribed from the HBB gene located on human chromosome 11 at position 11p15.5.
- Fetched live at run time. Nothing is saved in `data/raw/`, so the version is recorded here and in
  [SOURCES.md](../data/SOURCES.md).

## Result

- Translated protein begins `MVHLTPEEKSAVTALWGKVNVDEVGG…`.
- Index 6 is `E`, so the output was **"Healthy hemoglobin detected"**. This is the expected result,
  because the reference sequence is the normal (HbA) allele.

## Numbering: Glu6Val vs p.Glu7Val

- The translated sequence still contains the initiator Met at index 0.
- In the cell the Met is removed after translation, so the mature protein starts at Val. The legacy
  hemoglobin name "Glu6Val" numbers the mature protein, which makes residue 6 = index 6.
- HGVS counts the Met as residue 1, so the same variant is **p.Glu7Val**.
- Both names appear in the literature. Check which convention a source uses before comparing positions.

## Follow-up (done)

- CDS extraction and translation moved to `translate_cds` in [translation.py](../../../../src/cb/seq/translation.py),
  with `check_residue` for the variant call.
- `check_residue` raises if index 6 is neither `E` nor `V`, so a wrong numbering offset fails loudly.
- [test_translation.py](../../../../tests/seq/test_translation.py) covers the positive case (GAG → GTG in a
  synthetic CDS is reported as the variant), the reference case and the offset failure. Tests run offline.

## Positive control (run on the live record)

- Substituting GAG → GTG at codon index 6 of the real `NM_000518.5` CDS gives `MVHLTPVEKS…`
  (E7 → V in place), and the check prints **"Sickle cell variant detected: Glu6Val"**.
- The reference run still prints "Healthy hemoglobin detected", so both branches now give the right answer.
