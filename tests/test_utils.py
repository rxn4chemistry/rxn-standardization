from rxn_standardization.utils import (
    process_input,
    process_token,
    remove_stereochemistry,
)


def test_process_token() -> None:
    tokens = [
        "[Ag+]",
        "[Rh+3]",
        "[O-2]",
        "[Cl-]",
    ]
    expected = [
        "[ Ag + ]",
        "[ Rh +3 ]",
        "[ O -2 ]",
        "[ Cl - ]",
    ]

    for tok, exp in zip(tokens, expected):
        assert process_token(tok) == exp


def test_process_input() -> None:
    smiles = [
        "[O-] [Cl+3] ( [O-] ) ( [O-] ) O [Ag]",
        "[Cu+] ~ [Cu+] ~ [O-2]",
    ]
    expected = [
        "[ O - ] [ Cl +3 ] ( [ O - ] ) ( [ O - ] ) O [ Ag ]",
        "[ Cu + ] ~ [ Cu + ] ~ [ O -2 ]",
    ]

    for smi, exp in zip(smiles, expected):
        assert process_input(smi) == exp


def test_remove_stereochemistry() -> None:
    smiles = [
        "Cc1ccc([C@@H]2NNC(=O)[C@H]2NC(=O)c2ccccc2)cc1",
        "C[C@@H]1C2CC3(CC(=O)O2)[C@H](C)C[C@@H](O)[C@]3(O)[C@]1(C)CO",
    ]
    expected = [
        "Cc1ccc(C2NNC(=O)C2NC(=O)c2ccccc2)cc1",
        "CC1C2CC3(CC(=O)O2)C(C)CC(O)C3(O)C1(C)CO",
    ]

    for smi, exp in zip(smiles, expected):
        assert remove_stereochemistry(smi) == exp
