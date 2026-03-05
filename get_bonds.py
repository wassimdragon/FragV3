from rdkit import Chem
smiles = [
    "C1=C(C(=O)NC(=O)N1)I",
    "C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23",
    "C1=CC=C(C=C1)I",
    "O=[N+](C1=C(I)NC(Br)=N1)[O-]"
]
for s in smiles:
    m = Chem.MolFromSmiles(s)
    m = Chem.AddHs(m)
    print(f"{s[:15]}... total bonds: {m.GetNumBonds()}")
