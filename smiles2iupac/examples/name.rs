fn main() {
    for smi in ["CCCC","CC(C)C","CCCCO","CC(C)O","CC(=O)C","CCC=O","CCC(=O)O","CCN",
                "CCCC=C","CC#C","ClCCCl","OCCO","CC(=O)CC(=O)C","CCOC","FC(F)(F)C","CC(C)CC(C)C"] {
        match smiles2iupac::smiles_to_iupac(smi) {
            Ok(n) => println!("{smi:16} {n}"),
            Err(e) => println!("{smi:16} ERR {e}"),
        }
    }
}
