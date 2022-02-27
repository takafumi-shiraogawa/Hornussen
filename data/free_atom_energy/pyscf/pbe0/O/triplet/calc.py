from pyscf import gto, scf, cc

mol = gto.M(atom='O 0 0 0', basis='def2-TZVP', spin=2)
hfe = scf.UKS(mol)
hfe.xc = 'pbe0'
hfe.max_cycle = 20000
hfe.kernel()