from pyscf import gto, scf, cc

mol = gto.M(atom='B 0 0 0', basis='def2-TZVP', spin=1)
hfe = scf.UKS(mol)
hfe.xc = 'b3lyp'
hfe.max_cycle = 20000
hfe.kernel()