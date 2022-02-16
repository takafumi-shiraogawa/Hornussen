from pyscf import gto, scf, cc

mol = gto.M(atom='N 0 0 0', basis='def2-TZVP', spin=3)
hfe = scf.UHF(mol)
hfe.max_cycle = 20000
hfe.kernel()
mycc = cc.UCCSD(hfe).run()