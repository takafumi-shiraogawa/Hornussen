from pyscf import gto, scf, cc

mol = gto.M(atom='F 0 0 0', basis='def2-TZVP', spin=1)
hfe = scf.UHF(mol)
hfe.max_cycle = 20000
hfe.kernel(verbose=0)
mycc = cc.UCCSD(hfe).run()