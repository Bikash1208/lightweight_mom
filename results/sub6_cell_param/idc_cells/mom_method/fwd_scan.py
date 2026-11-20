from fwd_beta import *
psi0 = layout_psi(1.10, 1.10)
ref = peaks(psi0, 1.198, 1.298)
print('ref R %.2f L %.2f' % (ref['R'][0], ref['L'][0]))
for A in (1.10, 1.14, 1.18):
    for B in (0.90, 0.95, 1.00, 1.05):
        o = peaks(psi0, A, B)
        print('A %.2f B %.2f: R delta %+.2f th %.1f | L delta %+.2f th %.1f   (target R -0.58 th ~31, L -2.64 th ~26)' % (A, B, o['R'][0]-ref['R'][0], o['R'][1], o['L'][0]-ref['L'][0], o['L'][1]), flush=True)
