from mom_family import *
for t in (0.787e-3, 1.575e-3, 2.54e-3, 3.175e-3):
    for (W, nm) in ((2.91, 'W2.91'), (5.0, 'W5.0')):
        run('disc 5.45 L8.72 %s, t=%.3f mm' % (nm, t*1e3), disc_slot(5.45, 8.72, W), t)
run('dogbone t=1.575', dogbone(5.45, 8.9, 1.5, 4.5, 1.5), 1.575e-3)
