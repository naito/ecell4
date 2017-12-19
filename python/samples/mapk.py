from ecell4 import *


with species_attributes():
    K | {'N': 120}
    KK | {'N': 30}
    PP | {'N': 30}

kon1, koff1, kcat1, kon2, koff2, kcat2 = (
    4.483455086786913e-20, 1.35, 1.5,
    9.299017957780264e-20, 1.73, 15.0)

with reaction_rules():
    (K + KK == K_KK | (kon1, koff1)
        > Kp + KK | kcat1
        == Kp_KK | (kon2, koff2)
        > Kpp + KK | kcat2)

    (Kpp + PP == Kpp_PP | (kon1, koff1)
        > Kp + PP | kcat1
        == Kp_PP | (kon2, koff2)
        > K + PP | kcat2)

m = get_model()

y0 = {}
for sp in m.species_attributes():
    y0[sp.serial()] = sp.get_attribute('N')
print(y0)

obs = run_simulation(4, volume=1e-18, y0=y0, model=m, solver='gillespie', rndseed=0, return_type='observer')
obs.save('simple.csv')
