# exported from PySB model 'model'

import numpy as np
from pysb.integrate import Solver
from pysb.simulator import ScipyOdeSimulator
import pylab as pl
import matplotlib.pyplot as plt
from pysb import Model, Monomer, Parameter, Expression, Compartment, Rule, Observable, Initial, MatchOnce, Annotation, ANY, WILD

Model()

Monomer('A', ['B', 'C'], {'C': ['off', 'on'], 'B': ['off', 'on']})
Monomer('C', ['C'], {'C': ['off', 'on']})
Monomer('B', ['A', 'C'], {'A': ['off', 'on'], 'C': ['off', 'on']})

Parameter('B_signal_A_0_2kc', 1.0)
Parameter('B_signal_A_1_2kc', 1.0)
Parameter('B_signal_A_2_2kc', 1.0)
Parameter('B_signal_A_3_2kc', 1.0)
Parameter('C_signal_A_0_2kc', 1.0)
Parameter('C_signal_A_1_2kc', 1.0)
Parameter('C_signal_C_0_2kc', 1.0)
Parameter('C_signal_C_1_2kc', 1.0)
Parameter('A_signal_B_0_2kc', 1.0)
Parameter('A_signal_B_1_2kc', 1.0)
Parameter('A_signal_B_2_2kc', 1.0)
Parameter('A_signal_B_3_2kc', 1.0)
Parameter('C_signal_B_0_2kc', 1.0)
Parameter('C_signal_B_1_2kc', 1.0)
Parameter('A_0_0', 0.0)
Parameter('A_1_0', 0.3333333333333333)
Parameter('A_2_0', 0.3333333333333333)
Parameter('A_3_0', 0.3333333333333333)
Parameter('C_0_0', 1.0)
Parameter('C_1_0', 0.0)
Parameter('B_0_0', 0.0)
Parameter('B_1_0', 0.0)
Parameter('B_2_0', 0.0)
Parameter('B_3_0', 1.0)

Observable('A_0_obs', A(B='on', C='on'))
Observable('A_1_obs', A(B='on', C='off'))
Observable('A_2_obs', A(B='off', C='on'))
Observable('A_3_obs', A(B='off', C='off'))
Observable('C_0_obs', C(C='on'))
Observable('C_1_obs', C(C='off'))
Observable('B_0_obs', B(A='on', C='on'))
Observable('B_1_obs', B(A='on', C='off'))
Observable('B_2_obs', B(A='off', C='on'))
Observable('B_3_obs', B(A='off', C='off'))
Observable('A_obs', A())
Observable('C_obs', C())
Observable('B_obs', B())

Rule('B_signal_A_0', B(A='on', C='on') + A(B='off') >> B(A='on', C='on') + A(B='on'), B_signal_A_0_2kc)
Rule('B_signal_A_1', B(A='on', C='off') + A(B='off') >> B(A='on', C='off') + A(B='on'), B_signal_A_1_2kc)
Rule('B_signal_A_2', B(A='off', C='on') + A(B='off') >> B(A='off', C='on') + A(B='on'), B_signal_A_2_2kc)
Rule('B_signal_A_3', B(A='off', C='off') + A(B='on') >> B(A='off', C='off') + A(B='off'), B_signal_A_3_2kc)
Rule('C_signal_A_0', C(C='on') + A(C='off') >> C(C='on') + A(C='on'), C_signal_A_0_2kc)
Rule('C_signal_A_1', C(C='off') + A(C='on') >> C(C='off') + A(C='off'), C_signal_A_1_2kc)
Rule('C_signal_C_0', C(C='on') + C(C='off') >> C(C='on') + C(C='on'), C_signal_C_0_2kc)
Rule('C_signal_C_1', C(C='off') + C(C='on') >> C(C='off') + C(C='off'), C_signal_C_1_2kc)
Rule('A_signal_B_0', A(B='on', C='on') + B(A='off') >> A(B='on', C='on') + B(A='on'), A_signal_B_0_2kc)
Rule('A_signal_B_1', A(B='on', C='off') + B(A='on') >> A(B='on', C='off') + B(A='off'), A_signal_B_1_2kc)
Rule('A_signal_B_2', A(B='off', C='on') + B(A='on') >> A(B='off', C='on') + B(A='off'), A_signal_B_2_2kc)
Rule('A_signal_B_3', A(B='off', C='off') + B(A='on') >> A(B='off', C='off') + B(A='off'), A_signal_B_3_2kc)
Rule('C_signal_B_0', C(C='on') + B(C='off') >> C(C='on') + B(C='on'), C_signal_B_0_2kc)
Rule('C_signal_B_1', C(C='off') + B(C='on') >> C(C='off') + B(C='off'), C_signal_B_1_2kc)

Initial(A(B='on', C='on'), A_0_0)
Initial(A(B='on', C='off'), A_1_0)
Initial(A(B='off', C='on'), A_2_0)
Initial(A(B='off', C='off'), A_3_0)
Initial(C(C='on'), C_0_0)
Initial(C(C='off'), C_1_0)
Initial(B(A='on', C='on'), B_0_0)
Initial(B(A='on', C='off'), B_1_0)
Initial(B(A='off', C='on'), B_2_0)
Initial(B(A='off', C='off'), B_3_0)


time = [0, 1000]
model_solver = ScipyOdeSimulator(model, tspan=time, integrator='vode', integrator_options={'atol': 1e-12, 'rtol': 1e-12})

print
for i, each in enumerate(model.odes):
    print i, each
print

solved = model_solver.run()

for each in solved.species:
    print each
