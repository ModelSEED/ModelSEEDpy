\* Problem: Unknown *\

Maximize
 obj: + R_3OAS140 - R_3OAS140_reverse_60acb

Subject To
 cpd11492: + rxn2_reverse_bd47f - rxn2 + R_3OAS140_reverse_60acb
 - R_3OAS140 = 0
 cpd00067: - rxn2_reverse_bd47f + rxn2 + R_3OAS140_reverse_60acb
 - R_3OAS140 = 0
 cpd11468: + rxn2_reverse_bd47f - rxn2 + R_3OAS140_reverse_60acb
 - R_3OAS140 = 0
 cpd00011: + rxn2_reverse_bd47f - rxn2 - R_3OAS140_reverse_60acb
 + R_3OAS140 = 0
 cpd11628: - rxn2_reverse_bd47f + rxn2 - R_3OAS140_reverse_60acb
 + R_3OAS140 = 0
 cpd15377: - R_3OAS140_reverse_60acb + R_3OAS140 = 0
 cpd11492_conc: + cpd11492 >= 0
 cpd00067_conc: + cpd00067 >= 0
 cpd11468_conc: + cpd11468 >= 0
 cpd00011_conc: + cpd00011 >= 0
 cpd11628_conc: + cpd11628 >= 0
 cpd15377_conc: + cpd15377 >= 0
 R_3OAS140_kinetics: - R_3OAS140_reverse_60acb + R_3OAS140
 = 1.4455859694423e-16
 rxn2_kinetics: - rxn2_reverse_bd47f + rxn2 = 1.43228294999153e-17

Bounds
 0 <= R_3OAS140 <= 1000
 R_3OAS140_reverse_60acb = 0
 0 <= rxn2 <= 1000
 rxn2_reverse_bd47f = 0

End
