k.<r> = QuadraticField(165)
print("disc", k.discriminant(), "int basis", k.ring_of_integers().basis())
for P in k.primes_above(2): print("P2:", P, "e",P.ramification_index(),"f",P.residue_class_degree(), "gens", P.gens(), "unif", k.uniformizer(P))
om = (1+r)/2
print("om^2 =", om^2, " == om+41?", om^2 == om+41)
print("class number", k.class_number())
