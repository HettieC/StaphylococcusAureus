using PyCall
using Conda

Conda.add("equilibrator-api"; channel="conda-forge") # note: use a dash and not an underscore in "equilibrator-api"
eq = pyimport("equilibrator_api") # note: use an underscore and not a dash in "equilibrator_api"

using eQuilibrator
using Unitful

temp = 30u"°C"
i_strength = 150.0u"mM"
ph = 7.9
pmg = 2.0

equilibrator = eQuilibrator.Equilibrator(;pH=ph, pMg=pmg, temperature=temp, ionic_strength=i_strength);

rxn_string = "bigg.metabolite:atp + bigg.metabolite:h2o = bigg.metabolite:adp + bigg.metabolite:pi"

physiological_dg_prime(equilibrator, rxn_string)
# -46.26 ± 0.3 kJ mol^-1

standard_dg_prime(equilibrator, rxn_string)
# -29.14 ± 0.3 kJ mol^-1

dg_prime(equilibrator, rxn_string) # equilibrator_api default abundances/concentrations
# -29.14 ± 0.3 kJ mol^-1

concens = Dict("bigg.metabolite:atp"=>1u"mM", "bigg.metabolite:adp"=>100u"μM", "bigg.metabolite:pi"=>0.005u"M")
dg_prime(equilibrator, rxn_string; concentrations=concens) # user specified concentrations
# -47.98 ± 0.3 kJ mol^-1

ln_reversibility_index(equilibrator, rxn_string)
# -12.447 ± 0.082
