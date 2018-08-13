"""Mass fraction compositions for reactor materials.
"""

# fuel fraction in cermet 

rhos = {'UO2': 10.97, 'W' : 19.3, 'UN' : 11.3}

mats = {'UO2' : {8016  : 1.5204e-01,
                 92000 : 8.4796e-01
                },

        'UN'  : {92000 : 0.94441, 
                 7015 : 0.05559
                },

        'CO2' : {6000 : 2.7291e-01,
                 8016 : 7.2709e-01
                },

        'W'   : {74180 : 1.1746e-03,
                 74182 : 2.6227e-01,
                 74183 : 1.4241e-01,
                 74184 : 3.0658e-01,
                 74186 : 2.8757e-01
                },

        'H2O' : {1001 : 1.1187e-01,
                 1002 : 2.5713e-05,
                 8016 : 8.8811e-01,
                }
       }


def enrich_fuel(enrich, fuel):
    """Enrich Uranium fuel and mix with specified compound.

    Arguments:
    ----------
        enrich (float): U-235 mass conc.
        fuel (dict) (optional): mass composition of fuel
    Returns:
    --------
        fuel_comp (dict): isotopic mass vector of Uranium fuel compound.
    """
    fuel_comp = {}
    # build fuel comp. starting with bounded compound
    for isotope in mats[fuel]:
        if isotope == 92000:
            # add enriched Uranium to fuel composition
            fuel_comp.update({92235 : mats[fuel][92000]*enrich,
                              92238 : mats[fuel][92000]*(1-enrich) })
        else:
            fuel_comp.update({isotope : mats[fuel][isotope]})
    
    mmO = 32
    mmN = 14.0067
    mmU25 = 235.04
    mmU28 = 238.05

    if fuel == 'UN':
        mmfuel = mmN + 1 / ((enrich / mmU25) + ((1-enrich) / mmU28))
    if fuel == 'UO2':
        mmfuel = mmO + 1 / ((enrich / mmU25) + ((1-enrich) / mmU28))
    
    return fuel_comp, mmfuel
