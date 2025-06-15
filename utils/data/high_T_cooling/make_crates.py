import os
import numpy as np

elements = ["carbon","nitrogen","oxygen","neon","magnesium","silicon","sulphur","iron"]
atomic_number = [6, 7, 8, 10, 12, 14, 16, 26]
all_elements = ["HELIUM", "LITHIUM", "MAGNESIUM", "BERYLLIUM", "CARBON", "NITROGEN", "SULPHUR", "OXYGEN", "NEON", "SILICON", "IRON", "BORON", "FLUORINE", "SODIUM", "ALUMINIUM", "PHOSPHORUS", "CHLORINE", "ARGON", "POTASSIUM", "CALCIUM", "SCANDIUM", "TITANIUM", "VANADIUM", "CHROMIUM", "MANGANESE", "COBALT",  "NICKEL", "COPPER", "ZINC"]
pos_in_file = [8, 9, 10, 12, 14, 16, 18, 28]

eden = 0.0

def make_input(el,an,istate,eden):

    to_write = [
        "constant temperature 4.2 vary\n",
        "grid 3 9 0.05\n",
        "stop zone 1\n",
        "set dr 0\n",
        "hden -15.0\n",
        f"set eden {eden} log\n",
        'save cooling each "cooling.dat" no hash\n',
        "element hydrogen ionization 1.0 0.0\n",
        "no ionization reevaluation\n",
        "no collisional ionization\n",
        "no photoionization\n",
        "no charge transfer\n",
        "no auger effect\n",
        "no radiation pressure\n",
        "no molecules\n",
        #"no free free\n"
        #"no compton effect\n"
        #"no induced processes\n"
        f"element {el} abundance 15.0\n",
    ]

    # Set the ionization state of the metal
    zero_list = ["0.0"]*(an+1)
    zero_list[istate-1] = "1.0"
    to_write.append(f"element {el} ionization {' '.join(zero_list)}\n")
    for e in all_elements:
        if e.lower() != el:
            to_write.append(f"element {e} off\n")
    with open("cloudy_tmp.in","w") as f:
        for l in to_write:
            f.write(l)

# # Loop over elements
for el in elements:
    os.mkdir(el.upper())


counter = 0
for el,j in zip(elements,atomic_number):

    #Run the simulations
    for i in range(j+1):
        make_input(el,j,i+1,eden)
        os.system("../cloudy/source/cloudy.exe -r cloudy_tmp")
        os.system(f"mv cloudy_tmp.in {el}_{i+1}.in")
        os.system(f"mv cooling.dat {el}_{i+1}_cooling.dat")

    counter += 1

    os.system(f"mv {el}* ./{el.upper()}")


