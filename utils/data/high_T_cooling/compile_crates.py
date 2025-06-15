import glob

elements = ["carbon","nitrogen","oxygen","neon","magnesium","silicon","sulphur","iron"]
atomic_number = [6, 7, 8, 10, 12, 14, 16, 26]

for an,el in zip(atomic_number,elements):
    cooling_table = np.zeros((an+1,121)) # C is row oriented to store ions in rows
    for i in range(1,an+2):
        dat = np.loadtxt(f"./{el.upper()}/{el}_{i}_cooling.dat")
        cooling_table[i-1,:] = np.log10(dat[:,2])

    np.savetxt(f"./{el.upper()}/all_cool.dat",cooling_table)

