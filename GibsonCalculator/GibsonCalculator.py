import tkinter as tk
import csv

class DNA_Ligation_GUI:
    def __init__(self, master):
        self.master = master
        master.title("DNA Ligation Calculator")

        # Create labels and input fields for DNA size and concentration


        self.upstream_label = tk.Label(master, text="Size (bp)")
        self.upstream_label.grid(row=0, column=1, sticky=tk.W)

        self.conc_label = tk.Label(master, text="Concentration (ng/ul):")
        self.conc_label.grid(row=0, column=2, sticky=tk.W)

        self.upstream_label = tk.Label(master, text="Upstream DNA")
        self.upstream_label.grid(row=1, column=0, sticky=tk.W)
        self.upstream_size = tk.Entry(master)
        self.upstream_size.grid(row=1, column=1)
        self.upstream_size.insert(0, "500")


        self.upstream_conc = tk.Entry(master)
        self.upstream_conc.grid(row=1, column=2)

        self.apramycin_label = tk.Label(master, text="Apramycin resistance gene DNA")
        self.apramycin_label.grid(row=2, column=0, sticky=tk.W)
        self.apramycin_size = tk.Entry(master)
        self.apramycin_size.grid(row=2, column=1)
        self.apramycin_size.insert(0, "1400")
        self.apramycin_conc = tk.Entry(master)
        self.apramycin_conc.grid(row=2, column=2)

        self.downstream_label = tk.Label(master, text="Downstream DNA")
        self.downstream_label.grid(row=3, column=0, sticky=tk.W)
        self.downstream_size = tk.Entry(master)
        self.downstream_size.grid(row=3, column=1)
        self.downstream_size.insert(0, "500")
        self.downstream_conc = tk.Entry(master)
        self.downstream_conc.grid(row=3, column=2)

        self.vector_label = tk.Label(master, text="Vector DNA")
        self.vector_label.grid(row=4, column=0, sticky=tk.W)
        self.vector_size = tk.Entry(master)
        self.vector_size.grid(row=4, column=1)
        self.vector_size.insert(0, "4466")
        self.vector_conc = tk.Entry(master)
        self.vector_conc.grid(row=4, column=2)

        # Create labels and input field for molar ratio
        self.ratio_label = tk.Label(master, text="Molar ratio:")
        self.ratio_label.grid(row=5, column=0, sticky=tk.W)
        self.ratio = tk.Entry(master)
        self.ratio.grid(row=5, column=1)
        self.ratio.insert(0, "3")

        self.vector_mass_label = tk.Label(master, text="Vector Mass to add (ng)")
        self.vector_mass_label.grid(row=6, column=0, sticky=tk.W)
        self.vector_mass = tk.Entry(master)
        self.vector_mass.grid(row=6, column=1)
        self.vector_mass.insert(0, "50")
        self.name_label = tk.Label(master, text="Knockout name")
        self.name_label.grid(row=7, column=0, sticky=tk.W)
        self.name = tk.Entry(master)
        self.name.grid(row=7, column=1)

        # Create button to calculate required volumes
        self.calculate_button = tk.Button(master, text="Calculate", command=self.calculate_volumes)
        self.calculate_button.grid(row=9, column=1)

        # Create export button to export calculated values to pdf
        self.export_button = tk.Button(master, text="Export", command=self.export)
        self.export_button.grid(row=9, column=2)

        # Create labels to display calculated volumes
        self.upstream_vol_label = tk.Label(master, text="Upstream DNA volume (ul):")
        self.upstream_vol_label.grid(row=10, column=0, sticky=tk.W)
        self.upstream_vol = tk.Label(master, text="")
        self.upstream_vol.grid(row=10, column=1)

        self.apramycin_vol_label = tk.Label(master, text="Apramycin resistance gene DNA volume (ul):")
        self.apramycin_vol_label.grid(row=11, column=0, sticky=tk.W)
        self.apramycin_vol = tk.Label(master, text="")
        self.apramycin_vol.grid(row=11, column=1)

        self.downstream_vol_label = tk.Label(master, text="Downstream DNA volume (ul):")
        self.downstream_vol_label.grid(row=12, column=0, sticky=tk.W)
        self.downstream_vol = tk.Label(master, text="")
        self.downstream_vol.grid(row=12, column=1)

        self.vector_vol_label = tk.Label(master, text="Vector DNA volume (ul):")
        self.vector_vol_label.grid(row=13, column=0, sticky=tk.W)
        self.vector_vol = tk.Label(master, text="")
        self.vector_vol.grid(row=13, column=1)

    def calculate_volumes(self):
        # Get values from input fields
        upstream_size = float(self.upstream_size.get())
        upstream_conc = float(self.upstream_conc.get())
        apramycin_size = float(self.apramycin_size.get())
        apramycin_conc = float(self.apramycin_conc.get())
        downstream_size = float(self.downstream_size.get())
        downstream_conc = float(self.downstream_conc.get())
        vector_mass = float(self.vector_mass.get())
        vector_size = float(self.vector_size.get())
        vector_conc = float(self.vector_conc.get())
        ratio = float(self.ratio.get())
        name = self.name.get()

        # Calculate required volumes based on molar ratio and DNA concentration required
        #mass insert (g) = desired insert/vector molar ratio x mass of vector (g) x ratio of insert to vector lengths
        upstream_mass = vector_mass *  upstream_size / vector_size * ratio
        upstream_vol = upstream_mass / upstream_conc
        self.upstream_vol.config(text="{:.2f}".format(upstream_vol))

        apramycin_mass = vector_mass *  apramycin_size / vector_size * ratio
        apramycin_vol = apramycin_mass / apramycin_conc
        self.apramycin_vol.config(text="{:.2f}".format(apramycin_vol))

        downstream_mass = vector_mass *  downstream_size / vector_size * ratio
        downstream_vol = downstream_mass / downstream_conc
        self.downstream_vol.config(text="{:.2f}".format(downstream_vol))


        vector_vol = vector_mass / vector_conc
        self.vector_vol.config(text="{:.2f}".format(vector_vol))


    def save_to_csv(self):
        # Get values from input fields
        upstream_size = float(self.upstream_size.get())
        upstream_conc = float(self.upstream_conc.get())
        apramycin_size = float(self.apramycin_size.get())
        apramycin_conc = float(self.apramycin_conc.get())
        downstream_size = float(self.downstream_size.get())
        downstream_conc = float(self.downstream_conc.get())
        vector_mass = float(self.vector_mass.get())
        vector_size = float(self.vector_size.get())
        vector_conc = float(self.vector_conc.get())
        ratio = float(self.ratio.get())
        name = self.name.get()

        # Calculate required volumes based on molar ratio and DNA concentration required
        #mass insert (g) = desired insert/vector molar ratio x mass of vector (g) x ratio of insert to vector lengths
        upstream_mass = vector_mass *  upstream_size / vector_size * ratio
        upstream_vol = upstream_mass / upstream_conc
        apramycin_mass = vector_mass *  apramycin_size / vector_size * ratio
        apramycin_vol = apramycin_mass / apramycin_conc
        downstream_mass = vector_mass *  downstream_size / vector_size * ratio
        downstream_vol = downstream_mass / downstream_conc
        vector_vol = vector_mass / vector_conc
        gibson_vol = upstream_vol + apramycin_vol + downstream_vol + vector_vol

        # Write data to CSV file
        with open('GibsonCalculator.csv', mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Name', 'Upstream (ul)', 'Apramycin (ul)', 'Downstream (ul)', 'Vector (ul)', 'Gibson (ul)'])
            # Write data row
            writer.writerow([name, upstream_vol, apramycin_vol, downstream_vol, vector_vol, gibson_vol])
            # Write blank row
            writer.writerow([])

    def export(self):
        self.save_to_csv()

root = tk.Tk()
gui = DNA_Ligation_GUI(root)
root.mainloop()
