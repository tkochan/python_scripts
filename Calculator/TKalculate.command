#!/usr/bin/env python3

import tkinter as tk
from tkinter import ttk
# Create a new Tkinter window
root = tk.Tk()

root.title("Math is hard!")

# Create labels for the input boxes
label_v1 = tk.Label(root, text="Volume of stock to add:")
label_c1 = tk.Label(root, text="Concentration of stock:")
label_v2 = tk.Label(root, text="Final Volume:")
label_c2 = tk.Label(root, text="Final Concentration:")

# Add the labels and input boxes to the window
label_v1.grid(row=0, column=0)
entry_v1 = tk.Entry(root)
entry_v1.grid(row=0, column=1)
v1_dropdown = ttk.Combobox(root, values=["L", "mL", "uL"])
v1_dropdown.grid(row=0, column=2)
v1_dropdown.current(0)

label_c1.grid(row=1, column=0)
entry_c1 = tk.Entry(root)
entry_c1.grid(row=1, column=1)
c1_dropdown = ttk.Combobox(root, values=["g/L", "mg/L", "ug/L","g/mL", "mg/mL", "ug/mL","g/uL", "mg/mL", "ug/uL", "%w/v"])
c1_dropdown.grid(row=1, column=2)
c1_dropdown.current(0)

label_v2.grid(row=2, column=0)
entry_v2 = tk.Entry(root)
entry_v2.grid(row=2, column=1)
v2_dropdown = ttk.Combobox(root, values=["L", "mL", "uL"])
v2_dropdown.grid(row=2, column=2)
v2_dropdown.current(0)

label_c2.grid(row=3, column=0)
entry_c2 = tk.Entry(root)
entry_c2.grid(row=3, column=1)
c2_dropdown = ttk.Combobox(root, values=["g/L", "mg/L", "ug/L","g/mL", "mg/mL", "ug/mL","g/uL", "mg/mL", "ug/uL", "%w/v"])
c2_dropdown.grid(row=3, column=2)
c2_dropdown.current(0)

# Create an "Enter" button and a label to display the result
result_labela = tk.Label(root, text="")
result_labela.grid(row=4, columnspan=6)
#wv_unit_factor = 10
#Create dictionaries for volume and concentration conversion factors
#volume_unit_factors = {"L": 1, "mL": 1e-3, "uL": 1e-6}
#concentration_unit_factors = {"g/L": 1, "mg/L": 1e-3, "ug/L": 1e-6,
#                        "g/mL": 1e3, "mg/mL": 1, "ug/mL": 1e-3,
#                        "g/uL": 1e6, "mg/uL": 1e3, "ug/uL": 1,  "%w/v": wv_unit_factor }
#Create dictionaries for volume and concentration conversion factors
volume_unit_factors = {"L": 1, "mL": 1e-3, "uL": 1e-6}

# Define concentration unit factors based on user input


concentration_unit_factors = {"g/L": 1, "mg/L": 1e-3, "ug/L": 1e-6,
                        "g/mL": 1e3, "mg/mL": 1, "ug/mL": 1e-3,
                        "g/uL": 1e6, "mg/uL": 1e3, "ug/uL": 1, "%w/v":10}

def calculate_missing_input():
    # Get the values from the input boxes
    v1 = entry_v1.get()
    c1 = entry_c1.get()
    v2 = entry_v2.get()
    c2 = entry_c2.get()
    v1_unit = v1_dropdown.get()
    c1_unit = c1_dropdown.get()
    v2_unit = v2_dropdown.get()
    c2_unit = c2_dropdown.get()

    # Convert the input values to floats
    try:
        v1 = float(v1) if v1 else None
        c1 = float(c1) if c1 else None
        v2 = float(v2) if v2 else None
        c2 = float(c2) if c2 else None
    except ValueError:
        result_label.config(text="Error: Please enter valid numbers in at least 3 fields")
        error_button.config(state="normal")
        return

    # Convert the input values to a common unit
    v1_conv = v1 * volume_unit_factors[v1_unit] if v1 else None
    c1_conv = c1 * concentration_unit_factors[c1_unit] if c1 else None
    v2_conv = v2 * volume_unit_factors[v2_unit] if v2 else None
    c2_conv = c2 * concentration_unit_factors[c2_unit] if c2 else None

    # Calculate the missing value
    if not v1 and c1 and v2 and c2:
        V1_calc = v2_conv*c2_conv/c1_conv
        if V1_calc < 1e-3:
            V1_calc2 = V1_calc * 1000000
            result_labela.config(text=f"V1 = {V1_calc2:.3g} uL")
        elif V1_calc < 1:
            V1_calc2 = V1_calc * 1000
            result_labela.config(text=f"V1 = {V1_calc2:.3g} mL")
        elif V1_calc >= 1:
            result_labela.config(text=f"V1 = {V1_calc:.3g} L")
    elif v1 and not c1 and v2 and c2:
        result_labela.config(text=f"C1 = {v2_conv*c2_conv/v1_conv} mg/mL")


    elif v1 and c1 and not v2 and c2:
        V2_calc= v1_conv*c1_conv/c2_conv
        if V2_calc < 1e-3:
            V2_calc2 = V2_calc * 1000000
            result_labela.config(text=f"V2 = {V2_calc2:.3g} uL")
        elif V2_calc < 1:
            V2_calc2 = V2_calc * 1000
            result_labela.config(text=f"V2 = {V2_calc2:.3g} mL")
        elif V2_calc >= 1:
            result_labela.config(text=f"V2 = {V2_calc:.3g} L")
    elif v1 and c1 and v2 and not c2:
        result_labela.config(text=f"C2 = {(v1_conv*c1_conv/v2_conv)*1000} ug/mL")

#    error_button.config(state="disabled")
# Create a button to trigger the calculation
calculate_button1 = ttk.Button(root, text="Calculate", command=calculate_missing_input)
calculate_button1.grid(column=1, row=5, padx=5, pady=5, sticky="W")

# Create a button to clear the error message and allow the calculation to be re-run
#error_button = tk.Button(root, text="OK", command=lambda: [result_label.config(text=""), error_button.config(state="disabled")])
#error_button.grid(row=5, column=1)
#error_button.config(state="disabled")

#Create a title for the second calculator
calculator_label = ttk.Label(root, text="")
calculator_label.grid(column=0, row=7, padx=5, pady=5, sticky="E")
calculator_label2 = ttk.Label(root, text="")
calculator_label2.grid(column=1, row=8, padx=5, pady=5, sticky="E")

# Create a label for units of mass
mass_unit_label = ttk.Label(root, text="Mass Unit:")
mass_unit_label.grid(column=0, row=9, padx=5, pady=5, sticky="E")

# Create a label for units of volume
volume_unit_label = ttk.Label(root, text="Volume Unit:")
volume_unit_label.grid(column=0, row=10, padx=5, pady=5, sticky="E")

# Create an entry field for mass input
mass_entry = tk.Entry(root)
mass_entry.grid(column=1, row=9, padx=5, pady=5, sticky="W")

# Create an entry field for volume input
volume_entry = tk.Entry(root)
volume_entry.grid(column=1, row=10, padx=5, pady=5, sticky="W")


# Create a drop down menu for units of mass
mass_unit_var = tk.StringVar()
mass_unit_dropdown = ttk.Combobox(root, width=10, textvariable=mass_unit_var)
mass_unit_dropdown["values"] = ["µg", "mg", "g"]
mass_unit_dropdown.current(0)
mass_unit_dropdown.grid(column=2, row=9, padx=5, pady=5, sticky="W")

# Create a drop down menu for units of volume
volume_unit_var = tk.StringVar()
volume_unit_dropdown = ttk.Combobox(root, width=10, textvariable=volume_unit_var)
volume_unit_dropdown["values"] = ["mL", "L", "µL"]
volume_unit_dropdown.current(0)
volume_unit_dropdown.grid(column=2, row=10, padx=5, pady=5, sticky="W")

# Create a button to calculate the %w/v
def calculate():
    try:
        # Convert mass to grams
        mass = float(mass_entry.get())
        mass_unit = mass_unit_var.get()
        if mass_unit == "mg":
            mass /= 1000
        elif mass_unit == "µg":
            mass /= 1000000

        # Convert volume to milliliters
        volume = float(volume_entry.get())
        volume_unit = volume_unit_var.get()
        if volume_unit == "L":
            volume *= 1000
        elif volume_unit == "µL":
            volume /= 1000

        # Calculate %w/v
        percent_wv = (mass / volume) * 100
        result_label.config(text="%.2f%% w/v" % percent_wv)
    except ValueError:
        result_label.config(text="Invalid input")

# Create a button to calculate the %w/v
calculate_button = ttk.Button(root, text="Calculate", command=calculate)
calculate_button.grid(column=1, row=12, padx=5, pady=5, sticky="W")

# Create a label for the result
result_label = ttk.Label(root, text="")
result_label.grid(column=1, row=11, padx=5, pady=5, sticky="E")






# Start the Tkinter event loop
root.mainloop()
