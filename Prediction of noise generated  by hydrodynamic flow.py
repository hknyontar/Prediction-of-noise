import math
import logging

def calculate_valve_flow():
    logging.basicConfig(level=logging.INFO)

    # User Inputs
    nominal_valve_size = float(input("Nominal Valve Size : "))
    rated_cv = float(input("Rated Cv: "))
    required_cv = float(input("Required Cv: "))
    seat_diameter = float(input("Seat Diameter (mm): ")) / 1000
    fl = float(input("Liquid pressure recovery factor (FL): "))
    fd = float(input("Valve style modifier (Fd): "))

    #Pipe
    inlet_pipe_diameter = float(input("Internal pipe diameter (mm): ")) / 1000
    pipe_wall_thickness = float(input("Pipe wall thickness (mm): ")) / 1000
    speed_of_sound_pipe = float(input("Speed of sound in pipe (m/s): "))
    density_pipe_material = float(input("Density of pipe material (kg/m3): "))

    #Other
    speed_of_sound_air = float(input("Speed of sound in air (m/s): "))
    density_air = float(input("Density of air (kg/m3): "))

    #Calculation Inputs
    mass_flow_rate = float(input("Mass Flow Rate (kg/s): "))
    inlet_pressure = float(input("Valve inlet absolute pressure (bar): ")) * 1e5
    outlet_pressure = float(input("Valve outlet absolute pressure (bar): ")) * 1e5
    vapour_pressure = float(input("Vapour pressure of liquid (bar): ")) * 1e5
    density_liquid = float(input("Density of liquid (kg/m3): "))
    speed_of_sound_liquid = float(input("Speed of sound in liquid (m/s): "))

    n34 = float(input("Numerical constant (N34): "))
    n14 = float(input("Numerical constant (N14): "))
    rw = float(input("Acoustic power ratio (rw): "))
    fi = float(input("Frequency (fi): "))

    # Calculation
    differential_pressure_ratio = (inlet_pressure - outlet_pressure) / (inlet_pressure - vapour_pressure)
    #Pressure differential for Uvc calculation
    dp_c1 = differential_pressure_ratio * (inlet_pressure - vapour_pressure)
    dp_c2 = (fl ** 2) * (inlet_pressure - vapour_pressure)

    #Differential pressure ratio of incipient cavitation noise
    xfz = 0.9 / math.sqrt(1 + 3 * fd * math.sqrt(required_cv / (n34 * fl)))

    #Differential pressure ratio corrected for inlet pressure
    xfz_p1 = xfz * (600000 / inlet_pressure) ** 0.125

    #Jet diameter
    djet = n14 * fd * math.sqrt(required_cv * fl)

    #Vena contracta velocity  
    uvc = (1 / fl) * math.sqrt((2 * dp_c1) / density_liquid)

    #Mechanical stream power
    wm = (mass_flow_rate * (uvc ** 2) * (fl ** 2)) / 2

    #Flow condition
    dp_cond1 = inlet_pressure - outlet_pressure
    dp_cond2 = xfz_p1 * (inlet_pressure - vapour_pressure)
    flow_condition = "Turbulent" if dp_cond1 < dp_cond2 else "Cavitating"

    if flow_condition == "Turbulent":
        eta_turbulent = (10 ** -4) * (uvc / speed_of_sound_liquid)    #Acoustic efficiency factor (turbulent)
        wa = eta_turbulent * wm * rw               #Sound power (turbulent region)
        Nstr = (((0.036) * (fl**2) * (required_cv) * (fd**0.75)) / ((n34) * (xfz_p1**1.5) * ((nominal_valve_size)/1000) * (seat_diameter))) * ((1 / (inlet_pressure - vapour_pressure))**0.57)  #Strouhal number of jet
        fp_turb = (Nstr * uvc) / (djet) #Peak sound frequency (turbulent)
        fp_cav = 6 * fp_turb * (((1 - differential_pressure_ratio) / (1 - xfz_p1)) ** 2) * ((xfz_p1 / differential_pressure_ratio) ** 2.5) #Peak sound frequency (cavitating)
        fr = speed_of_sound_pipe / (math.pi * inlet_pipe_diameter) #Ring frequency
        TLfr = -10 - 10 * math.log10((speed_of_sound_pipe * density_pipe_material * pipe_wall_thickness) / (speed_of_sound_air * density_air * inlet_pipe_diameter)) #Transmission loss at ring frequency
        delta_TLfp_turb = -20 * math.log10((fr / fp_turb) + ((fp_turb / fr)**1.5)) #Transmission loss corrected for fp,turb
        Tlturb = TLfr + delta_TLfp_turb #Overall transmission loss for turbulent flow
        Lpi_internal = 10 * math.log10(3.2 * 10**9 * wa * density_liquid * speed_of_sound_liquid / inlet_pipe_diameter**2)
        Lp_Ae_1m_turb = Lpi_internal + Tlturb - 10 * math.log10((inlet_pipe_diameter + (2 * pipe_wall_thickness) + 2) / (inlet_pipe_diameter + (2 * pipe_wall_thickness))) #External sound pressure level (turbulent)
        Fturb_fi = -10 * math.log10((0.25 * ((fi / fp_turb)**3) + ((fi / fp_turb)**-1))) - 3.1 #Frequency distribution function (turbulent)
        Fcav_fi = -10 * math.log10((0.25 * ((fi / fp_cav)**1.5) + ((fi / fp_cav)**-1.5))) - 3.5 #Frequency distribution function (cavitating)
        Lpi_fi_turb = Lpi_internal + Fturb_fi #Internal sound pressure level at fi (turbulent)
        delta_TLf_fi = -20 * math.log10((fr / fi) + ((fi / fr)**1.5)) #Transmission loss corrected for fi
        TL_fi = TLfr + delta_TLf_fi #Transmission loss at fi
        Lpe_1m_fi = Lpi_fi_turb + TL_fi - 10 * math.log10((inlet_pipe_diameter + (2 * pipe_wall_thickness) + 2) / (inlet_pipe_diameter + (2 * pipe_wall_thickness))) #External sound pressure level at fi

        #Results
        print(f"\n--- Results ---")
        print(f"Flow condition: {flow_condition}")
        print(f"Differential pressure ratio: {differential_pressure_ratio:.6f}")
        print(f"Pressure drop (dp_c1): {dp_c1:.3f} Pa")
        print(f"Pressure drop (dp_c2): {dp_c2:.3f} Pa")
        print(f"Critical pressure ratio (xfz): {xfz:.6f}")
        print(f"Corrected critical pressure ratio (xfz_p1): {xfz_p1:.6f}")
        print(f"Jet diameter (djet): {djet:.6f} m")
        print(f"Velocity coefficient (uvc): {uvc:.6f} m/s")
        print(f"Kinetic power per unit mass (wm): {wm:.6f} W/kg")
        print(f"Condition pressure drop 1 (dp_cond1): {dp_cond1:.3f} Pa")
        print(f"Condition pressure drop 2 (dp_cond2): {dp_cond2:.3f} Pa")
        print(f"Turbulent viscosity (eta_turbulent): {eta_turbulent:.6f}")
        print(f"Sound power in turbulent region (wa): {wa:.6f} W")
        print(f"Strouhal number of jet (Nstr): {Nstr:.6f}")
        print(f"Peak sound frequency (turbulent) (fp_turb): {fp_turb:.6f} Hz")
        print(f"Peak sound frequency (cavitating) (fp_cav): {fp_cav:.6f} Hz")
        print(f"Ring frequency (fr): {fr:.6f} Hz")
        print(f"Transmission loss at ring frequency (TLfr): {TLfr:.2f} dB")
        print(f"Transmission loss corrected for fp_turb (delta_TLfp_turb): {delta_TLfp_turb:.2f} dB")
        print(f"Overall transmission loss for turbulent flow (Tlturb): {Tlturb:.2f} dB")
        print(f"Internal sound pressure level (Lpi_internal): {Lpi_internal:.2f} dB")
        print(f"External sound pressure level (turbulent) (Lp_Ae_1m_turb): {Lp_Ae_1m_turb:.2f} dB")
        print(f"Frequency distribution function (turbulent) (Fturb_fi): {Fturb_fi:.2f} dB")
        print(f"Frequency distribution function (cavitating) (Fcav_fi): {Fcav_fi:.2f} dB")
        print(f"Internal sound pressure level at fi (turbulent) (Lpi_fi_turb): {Lpi_fi_turb:.2f} dB")
        print(f"Transmission loss corrected for fi (delta_TLf_fi): {delta_TLf_fi:.2f} dB")
        print(f"Transmission loss at fi (TL_fi): {TL_fi:.2f} dB")
        print(f"External sound pressure level at fi (Lpe_1m_fi): {Lpe_1m_fi:.2f} dB")

    else:
        eta_turbulent = (10 ** -4) * (uvc / speed_of_sound_liquid)
        eta_cavitating = (0.32) * (eta_turbulent) * (math.sqrt((dp_cond1 / dp_c1) * (1 / xfz_p1))) * (math.exp(5 * xfz_p1)) * (((1 - xfz_p1) / (1 - differential_pressure_ratio)) ** 0.5) * ((differential_pressure_ratio / xfz_p1) ** 5) * ((differential_pressure_ratio - xfz_p1) ** 1.5)
        wa = (eta_turbulent + eta_cavitating) * wm * rw               #Sound power (cavitating region)
        Nstr = (((0.036) * (fl**2) * (required_cv) * (fd**0.75)) / ((n34) * (xfz_p1**1.5) * ((nominal_valve_size)/1000) * (seat_diameter))) * ((1 / (inlet_pressure - vapour_pressure))**0.57)   #Strouhal number of jet
        fp_turb = (Nstr * uvc) / (djet) #Peak sound frequency (turbulent)
        fp_cav = 6 * fp_turb * (((1 - differential_pressure_ratio) / (1 - xfz_p1)) ** 2) * ((xfz_p1 / differential_pressure_ratio) ** 2.5) #Peak sound frequency (cavitating)
        fr = speed_of_sound_pipe / (math.pi * inlet_pipe_diameter) #Ring frequency
        TLfr = -10 - 10 * math.log10((speed_of_sound_pipe * density_pipe_material * pipe_wall_thickness) / (speed_of_sound_air * density_air * inlet_pipe_diameter)) #Transmission loss at ring frequency
        delta_TLfp_turb = -20 * math.log10((fr / fp_turb) + ((fp_turb / fr)**1.5)) #Transmission loss corrected for fp,turb
        Tlturb = TLfr + delta_TLfp_turb #Overall transmission loss for turbulent flow
        Lpi_internal = 10 * math.log10(3.2 * 10**9 * wa * density_liquid * speed_of_sound_liquid / inlet_pipe_diameter**2)
        Lp_Ae_1m_turb = Lpi_internal + Tlturb - 10 * math.log10((inlet_pipe_diameter + (2 * pipe_wall_thickness) + 2) / (inlet_pipe_diameter + (2 * pipe_wall_thickness))) #External sound pressure level (turbulent)
        Tlcav = Tlturb + 10 * math.log10( ((250) * (fp_cav**1.5) * (eta_cavitating)) / ((fp_turb**2) * (eta_turbulent + eta_cavitating))) #Overall transmission loss (cavitation)
        Lp_Ae_1m_cav = Lpi_internal + Tlcav - 10 * math.log10((inlet_pipe_diameter + (2 * pipe_wall_thickness) + 2) / (inlet_pipe_diameter + (2 * pipe_wall_thickness))) #External sound pressure level (cavitating)
        Fturb_fi = -10 * math.log10((0.25 * ((fi / fp_turb)**3) + ((fi / fp_turb)**-1))) - 3.1 #Frequency distribution function (turbulent)
        Fcav_fi = -10 * math.log10((0.25 * ((fi / fp_cav)**1.5) + ((fi / fp_cav)**-1.5))) - 3.5 #Frequency distribution function (cavitating)
        Lpi_fi_turb = Lpi_internal + Fturb_fi #Internal sound pressure level at fi (turbulent)
        Lpi_fi_cav = Lpi_internal + 10 * math.log10(((eta_turbulent / (eta_turbulent + eta_cavitating)) * (10**(0.1 * Fturb_fi)) + (eta_cavitating / (eta_cavitating + eta_turbulent)) * (10**(0.1 * Fcav_fi)))) #Internal sound pressure level at fi (cavitating)
        delta_TLf_fi = -20 * math.log10((fr / fi) + ((fi / fr)**1.5)) #Transmission loss corrected for fi
        TL_fi = TLfr + delta_TLf_fi #Transmission loss at fi
        Lpe_1m_fi = Lpi_fi_cav + TL_fi - 10 * math.log10((inlet_pipe_diameter + (2 * pipe_wall_thickness) + 2) / (inlet_pipe_diameter + (2 * pipe_wall_thickness))) #External sound pressure level at fi
  

        # Results
        print(f"\n--- Results ---")
        print(f"Flow condition: {flow_condition}")
        print(f"Differential pressure ratio: {differential_pressure_ratio:.6f}")
        print(f"Pressure drop (dp_c1): {dp_c1:.3f} Pa")
        print(f"Pressure drop (dp_c2): {dp_c2:.3f} Pa")
        print(f"Critical pressure ratio (xfz): {xfz:.6f}")
        print(f"Corrected critical pressure ratio (xfz_p1): {xfz_p1:.6f}")
        print(f"Jet diameter (djet): {djet:.6f} m")
        print(f"Velocity coefficient (uvc): {uvc:.6f} m/s")
        print(f"Kinetic power per unit mass (wm): {wm:.6f} W/kg")
        print(f"Condition pressure drop 1 (dp_cond1): {dp_cond1:.3f} Pa")
        print(f"Condition pressure drop 2 (dp_cond2): {dp_cond2:.3f} Pa")
        print(f"Turbulent viscosity (eta_turbulent): {eta_turbulent:.6e}")
        print(f"Cavitating viscosity (eta_cavitating): {eta_cavitating:.6e}")
        print(f"Sound power in cavitating region (wa): {wa:.6f} W")
        print(f"Strouhal number of jet (Nstr): {Nstr:.6f}")
        print(f"Peak sound frequency (turbulent) (fp_turb): {fp_turb:.6f} Hz")
        print(f"Peak sound frequency (cavitating) (fp_cav): {fp_cav:.6f} Hz")
        print(f"Ring frequency (fr): {fr:.6f} Hz")
        print(f"Transmission loss at ring frequency (TLfr): {TLfr:.2f} dB")
        print(f"Transmission loss corrected for fp_turb (delta_TLfp_turb): {delta_TLfp_turb:.2f} dB")
        print(f"Overall transmission loss for turbulent flow (Tlturb): {Tlturb:.2f} dB")
        print(f"Internal sound pressure level (Lpi_internal): {Lpi_internal:.2f} dB")
        print(f"External sound pressure level (turbulent) (Lp_Ae_1m_turb): {Lp_Ae_1m_turb:.2f} dB")
        print(f"Overall transmission loss (cavitation) (Tlcav): {Tlcav:.2f} dB")
        print(f"External sound pressure level (cavitating) (Lp_Ae_1m_cav): {Lp_Ae_1m_cav:.2f} dB")
        print(f"Frequency distribution function (turbulent) (Fturb_fi): {Fturb_fi:.2f} dB")
        print(f"Frequency distribution function (cavitating) (Fcav_fi): {Fcav_fi:.2f} dB")
        print(f"Internal sound pressure level at fi (turbulent) (Lpi_fi_turb): {Lpi_fi_turb:.2f} dB")
        print(f"Internal sound pressure level at fi (cavitating) (Lpi_fi_cav): {Lpi_fi_cav:.2f} dB")
        print(f"Transmission loss corrected for fi (delta_TLf_fi): {delta_TLf_fi:.2f} dB")
        print(f"Transmission loss at fi (TL_fi): {TL_fi:.2f} dB")
        print(f"External sound pressure level at fi (Lpe_1m_fi): {Lpe_1m_fi:.2f} dB")


if __name__ == "__main__":
    calculate_valve_flow()
