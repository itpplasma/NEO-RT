def get_torque(data):
    NM_IN_DYNCM = 1e-7
    CM3_IN_M3 = 1e-6

    dVds = data[1]
    Tco = data[3]
    Tctr = data[4]
    Tt = data[5]
    Ttot = Tco + Tctr + Tt  # Torque density in CGS over s
    Tphi_dens = (Ttot*NM_IN_DYNCM)/(dVds*CM3_IN_M3)  # Torque density in SI
    return Tphi_dens
