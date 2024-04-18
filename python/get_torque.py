def get_torque(data):
    NM_IN_DYNCM = 1e-7
    M3_IN_CM3 = 1e-6

    s = data[0]
    dVds = data[1]
    M_t = data[2]
    Tco = data[3]
    Tctr = data[4]
    Tt = data[5]
    Ttot = Tco + Tctr + Tt  # Torque density in CGS over s
    Tphi_dens = Ttot*NM_IN_DYNCM*dVds*M3_IN_CM3  # Torque density in SI
    return Tphi_dens
