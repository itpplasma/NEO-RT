# Plot canonical_freqs_vs_eta_t.dat (trapped)
set terminal wxt 0 title "Trapped particles"
set xlabel "eta [1/G]"
set ylabel "Frequency [rad/s]"
plot "canonical_freqs_vs_eta_t.dat" using 1:2 with lines title "omega_b", \
     "canonical_freqs_vs_eta_t.dat" using 1:3 with lines title "Omega_tor"

# Plot canonical_freqs_vs_eta_pco.dat (co-passing, sign_vpar=+1)
set terminal wxt 1 title "Co-passing (sign_vpar=+1)"
set xlabel "eta [1/G]"
set ylabel "Frequency [rad/s]"
plot "canonical_freqs_vs_eta_pco.dat" using 1:2 with lines title "omega_b", \
     "canonical_freqs_vs_eta_pco.dat" using 1:3 with lines title "Omega_tor"

# Plot canonical_freqs_vs_eta_pct.dat (co-passing, sign_vpar=-1)
set terminal wxt 2 title "Co-passing (sign_vpar=-1)"
set xlabel "eta [1/G]"
set ylabel "Frequency [rad/s]"
plot "canonical_freqs_vs_eta_pct.dat" using 1:2 with lines title "-omega_b", \
     "canonical_freqs_vs_eta_pct.dat" using 1:3 with lines title "Omega_tor"

pause -1 "Press Enter to exit"
