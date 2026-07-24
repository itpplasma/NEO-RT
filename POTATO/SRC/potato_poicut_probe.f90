program potato_poicut_probe
    use potato_input_mod, only: read_potato_input, rho_pol, rho_pol_max, &
        npoicut, edge_extension
    use field_eq_mod, only: allow_sol
    use poicut_mod, only: npc, rpc_arr, zpc_arr, rmagaxis, zmagaxis, &
        Rbou_lfs, Zbou_lfs, Rbou_hfs, Zbou_hfs

    implicit none

    call read_potato_input('potato.in')
    allow_sol = edge_extension

    call find_poicut(rho_pol_max, npoicut)
    call rhopol_boundary(rho_pol)

    write(*, '(a,1x,*(es24.16e3,1x))') 'POTATO_POICUT', &
        rpc_arr(0), zpc_arr(0), rpc_arr(npc), zpc_arr(npc), &
        rmagaxis, zmagaxis, Rbou_hfs, Zbou_hfs, Rbou_lfs, Zbou_lfs

end program potato_poicut_probe
