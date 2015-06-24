# Ice shelf version of the Rutgers version of the Regional Oceanic Modeling System (ROMS)

### Developer: Benjamin K. Galton-Fenzi

Ice-shelf model OPTIONS:                                                 
                                                                         
- ICESHELF            use if including ice shelf cavities                  
- ICESHELF_MORPH      use if allow ice draft to evolve                     
- ICESHELF_2EQN_VBC   use to activate 2-equation ice/ocean thermodynamics            
- ICESHELF_3EQN_VBC   use to activate 3-equation ice/ocean thermodynamics  
- ICESHELF_TEOS10     use for teos10 in situ freezing point temperature    
- ANA_SEAICE          use to prescribe simple open ocean sea ice model     
- SEAICE_CLIMA        use to prescribe seasonal climatology surface fluxes 
- SEAICE_WINTER       use to prescribe constant winter surface fluxes      

Tide model OPTIONS:

- NODAL_TIDES         use for astronomical tides with nodal corrections
- ANA_TIDES           use if analytical astronomical tides                

Includes the following test cases for:

- ICETEST: following ISOMIP specifications
- ICESHELF3D_TOY: Small domain (100x200 km and 500 deep)
- ICESHELF2D
- ICESHELF_TIDES: ICESHELF3D_TOY with tidal forcing
- ICECLIFF2D_TOY: Small toy ice cliff test case

-
-
