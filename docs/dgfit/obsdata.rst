#############
Observed Data
#############

The observed data to be fit is given in ascii data files for each type of data.
Observed data types supported include:

* avnhi: Dust-to-Gas as A(V)/N(HI)
* ext: Extinction curves given as A(lambda)/A(V)
* ir_emis: Infrared Emission given in MJy/sr/(10^20 H atoms) units
* abund: Atomic Dust Abundances in # atoms/(10^6 H atoms) units
* scat_a: Dust scattering albedo
* scat_g: Dust scattering phase function asymmetry (g)

A master obsdata file specifies all the possible observed data for a fit.

An example for the Milky Way Average (e.g., R(V)=3.1) is:

::

    # type filename
    avnhi   data/MW_diffuse_Gordon09_avnhi.dat
    ext     data/MW_diffuse_Gordon23_ext.dat
    abund   data/MW_diffuse_Jenkins09_abundances.dat
    ir_emis data/MW_diffuse_Compiegne11_ir_emission.dat
    scat_a  data/MW_diffuse_DGL_Gordon04_albedo.dat
    scat_g  data/MW_diffuse_DGL_Gordon04_g.dat

Examples
========

All examples are for the Milky Way Diffuse Average (R(V) = 3.1).

avnhi
-----

This file is always needed to run DGFit.

::

    # from Gordon et al. (2009)
    # average for R(V) = 3.1 curves
    #
    # Av_to_NHI  unc
    5.7e-22   0.2e-22

ext
---

Only the first few lines of data are shown.

::

    # wave A(l)/A(V) unc
    0.09199999999999997 6.130568301067063 0.06130568301067063
    0.09474547805176656 5.719941867413974 0.05719941867413974
    0.09757288707888888 5.349876545525296 0.053498765455252964
    0.10048467207804745 5.016544683743152 0.05016544683743152
    0.10348335100988702 4.716475274654454 0.04716475274654454
    0.10657151697641855 4.446523390036882 0.04446523390036882
    0.10975184046339928 4.203842638027271 0.04203842638027271
    0.11302707164963025 3.9858604940924947 0.03985860494092495

abund
-----

::

    # atom  abund  abund_unc  total_abund  total_abund_unc
    C    83.41   8.341   288.4  28.8
    O   109.26  10.926   575.4  57.5
    Mg   31.90   3.190    41.7   4.2
    Si   31.24   3.124    40.7   4.1
    Fe   33.34   3.334    34.7   3.5
    # abundances from Jenkins (2009) for F_* = 0.36   gas/dust = 150
    #   appropriate for the diffuse ISM (esp. IR emission)
    #   see Gordon et al. (2014)
    #
    # Units are # atoms per 10^6 H atoms

ir_emis
-------

Only the first few lines of data are shown.

::

    #INSTRU  FILTER    WAVE     SPEC         ERROR        UNIT   
    FIRAS   SPECTRUM  1200.00  0.0130615    0.00313187   MJy/sr  
    FIRAS   SPECTRUM  1102.60  0.0191160    0.00629346   MJy/sr  
    FIRAS   SPECTRUM  1050.10  0.0169294    0.0135306    MJy/sr  
    FIRAS   SPECTRUM  1002.37  0.0130511    0.00842313   MJy/sr  
    FIRAS   SPECTRUM  958.786  0.0235355    0.00740703   MJy/sr  
    FIRAS   SPECTRUM  918.836  0.0376609    0.0143819    MJy/sr  
    FIRAS   SPECTRUM  882.083  0.0352799    0.0149167    MJy/sr  
    FIRAS   SPECTRUM  848.157  0.0284347    0.00481095   MJy/sr  
    FIRAS   SPECTRUM  816.744  0.0309234    0.00351729   MJy/sr  

scat_a
------

Only the first few lines of data are shown.

::

    # wave albedo unc ref
    0.09000000000000001 0.28 0.04 "Sujatha et al. 2007 (DGL/FUSE/Voyager)"
    0.1 0.28 0.04 "Sujatha et al. 2007 (DGL/FUSE/Voyager)"
    0.10500000000000001 0.55 0.25 "Murthy et al 1993 (DGL/Voyager)"
    0.11 0.28 0.04 "Sujatha et al. 2007 (DGL/FUSE/Voyager)"
    0.11 0.4 0.1 "Sujatha et al. 2005 (DGL/Voyager)"
    0.125 0.4 0.2 "Shalima & Murthy 2004 (DGL/Voyager)"
    0.15 0.57 0.15 "Petersohn 1997 (DGL/DE)"


scat_g
------

Only the first few lines of data are shown.

::

    0.09000000000000001 0.61 0.07 "Sujatha et al. 2007 (DGL/FUSE/Voyager)"
    0.1 0.61 0.07 "Sujatha et al. 2007 (DGL/FUSE/Voyager)"
    0.11 0.61 0.07 "Sujatha et al. 2007 (DGL/FUSE/Voyager)"
    0.11 0.55 0.25 "Sujatha et al. 2005 (DGL/Voyager)"
    0.15 0.89 0.05 "Petersohn 1997 (DGL/DE)"
    0.155 0.7 0.1 "Lillie & Witt 1976 (DGL/OAO-2)"
    0.155 0.51 0.19 "Sujatha et al. 2010 (DGL/GALEX)"
