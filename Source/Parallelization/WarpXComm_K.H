#ifndef WARPX_COMM_K_H_
#define WARPX_COMM_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_x (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Bxa,
                            amrex::Array4<amrex::Real const> const& Bxf,
                            amrex::Array4<amrex::Real const> const& Bxc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;
    Bxa(j,k,l) = owx * Bxc(jg,kg,lg) + wx * Bxc(jg+1,kg,lg) + Bxf(j,k,l);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_y (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Bya,
                            amrex::Array4<amrex::Real const> const& Byf,
                            amrex::Array4<amrex::Real const> const& Byc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    // Note that for 2d, l=0, because the amrex convention is used here.

#if (AMREX_SPACEDIM == 3)
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Bya(j,k,l) = owy * Byc(jg,kg,lg) + wy * Byc(jg,kg+1,lg) + Byf(j,k,l);
#else
    Bya(j,k,l) = Byc(jg,kg,lg) + Byf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_z (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Bza,
                            amrex::Array4<amrex::Real const> const& Bzf,
                            amrex::Array4<amrex::Real const> const& Bzc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    // Note that for 2d, l=0, because the amrex convention is used here.

#if (AMREX_SPACEDIM == 3)
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;
    Bza(j,k,l) = owz * Bzc(jg,kg,lg) + owz * Bzc(jg,kg,lg+1) + Bzf(j,k,l);
#else
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Bza(j,k,l) = owy * Bzc(jg,kg,lg) + owy * Bzc(jg,kg+1,lg) + Bzf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_x (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Exa,
                            amrex::Array4<amrex::Real const> const& Exf,
                            amrex::Array4<amrex::Real const> const& Exc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 3)
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;
    Exa(j,k,l) = owy * owz * Exc(jg  ,kg  ,lg  )
        +         wy * owz * Exc(jg  ,kg+1,lg  )
        +        owy *  wz * Exc(jg  ,kg  ,lg+1)
        +         wy *  wz * Exc(jg  ,kg+1,lg+1)
        + Exf(j,k,l);
#else
    Exa(j,k,l) = owy * Exc(jg,kg,lg) + wy * Exc(jg,kg+1,lg) + Exf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_y (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Eya,
                            amrex::Array4<amrex::Real const> const& Eyf,
                            amrex::Array4<amrex::Real const> const& Eyc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

#if (AMREX_SPACEDIM == 3)
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;
    Eya(j,k,l) = owx * owz * Eyc(jg  ,kg  ,lg  )
        +         wx * owz * Eyc(jg+1,kg  ,lg  )
        +        owx *  wz * Eyc(jg  ,kg  ,lg+1)
        +         wx *  wz * Eyc(jg+1,kg  ,lg+1)
        + Eyf(j,k,l);
#else
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Eya(j,k,l) = owx * owy * Eyc(jg  ,kg  ,lg)
        +         wx * owy * Eyc(jg+1,kg  ,lg)
        +        owx *  wy * Eyc(jg  ,kg+1,lg)
        +         wx *  wy * Eyc(jg+1,kg+1,lg)
        + Eyf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_z (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Eza,
                            amrex::Array4<amrex::Real const> const& Ezf,
                            amrex::Array4<amrex::Real const> const& Ezc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

#if (AMREX_SPACEDIM == 3)
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Eza(j,k,l) = owx * owy * Ezc(jg  ,kg  ,lg  )
        +         wx * owy * Ezc(jg+1,kg  ,lg  )
        +        owx *  wy * Ezc(jg  ,kg+1,lg  )
        +         wx *  wy * Ezc(jg+1,kg+1,lg  )
        + Ezf(j,k,l);
#else
    Eza(j,k,l) = owx * Ezc(jg,kg,lg) + wx * Ezc(jg+1,kg,lg) + Ezf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf,
                               amrex::Array4<amrex::Real const> const& Bxc,
                               amrex::Array4<amrex::Real const> const& Bxg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Bxg(jg  ,kg  ,0)
        +     owx *  wy * Bxg(jg  ,kg+1,0)
        +      wx * owy * Bxg(jg+1,kg  ,0)
        +      wx *  wy * Bxg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wy = 0.5-wy;  owy = 1.0-wy;
    Real bc = owx * owy * Bxc(jg  ,kg  ,0)
        +     owx *  wy * Bxc(jg  ,kg-1,0)
        +      wx * owy * Bxc(jg+1,kg  ,0)
        +      wx *  wy * Bxc(jg+1,kg-1,0);

    // interp from fine staggered to fine nodal
    Real bf = 0.5*(Bxf(j,k-1,0) + Bxf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);

#endif

    Bxa(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf,
                               amrex::Array4<amrex::Real const> const& Byc,
                               amrex::Array4<amrex::Real const> const& Byg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Byg(jg  ,kg  ,0)
        +     owx *  wy * Byg(jg  ,kg+1,0)
        +      wx * owy * Byg(jg+1,kg  ,0)
        +      wx *  wy * Byg(jg+1,kg+1,0);

    // interp from coarse stagged (cell-centered for By) to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    wy = 0.5-wy;  owy = 1.0-wy;
    Real bc = owx * owy * Byc(jg  ,kg  ,0)
        +     owx *  wy * Byc(jg  ,kg-1,0)
        +      wx * owy * Byc(jg-1,kg  ,0)
        +      wx *  wy * Byc(jg-1,kg-1,0);

    // interp form fine stagger (cell-centered for By) to fine nodal
    int joff = (j == jg*2) ? -1 : 1;
    int koff = (k == kg*2) ? -1 : 1;
    Real bf = 0.25*(Byf(j,k,0) + Byf(j+joff,k,0) + Byf(j,k+koff,0) + Byf(j+joff,k+koff,0));

#else

#endif

    Bya(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf,
                               amrex::Array4<amrex::Real const> const& Bzc,
                               amrex::Array4<amrex::Real const> const& Bzg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Bzg(jg  ,kg  ,0)
        +     owx *  wy * Bzg(jg  ,kg+1,0)
        +      wx * owy * Bzg(jg+1,kg  ,0)
        +      wx *  wy * Bzg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    Real bc = owx * owy * Bzc(jg  ,kg  ,0)
        +     owx *  wy * Bzc(jg  ,kg+1,0)
        +      wx * owy * Bzc(jg-1,kg  ,0)
        +      wx *  wy * Bzc(jg-1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real bf = 0.5*(Bzf(j-1,k,0) + Bzf(j,k,0));

#else

#endif

    Bza(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf)
{
#if (AMREX_SPACEDIM == 2)
    Bxa(j,k,0) = 0.5*(Bxf(j,k-1,0) + Bxf(j,k,0));
#else
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf)
{
#if (AMREX_SPACEDIM == 2)
    int jg = amrex::coarsen(j,2);
    int kg = amrex::coarsen(k,2);
    int joff = (j == jg*2) ? -1 : 1;
    int koff = (k == kg*2) ? -1 : 1;
    Bya(j,k,0) = 0.25*(Byf(j,k,0) + Byf(j+joff,k,0) + Byf(j,k+koff,0) + Byf(j+joff,k+koff,0));
#else
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf)
{
#if (AMREX_SPACEDIM == 2)
    Bza(j,k,0) = 0.5*(Bzf(j-1,k,0) + Bzf(j,k,0));
#else
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf,
                               amrex::Array4<amrex::Real const> const& Exc,
                               amrex::Array4<amrex::Real const> const& Exg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Exg(jg  ,kg  ,0)
        +     owx *  wy * Exg(jg  ,kg+1,0)
        +      wx * owy * Exg(jg+1,kg  ,0)
        +      wx *  wy * Exg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    Real ec = owx * owy * Exc(jg  ,kg  ,0)
        +     owx *  wy * Exc(jg  ,kg+1,0)
        +      wx * owy * Exc(jg-1,kg  ,0)
        +      wx *  wy * Exc(jg-1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Exf(j-1,k,0) + Exf(j,k,0));

#else

#endif

    Exa(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf,
                               amrex::Array4<amrex::Real const> const& Eyc,
                               amrex::Array4<amrex::Real const> const& Eyg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal and coarse staggered to fine nodal
    Real eg = owx * owy * (Eyg(jg  ,kg  ,0) + Eyc(jg  ,kg  ,0))
        +     owx *  wy * (Eyg(jg  ,kg+1,0) + Eyc(jg  ,kg+1,0))
        +      wx * owy * (Eyg(jg+1,kg  ,0) + Eyc(jg+1,kg  ,0))
        +      wx *  wy * (Eyg(jg+1,kg+1,0) + Eyc(jg+1,kg+1,0));
    Real ec = 0.0;

    // interp from fine staggered to fine nodal
    Real ef = Eyf(j,k,0);

#else

#endif

    Eya(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf,
                               amrex::Array4<amrex::Real const> const& Ezc,
                               amrex::Array4<amrex::Real const> const& Ezg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Ezg(jg  ,kg  ,0)
        +     owx *  wy * Ezg(jg  ,kg+1,0)
        +      wx * owy * Ezg(jg+1,kg  ,0)
        +      wx *  wy * Ezg(jg+1,kg+1,0);

    // interp from coarse stagged to fine nodal
    wy = 0.5-wy;  owy = 1.0-wy;
    Real ec = owx * owy * Ezc(jg  ,kg  ,0)
        +     owx *  wy * Ezc(jg  ,kg-1,0)
        +      wx * owy * Ezc(jg+1,kg  ,0)
        +      wx *  wy * Ezc(jg+1,kg-1,0);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Ezf(j,k-1,0) + Ezf(j,k,0));

#else

#endif

    Eza(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf)
{
#if (AMREX_SPACEDIM == 2)
    Exa(j,k,0) = 0.5*(Exf(j-1,k,0) + Exf(j,k,0));
#else

#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf)
{
#if (AMREX_SPACEDIM == 2)
    Eya(j,k,0) = Eyf(j,k,0);
#else

#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf)
{
#if (AMREX_SPACEDIM == 2)
    Eza(j,k,0) =  0.5*(Ezf(j,k-1,0) + Ezf(j,k,0));
#else

#endif
}

#endif
