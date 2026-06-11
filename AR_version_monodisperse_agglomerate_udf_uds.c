/*============================================================================*/
/*                                                                            */
/*   Soot Formation UDF for ANSYS Fluent                                      */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*   Author      : Arash Fakharnezhad                                         */
/*   Affiliation : The University of Melbourne                                */
/*   Year        : 2026                                                       */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*   Description:                                                             */
/*   This User-Defined Function (UDF) implements a physics-based soot         */
/*   formation model for combustion simulations. The model accounts for:      */
/*       - Nucleation                                                         */
/*       - Coagulation                                                        */
/*       - Surface growth                                                     */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*   Methodology:                                                             */
/*   - Transport equations are solved using User-Defined Scalars (UDS):       */
/*       UDS-0 : Soot number density N[#/kg]                                  */
/*       UDS-1 : Soot carbon concentration C[kmol/kg]                         */
/*                                                                            */
/*----------------------------------------------------------------------------*/
/*   References:                                                              */
/*   - Appel et al. (2000). Combustion and Flame                              */
/*   - Brookes et al. (1999). Combustion and Flame                            */
/*   - Fakharnezhad et al. (2026). Powder Technology                          */
/*   - Fakharnezhad et al. (2025). Aerosol Research                           */
/*                                                                            */
/*============================================================================*/

#include "udf.h"
#include "mem.h"
#include "prop.h"
#include "sg_mem.h"
#include "materials.h"



/* ================= CONSTANTS ================= */

/* Mathematical constant */
real pi = 3.14159265359;

/* Physical constants */
real N_a = 6.022137e+26;       /* Avogadro's number, 1/kmol */
real R   = 8314.46261815324;   /* Gas constant, m3·Pa/(K·kmol) */
real KB  = 1.3806503e-23;      /* Boltzmann constant, m2·kg/s2·K */

/* Soot properties */
real rho_s = 1500.0;           /* Soot bulk density, kg/m3 */
real MW_C  = 12.0;             /* Molecular weight of C atom, kg/kmol */
real Xsoot = 2.3e19;           /* Surface radical density, 1/m2 */

/* Nucleation parameters */

/* Moss–Brookes nucleation model */
real C_a  = 54.0;        /* Nucleation pre-factor (1/s) */
real T_a  = 21100.0;     /* Activation temperature (K) */
real NC_P = 36;          /* Carbon atoms per nucleus (-) */

/* MD-based nucleation model */
real k_nucl_MD = 1.299e+32;   /* Pre-factor (m^7.5 / (kmol^3.5·s)) */
real n_nucl_MD = 3.5;         /* Reaction order (-) */



/*   Coagulation / Collision Kernel Parameters                                */

/* Enhancement factors (free molecular regime) */
real gamma_poly = 1.35;   /* Polydispersity enhancement factor (-) */
real gamma_vdW  = 2.2;    /* van der Waals enhancement factor (-) */

/* Species molecular weights */
real M_F  = 26.03788;          /* C2H2, kg/kmol */
real M_H2 = 2.016;             /* H2, kg/kmol */
real M_H  = 1.00784;           /* H, kg/kmol */

/* Reference conditions */
real p_ref = 101325.0;         /* Reference pressure, Pa */


/*----------------------------------------------------------------------------*/
/*   Nucleation Source Term for Soot Number Density, N[#/kg] (UDS-0)          */
/*----------------------------------------------------------------------------*/
/*   Available models:                                                        */
/*   - Moss–Brookes nucleation model (default)                                */
/*   - MD-based nucleation model (optional)                                   */

DEFINE_SOURCE(N_nucl, c, t, dS, eqn)
{
    real P, rho, T;
    real X_F;              /* Molar fraction of C2H2 */
    real conc_C2H2;        /* Concentration of C2H2 (kmol/m3) */
    real nuc_rate_N;       /* Nucleation rate (#/m3/s) */

    /* === Cell properties === */
    rho = C_R(c, t);              /* Gas density (kg/m3) */
    T   = C_T(c, t);              /* Temperature (K) */
    P   = C_P(c, t) + p_ref;      /* Pressure (Pa) */

	/* === Molar fraction of soot precursor (C2H2) === */
	/* C2H2 is the soot precursor (species index 17, must match mechanism) */
	/* Molar fraction computed from mass fraction (C_YI) using ideal gas relation */
    X_F = C_YI(c, t, 17) * rho * R * T / (P * M_F);

    /* === Convert to concentration (kmol/m3) === */
    conc_C2H2 = X_F * P / (R * T);

    /* Only one nucleation model should be active at a time */

    /* ========================================================= */
    /* === Moss–Brookes nucleation model (default) ============= */
    /* ========================================================= */
    
    nuc_rate_N = C_a * N_a * conc_C2H2 * exp(-T_a / T);

    /* ========================================================= */
    /* === MD-based nucleation model (uncomment to use) ======== */
    /* ========================================================= */
    /*
    nuc_rate_N = k_nucl_MD * pow(conc_C2H2, n_nucl_MD);
    */

    /* === Store for post-processing === */
    C_UDMI(c, t, 0) = nuc_rate_N;

    dS[eqn] = 0.0;

    return nuc_rate_N;
}



/*----------------------------------------------------------------------------*/
/*   Nucleation Source Term for Soot Carbon Concentration, C[kmol/kg] (UDS-1) */
/*----------------------------------------------------------------------------*/
/*   Available models:                                                        */
/*   - Moss–Brookes nucleation model (default)                                */
/*   - MD-based nucleation model (optional)                                   */

DEFINE_SOURCE(C_nucl, c, t, dS, eqn)
{
    real P, rho, T;
    real X_F;              /* Molar fraction of C2H2 */
    real conc_C2H2;        /* Concentration of C2H2 (kmol/m3) */
    real nuc_rate_N;       /* Nucleation rate (#/m3/s) */
    real nuc_rate_C;       /* Carbon nucleation rate (kmol/m3/s) */

    /* === Cell properties === */
    rho = C_R(c, t);              /* Gas density (kg/m3) */
    T   = C_T(c, t);              /* Temperature (K) */
    P   = C_P(c, t) + p_ref;      /* Pressure (Pa) */

	/* === Molar fraction of soot precursor (C2H2) === */
	/* C2H2 is the soot precursor (species index 17, must match mechanism) */
	/* Molar fraction computed from mass fraction (C_YI) using ideal gas relation */
    X_F = C_YI(c, t, 17) * rho * R * T / (P * M_F);

    /* === Convert to concentration (kmol/m3) === */
    conc_C2H2 = X_F * P / (R * T);

    /* Only one nucleation model should be active at a time */

    /* ========================================================= */
    /* === Moss–Brookes nucleation model (default) ============= */
    /* ========================================================= */
    
    nuc_rate_N = C_a * N_a * conc_C2H2 * exp(-T_a / T);

    /* ========================================================= */
    /* === MD-based nucleation model (uncomment to use) ======== */
    /* ========================================================= */
    /*
    nuc_rate_N = k_nucl_MD * pow(conc_C2H2, n_nucl_MD);
    */

    /* === Convert number-based nucleation to carbon formation === */
    
    nuc_rate_C = nuc_rate_N * NC_P / N_a;

    /* === Store for post-processing === */
    C_UDMI(c, t, 1) = nuc_rate_C;

    dS[eqn] = 0.0;

    return nuc_rate_C;
}



	
/*----------------------------------------------------------------------------*/
/*   Coagulation Source Term for Soot Number Density (UDS-0)                  */
/*----------------------------------------------------------------------------*/
/*   Description:                                                             */
/*   Computes soot number density change due to particle coagulation          */
/*   using a fractal aggregate model and collision kernel formulation.        */

DEFINE_SOURCE(N_coag, c, t, dS, eqn)
{
    real P, rho, T;
    real N_t, C_t;              /* Soot number density and carbon concentration */
    real coag_rate_N;           /* Coagulation rate (#/m3/s) */

    /* Aggregate properties */
    real V_t = 0.0, v_Agg = 0.0;
    real dp = 0.0, np_Agg = 0.0, m_Agg = 0.0;
    real dm = 0.0, dmdg = 0.0, dg = 0.0, dc = 0.0;

    /* Geometric parameters */
    real vol_equi_diam;
    real a_sphere_ref, v_sphere_ref, a_Agg, A_Ds;
    real soot_nucl_dia;

    /* Fractal model parameters */
    real Ds, y0, A1, t1, x;

    /* Collision kernel */
    real beta = 0.0;

    /* === Cell properties === */
    rho = C_R(c, t);             /* Gas density (kg/m3) */
    T   = C_T(c, t);             /* Temperature (K) */
    P   = C_P(c, t) + p_ref;     /* Pressure (Pa) */

    /* === Transported soot variables === */
    N_t = C_UDSI(c, t, 0);
    C_t = C_UDSI(c, t, 1);

    /* === Prevent negative values === */
    if (N_t < 0.0) N_t = 0.0;
    if (C_t < 0.0) C_t = 0.0;

    /* === Soot nucleus diameter === */
    soot_nucl_dia = pow(NC_P * MW_C * 6.0 / (rho_s * N_a * pi), 1.0 / 3.0);

    /* === Coagulation only if soot exists === */
    if (N_t > 0.0 && C_t > 0.0)
    {
        /* === Aggregate volume === */
        V_t  = C_t * MW_C / rho_s;   /* Total soot volume */
        v_Agg = V_t / N_t;           /* Volume per aggregate */

        /* === Equivalent diameter === */
        vol_equi_diam = pow((6.0 * v_Agg / pi), 1.0 / 3.0);

        a_sphere_ref = pi * pow(soot_nucl_dia, 2);
        v_sphere_ref = pi * pow(soot_nucl_dia, 3) / 6.0;

        /* === Fractal dimension (Ds) === */
        if (vol_equi_diam <= 6.67e-9)
        {
            Ds = 2.0;
        }
        else
        {
            y0 = 2.41143;
            A1 = -0.88858;
            t1 = 8.66525;
            x  = vol_equi_diam * 1e9;

            Ds = A1 * exp(-x / t1) + y0;
        }

        /* Limit Ds */
        if (Ds < 2.0) Ds = 2.0;
        if (Ds > 3.0) Ds = 3.0;

        /* === Aggregate geometry === */
        a_Agg = a_sphere_ref * pow((v_Agg / v_sphere_ref), (Ds / 3.0));
        A_Ds  = a_Agg * N_t;

        dp = 6.0 * V_t / A_Ds;

        /* Number of primary particles per aggregate */
        np_Agg = 6.0 * V_t / (N_t * pi * pow(dp, 3));

        /* === Aggregate mass === */
        m_Agg = rho_s * np_Agg * pi * pow(dp, 3) / 6.0;

        /* === Collision diameters === */
        dm = dp * pow(np_Agg, 0.45);   /* Mobility diameter */

        if (np_Agg > 1.8)
            dmdg = pow(np_Agg, -0.2) + 0.4;
        else
            dmdg = sqrt(5.0 / 3.0);

        dg = dm / dmdg;   /* Gyration diameter */
        dc = dg;          /* Collision diameter */

        /* === Collision kernel === */
        beta = gamma_poly * gamma_vdW * 4.0 * sqrt(pi * KB * T / m_Agg) * pow(dc, 2);

        /* === Coagulation rate === */

        coag_rate_N = -0.5 * beta * pow(rho * N_t, 2.0);

        dS[eqn] = -beta * pow(rho, 2.0) * N_t;  /* Linearized source term */
    }
    else
    {
        coag_rate_N = 0.0;
        dS[eqn] = 0.0;
    }

    /* === Store for post-processing === */
    C_UDMI(c, t, 2) = -coag_rate_N;     /* Coagulation magnitude */
    C_UDMI(c, t, 3) = rho * N_t;        /* Number density (#/m3) */
    C_UDMI(c, t, 4) = np_Agg;           /* Number of primary particles */
    C_UDMI(c, t, 5) = dp;               /* Primary particle diameter */
    C_UDMI(c, t, 6) = dm;               /* Mobility diameter */

    return coag_rate_N;
}


/*----------------------------------------------------------------------------*/
/*   Surface Growth Source Term for Soot Carbon Concentration (UDS-1)         */
/*----------------------------------------------------------------------------*/
/*   Description:                                                             */
/*   Computes soot mass growth via surface reactions using the HACA model     */
/*   (Hydrogen-Abstraction/Carbon-Addition mechanism).                        */

DEFINE_SOURCE(C_SG, c, t, dS, eqn)
{
    real SG_rate_C;   /* Surface growth rate (kmol/m3/s) */

    real P, rho, T;
    real X_F, X_H2, X_H;    /* Molar fractions */

    real N_t, C_t;

    /* Aggregate properties */
    real V_t = 0.0, v_Agg = 0.0;
    real dp = 0.0, np_Agg = 0.0, m_Agg = 0.0;

    /* Geometry */
    real vol_equi_diam;
    real a_sphere_ref, v_sphere_ref, soot_nucl_dia;
    real a_Agg, A_Ds;

    /* Fractal model */
    real Ds, y0, A1, t1, x;

    /* HACA parameters */
    real Ks, alpha, H_H2_ratio, Corrected_Xsoot;

    /* === Cell properties === */
    rho = C_R(c, t);          /* Gas density (kg/m3) */
    T   = C_T(c, t);          /* Temperature (K) */
    P   = C_P(c, t) + p_ref;  /* Pressure (Pa) */

    /* === Transported soot variables === */
    N_t = C_UDSI(c, t, 0);
    C_t = C_UDSI(c, t, 1);

    /* === Prevent negative values === */
    if (N_t < 0.0) N_t = 0.0;
    if (C_t < 0.0) C_t = 0.0;

    /* === Species molar fractions === */
    /* Species indices must match the reaction mechanism */

    /* C2H2 (soot precursor) */
    X_F  = C_YI(c, t, 17) * rho * R * T / (P * M_F);

    /* H2 */
    X_H2 = C_YI(c, t, 0) * rho * R * T / (P * M_H2);

    /* H */
    X_H  = C_YI(c, t, 1) * rho * R * T / (P * M_H);

    /* === Soot nucleus diameter === */
    soot_nucl_dia = pow(NC_P * MW_C * 6.0 / (rho_s * N_a * pi), 1.0 / 3.0);

    /* === Surface growth only if soot exists === */
    if (N_t > 0.0 && C_t > 0.0 && X_H2 > 0.0)
    {
        /* === Aggregate volume === */
        V_t  = C_t * MW_C / rho_s;
        v_Agg = V_t / N_t;

        /* === Equivalent diameter === */
        vol_equi_diam = pow((6.0 * v_Agg / pi), 1.0 / 3.0);

        a_sphere_ref = pi * pow(soot_nucl_dia, 2);
        v_sphere_ref = pi * pow(soot_nucl_dia, 3) / 6.0;

        /* === Fractal dimension (Ds) === */
        if (vol_equi_diam <= 6.67e-9)
        {
            Ds = 2.0;
        }
        else
        {
            y0 = 2.41143;
            A1 = -0.88858;
            t1 = 8.66525;
            x  = vol_equi_diam * 1e9;

            Ds = A1 * exp(-x / t1) + y0;
        }

        if (Ds < 2.0) Ds = 2.0;
        if (Ds > 3.0) Ds = 3.0;

        /* === Aggregate geometry === */
        a_Agg = a_sphere_ref * pow((v_Agg / v_sphere_ref), (Ds / 3.0));
        A_Ds  = a_Agg * N_t;

        dp = 6.0 * V_t / A_Ds;
        np_Agg = 6.0 * V_t / (N_t * pi * pow(dp, 3));

        /* === Aggregate mass === */
        m_Agg = rho_s * np_Agg * pi * pow(dp, 3) / 6.0;

        /* === HACA kinetic parameters === */
        Ks = 80.0e3 * pow(T, 1.56) * exp(-1912.4 / T);

        alpha = tanh(((12.65 - 0.00563 * T) /
                     log10((rho_s * pi * pow(dp, 3) * N_a / 6.0) / MW_C))
                     - 1.38 + 0.00068 * T);

        /* === H/H2 ratio correction === */
        H_H2_ratio = X_H / X_H2;
        
        Corrected_Xsoot = H_H2_ratio * 0.685 * Xsoot;

        /* Only one surface growth model should be active */

        /*=========================================================*/
        /*   HACA Surface Growth Model (default)                   */
        /*=========================================================*/
        /* C_dot = 2 * Ks * alpha * X_soot * A * [C2H2] * N * rho / N_A */
        SG_rate_C = 2.0 * Ks * alpha * Corrected_Xsoot *
                    a_Agg * (X_F * P / (R * T)) * N_t * rho / N_a;

        /*=========================================================*/
        /*   MD-based Surface Growth Model (optional)              */
        /*=========================================================*/
        /*
        SG_rate_C = exp(-2017.3/T + 3.946) * H_H2_ratio * 0.685 *
                    a_Agg * (X_F * P / (R * T)) * N_t * rho;
        */
    }

    /* === Store results === */
    C_UDMI(c, t, 7) = SG_rate_C;

    /* Soot volume fraction */
    C_UDMI(c, t, 8) = C_t * MW_C / rho_s * rho;

    /* === Source term === */
    dS[eqn] = 0.0;

    return SG_rate_C;
}



/*----------------------------------------------------------------------------*/
/*   Scalar Diffusivity for User-Defined Scalars (UDS)                        */
/*----------------------------------------------------------------------------*/
/*   Description:                                                             */
/*   Defines scalar diffusivity based on laminar viscosity and a constant     */
/*   Schmidt number assumption.                                               */

DEFINE_DIFFUSIVITY(scalar_diff, c, t, i)
{
    real D;   /* Scalar diffusivity (m^2/s) */

    /* === Diffusivity model === */
    /* D = mu / Sc, with Sc = 0.7 */
    D = C_MU_L(c, t) / 0.7;

    /* === Store for post-processing === */
    C_UDMI(c, t, 9) = D;              /* Scalar diffusivity (m^2/s) */

    return D;
}
