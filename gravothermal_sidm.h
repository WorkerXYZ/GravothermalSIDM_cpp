#ifndef GRAVOTHERMAL_SIDM_H
#define GRAVOTHERMAL_SIDM_H

#include <stddef.h>

#define GSIDM_PATH_MAX 4096
#define GSIDM_NAME_MAX 256
#define GSIDM_MODEL_NAME_MAX 128

#ifndef GSIDM_OPENMP_MIN_SHELLS
#define GSIDM_OPENMP_MIN_SHELLS 2048
#endif

typedef enum {
    GSIDM_PROFILE_NFW = 0,
    GSIDM_PROFILE_HERNQUIST = 1,
    GSIDM_PROFILE_ISOTHERMAL = 2,
    GSIDM_PROFILE_TNFW = 3,
    GSIDM_PROFILE_THERN = 4
} gsidm_profile_t;

typedef enum {
    GSIDM_SCATTER_CONSTANT = 0,
    GSIDM_SCATTER_POWERLAW = 1,
    GSIDM_SCATTER_YUKAWA = 2
} gsidm_scatter_kind_t;

typedef struct {
    gsidm_scatter_kind_t kind;      /* scattering model type */
    double powerlaw_index;          /* exponent n in f(v)=v^n */
    int yukawa_index;               /* Yukawa velocity moment index */
    int yukawa_order;               /* approximation order for Yukawa model */
    int t_channel_only;             /* use t-channel-only Yukawa approximation */
    char raw_name[GSIDM_MODEL_NAME_MAX]; /* original model name */
} gsidm_scattering_model_t;

typedef struct {
    char run_name[GSIDM_NAME_MAX];
    char output_dir[GSIDM_PATH_MAX];

    gsidm_profile_t profile;        /* initial density profile */
    double r_s;                     /* scale radius, kpc */
    double rho_s;                   /* scale density, Msun/pc^3 */
    double r_cut;                   /* exponential cutoff radius in r_s units */
    double a;                       /* LMFP conductivity coefficient */
    double b;                       /* SMFP conductivity coefficient */
    double C;                       /* heat conduction calibration factor */
    double gamma;                   /* adiabatic index */
    gsidm_scattering_model_t lmfp_model; /* velocity dependence in LMFP regime */
    gsidm_scattering_model_t smfp_model; /* velocity dependence in SMFP regime */
    double sigma_m_with_units;      /* cross section per unit mass, cm^2/g */
    double w_units;                 /* velocity scale in scattering model, km/s */
    int n_shells;                   /* number of spherical mass shells */
    double r_min;                   /* innermost shell radius in r_s units */
    double r_max;                   /* outer halo radius in r_s units */
    double p_rmax_factor;           /* pressure integral upper limit, times r_max */
    int n_adjustment_fixed;         /* fixed hydrostatic adjustment count */
    int n_adjustment_max;           /* max hydrostatic adjustment count */
    int flag_timestep_use_relaxation; /* include relaxation-time step criterion */
    int flag_timestep_use_energy;   /* include energy-change step criterion */
    int flag_hydrostatic_initial;   /* enforce hydrostatic equilibrium at t=0 */

    double t_end;                   /* final dimensionless time */
    double rho_factor_end;          /* stop when central rho reaches this factor */
    double Kn_end;                  /* stop when minimum Kn falls below this value */
    int has_save_frequency_rate;    /* enable saving by t/t_relax cadence */
    double save_frequency_rate;     /* snapshots per relaxation time */
    int has_save_frequency_timing;  /* enable saving by fractional time change */
    double save_frequency_timing;   /* fractional time change per snapshot */
    int has_save_frequency_density; /* enable saving by central density change */
    double save_frequency_density;  /* fractional central density change per snapshot */
    double t_epsilon;               /* heat conduction timestep safety factor */
    double r_epsilon;               /* hydrostatic adjustment tolerance */

    int use_openmp;
    int omp_threads;
    int resume_enabled;
    char resume_target[GSIDM_PATH_MAX];
} gsidm_config_t;

void gsidm_set_default_config(gsidm_config_t *config);
int gsidm_parse_args(gsidm_config_t *config, int argc, char **argv);
int gsidm_run(const gsidm_config_t *config);
void gsidm_print_help(const char *program_name);

#endif
