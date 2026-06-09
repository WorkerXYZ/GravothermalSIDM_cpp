#include "gravothermal_sidm.h"

#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define GSIDM_MKDIR(path) _mkdir(path)
#else
#include <sys/stat.h>
#include <sys/types.h>
#define GSIDM_MKDIR(path) mkdir((path), 0777)
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define GSIDM_PREFIX_DEFAULT "time"
#define GSIDM_PREFIX_DEBUG "broken"
#define GSIDM_CENTRAL_INDEX 3

#define GSIDM_MAGIC 0x475349444D42494EULL
#define GSIDM_VERSION 2ULL
#define GSIDM_VERSION_BASE 1ULL

typedef struct {
    char dir_data[GSIDM_PATH_MAX];
    char path_init[GSIDM_PATH_MAX];
    char path_index[GSIDM_PATH_MAX];
} gsidm_record_t;

/* quantities saved in snapshot files */
typedef struct {
    double rho_center;              /* initial central density, dimensionless */
    int n_adjustment;               /* hydrostatic adjustments in current step */
    int n_conduction;               /* total heat conduction steps */
    size_t n_save;                  /* snapshot counter */
    double t_epsilon;               /* heat conduction timestep safety factor */
    double r_epsilon;               /* hydrostatic adjustment tolerance */
    double t;                       /* current time, dimensionless */
    double t_before;                /* time before last conduction step */
    const double *m;                /* enclosed shell mass, dimensionless */
    const double *r;                /* shell outer radius, dimensionless */
    const double *rho;              /* shell-averaged density, dimensionless */
    const double *p;                /* shell-averaged pressure, dimensionless */
    const double *Mcy;              /* projected cylindrical mass, dimensionless */
    size_t n_shells;                /* number of spherical mass shells */
} gsidm_snapshot_t;

typedef struct {
    const gsidm_config_t *config;
    gsidm_record_t record;

    double scale_r;                 /* radius unit, m */
    double scale_rho;               /* density unit, kg/m^3 */
    double scale_m;                 /* mass unit, kg */
    double scale_u;                 /* specific energy unit, m^2/s^2 */
    double scale_p;                 /* pressure unit, kg/(m*s^2) */
    double scale_v;                 /* velocity unit, m/s */
    double scale_t;                 /* time unit, s */
    double scale_L;                 /* luminosity unit, kg*m^2/s^3 */
    double scale_sigma_m;           /* cross section per mass unit, m^2/kg */

    /* dimensionless cross section and velocity scale */
    double sigma_m;                 /* cross section per mass, dimensionless */
    double w;                       /* velocity scale in scattering model, dimensionless */
    double t_relax;                 /* relaxation time, dimensionless */

    double rho_center;              /* initial central density, dimensionless */
    int n_adjustment;               /* hydrostatic adjustments in current step */
    int n_conduction;               /* total heat conduction steps */
    size_t n_save;                  /* snapshot counter */
    double t_epsilon;               /* heat conduction timestep safety factor */
    double r_epsilon;               /* hydrostatic adjustment tolerance */
    double t;                       /* current time, dimensionless */
    double t_before;                /* time before last conduction step */

    /* quantities saved to file */
    double *r;                      /* shell outer radius, dimensionless */
    double *m;                      /* enclosed shell mass, dimensionless */
    double *Mcy;                    /* projected cylindrical mass, dimensionless */
    double *rho;                    /* shell-averaged density, dimensionless */
    double *p;                      /* shell-averaged pressure, dimensionless */

    /* derived quantities, recovered by update_derived_parameters */
    double *u;                      /* specific energy, dimensionless */
    double *v;                      /* 1d velocity dispersion, dimensionless */
    double *L;                      /* luminosity through shell boundary, dimensionless */
    double *Kn;                     /* Knudsen number */
    double *Kinv_lmfp;              /* inverse LMFP conductivity, dimensionless */
    double *Kinv_smfp;              /* inverse SMFP conductivity, dimensionless */

    double *r_ext;                  /* radius array with central boundary */
    double *m_ext;                  /* mass array with central boundary */
    double *a_arr;                  /* lower diagonal of hydrostatic matrix */
    double *b_arr;                  /* main diagonal of hydrostatic matrix */
    double *c_arr;                  /* upper diagonal of hydrostatic matrix */
    double *d_arr;                  /* RHS of hydrostatic matrix equation */
    double *delta_r;                /* shell radius change per adjustment */
    double *delta_uc;               /* specific energy change by conduction */
    double *delta_rho;              /* density change per adjustment */
    double *delta_ph;               /* pressure change per adjustment */
    double *keff;                   /* effective thermal conductivity, dimensionless */
    double *tdma_ba;                /* TDMA work array for main diagonal */
    double *tdma_dc;                /* TDMA work array for RHS */
    double *r2_ext;                 /* squared extended radius */
    double *r3_ext;                 /* cubed extended radius */
    double *r2_times_dr;            /* r^2 * delta_r */
    double *r3_shell;               /* cubed shell radius */
} gsidm_state_t;

typedef struct {
    const gsidm_config_t *config;
} gsidm_pressure_context_t;

static const double GSIDM_G = 6.67430e-11;                      /* m^3 kg^-1 s^-2 */
static const double GSIDM_SOLAR_MASS_KG = 1.988409870698051e30; /* kg */
static const double GSIDM_PARSEC_METERS = 3.0856775814913673e16; /* m */
static const double GSIDM_KPC_METERS = 1.0e3 * 3.0856775814913673e16; /* m */
static const double GSIDM_CM2_PER_G_TO_M2_PER_KG = 0.1;         /* (cm^2/g)->(m^2/kg) */
static const double GSIDM_KM_PER_S_TO_M_PER_S = 1000.0;         /* (km/s)->(m/s) */

static void gsidm_update_derived_parameters(gsidm_state_t *state);
static int gsidm_parallel_enabled(const gsidm_config_t *config);
static double gsidm_adaptive_simpson(
    double (*fn)(double, void *),
    void *ctx,
    double a,
    double b,
    double tol);

static int gsidm_string_equals(const char *left, const char *right) {
    return strcmp(left, right) == 0;
}

static void gsidm_trim(char *text) {
    size_t len;
    size_t start = 0;
    if (text == NULL) {
        return;
    }
    len = strlen(text);
    while (start < len &&
           (text[start] == ' ' || text[start] == '\t' || text[start] == '\r' || text[start] == '\n')) {
        start += 1;
    }
    if (start > 0) {
        memmove(text, text + start, len - start + 1);
    }
    len = strlen(text);
    while (len > 0 &&
           (text[len - 1] == ' ' || text[len - 1] == '\t' || text[len - 1] == '\r' || text[len - 1] == '\n')) {
        text[len - 1] = '\0';
        len -= 1;
    }
}

static void gsidm_copy_string(char *dst, size_t dst_size, const char *src) {
    if (dst_size == 0) {
        return;
    }
    snprintf(dst, dst_size, "%s", src);
}

static const char *gsidm_profile_name(gsidm_profile_t profile) {
    switch (profile) {
        case GSIDM_PROFILE_NFW:
            return "NFW";
        case GSIDM_PROFILE_HERNQUIST:
            return "Hernquist";
        case GSIDM_PROFILE_ISOTHERMAL:
            return "Isothermal";
        case GSIDM_PROFILE_TNFW:
            return "tNFW";
        case GSIDM_PROFILE_THERN:
            return "tHern";
        default:
            return "Unknown";
    }
}

static int gsidm_parse_profile(const char *text, gsidm_profile_t *profile) {
    if (gsidm_string_equals(text, "NFW")) {
        *profile = GSIDM_PROFILE_NFW;
        return 0;
    }
    if (gsidm_string_equals(text, "Hernquist")) {
        *profile = GSIDM_PROFILE_HERNQUIST;
        return 0;
    }
    if (gsidm_string_equals(text, "Isothermal")) {
        *profile = GSIDM_PROFILE_ISOTHERMAL;
        return 0;
    }
    if (gsidm_string_equals(text, "tNFW") || gsidm_string_equals(text, "TNFW")) {
        *profile = GSIDM_PROFILE_TNFW;
        return 0;
    }
    if (gsidm_string_equals(text, "tHern") || gsidm_string_equals(text, "THERN")) {
        *profile = GSIDM_PROFILE_THERN;
        return 0;
    }
    return -1;
}

static double gsidm_yukawa_born_viscosity_order1(
    double s,
    double par0,
    double par1,
    double par2,
    double par3,
    int t_channel_only) {
    double s2 = s * s + 1.0e-4;
    double s4 = s * s * s * s + 1.0e-8;
    double sn = t_channel_only ? s2 : s4;
    double g = 1.5 * log(1.0 + pow(sn * par1, par2)) / (par0 * par2 * s4);
    return pow(1.0 + 1.0 / pow(g, par3), -1.0 / par3);
}

static int gsidm_parse_scattering_model(
    const char *text,
    gsidm_scattering_model_t *model,
    char *error,
    size_t error_size) {
    char buffer[GSIDM_MODEL_NAME_MAX];
    char *token1;
    char *token2;

    memset(model, 0, sizeof(*model));
    gsidm_copy_string(model->raw_name, sizeof(model->raw_name), text);

    if (gsidm_string_equals(text, "constant")) {
        model->kind = GSIDM_SCATTER_CONSTANT;
        return 0;
    }

    if (strncmp(text, "powerlaw_", 9) == 0) {
        model->kind = GSIDM_SCATTER_POWERLAW;
        model->powerlaw_index = strtod(text + 9, NULL);
        return 0;
    }

    if (strstr(text, "YukawaBornViscosityApprox") != NULL) {
        gsidm_copy_string(buffer, sizeof(buffer), text);
        token1 = strchr(buffer, '_');
        if (token1 == NULL) {
            snprintf(error, error_size, "bad Yukawa model: %s", text);
            return -1;
        }
        *token1 = '\0';
        token1 += 1;
        token2 = strchr(token1, '_');
        if (token2 == NULL) {
            snprintf(error, error_size, "bad Yukawa model: %s", text);
            return -1;
        }
        *token2 = '\0';
        token2 += 1;

        if (token1[0] != 'K' || strncmp(token2, "order", 5) != 0) {
            snprintf(error, error_size, "bad Yukawa model: %s", text);
            return -1;
        }

        model->kind = GSIDM_SCATTER_YUKAWA;
        model->yukawa_index = atoi(token1 + 1);
        model->yukawa_order = atoi(token2 + 5);
        model->t_channel_only = strstr(text, "Tchannel") != NULL;
        return 0;
    }

    snprintf(error, error_size, "unknown scattering model: %s", text);
    return -1;
}

static double gsidm_eval_scattering_model(
    const gsidm_scattering_model_t *model,
    double x) {
    static const double params_t[4][4] = {
        {8.0, 0.303941, 0.74, 0.68},
        {24.0, 0.826198, 0.82, 0.76},
        {48.0, 1.36217, 0.83, 0.79},
        {80.0, 1.90106, 0.84, 0.81}
    };
    static const int params_t_index[4] = {3, 5, 7, 9};
    static const double params_f[4][4] = {
        {8.0, 0.0339848, 0.37, 0.63},
        {24.0, 0.251115, 0.41, 0.71},
        {48.0, 0.682602, 0.42, 0.74},
        {80.0, 1.32953, 0.43, 0.76}
    };
    int i;

    if (model->kind == GSIDM_SCATTER_CONSTANT) {
        return 1.0;
    }
    if (model->kind == GSIDM_SCATTER_POWERLAW) {
        return pow(x, model->powerlaw_index);
    }
    if (model->kind == GSIDM_SCATTER_YUKAWA) {
        const double (*params)[4] = model->t_channel_only ? params_t : params_f;
        if (model->yukawa_order == 1) {
            for (i = 0; i < 4; ++i) {
                if (params_t_index[i] == model->yukawa_index) {
                    return gsidm_yukawa_born_viscosity_order1(
                        x,
                        params[i][0],
                        params[i][1],
                        params[i][2],
                        params[i][3],
                        model->t_channel_only);
                }
            }
        } else if (model->yukawa_order == 2 && model->yukawa_index == 5) {
            double K5 = gsidm_yukawa_born_viscosity_order1(
                x, params[1][0], params[1][1], params[1][2], params[1][3], model->t_channel_only);
            double K7 = gsidm_yukawa_born_viscosity_order1(
                x, params[2][0], params[2][1], params[2][2], params[2][3], model->t_channel_only);
            double K9 = gsidm_yukawa_born_viscosity_order1(
                x, params[3][0], params[3][1], params[3][2], params[3][3], model->t_channel_only);
            return (28.0 * K5 * K5 + 80.0 * K5 * K9 - 64.0 * K7 * K7) /
                   (77.0 * K5 - 112.0 * K7 + 80.0 * K9);
        }
    }

    fprintf(stderr, "Unsupported scattering model at runtime: %s\n", model->raw_name);
    exit(1);
}

void gsidm_set_default_config(gsidm_config_t *config) {
    char error[128];

    memset(config, 0, sizeof(*config));

    /* set default input parameters */
    gsidm_copy_string(config->run_name, sizeof(config->run_name), "8e6-c45-s50-200");
    gsidm_copy_string(config->output_dir, sizeof(config->output_dir), "LensingPerturber");

    /* quantities to initialize halo */
    config->profile = GSIDM_PROFILE_NFW;     /* initial density profile */
    config->r_s = 0.09265;                   /* scale radius, kpc */
    config->rho_s = 0.2808;                  /* scale density, Msun/pc^3 */
    config->r_cut = 60.0;                    /* cutoff radius in r_s units */

    config->a = 4.0 / sqrt(M_PI);            /* LMFP conductivity coefficient */
    config->b = 25.0 * sqrt(M_PI) / 32.0;    /* SMFP conductivity coefficient */
    config->C = 0.753;                       /* heat conduction calibration factor */
    config->gamma = 5.0 / 3.0;               /* adiabatic index */
    if (gsidm_parse_scattering_model("constant", &config->lmfp_model, error, sizeof(error)) != 0) {
        fprintf(stderr, "%s\n", error);
        exit(1);
    }
    if (gsidm_parse_scattering_model("constant", &config->smfp_model, error, sizeof(error)) != 0) {
        fprintf(stderr, "%s\n", error);
        exit(1);
    }
    config->sigma_m_with_units = 50.0;       /* cross section per mass, cm^2/g */
    config->w_units = 1.0;                   /* velocity scale, km/s */

    config->n_shells = 200;                  /* number of spherical mass shells */
    config->r_min = 0.005;                   /* innermost radius in r_s units */
    config->r_max = 60.0;                    /* outer radius in r_s units */
    config->p_rmax_factor = 10.0;            /* pressure integral upper limit factor */
    config->n_adjustment_fixed = -1;         /* fixed hydrostatic adjustment count */
    config->n_adjustment_max = -1;           /* max hydrostatic adjustment count */
    config->flag_timestep_use_relaxation = 1; /* use relaxation timestep criterion */
    config->flag_timestep_use_energy = 1;    /* use energy timestep criterion */
    config->flag_hydrostatic_initial = 0;    /* adjust initial condition if nonzero */

    /* runtime quantities */
    config->t_end = 2000.0;                  /* final dimensionless time */
    config->rho_factor_end = INFINITY;       /* central density stop factor */
    config->Kn_end = 0.0;                    /* minimum Knudsen-number stop value */
    config->has_save_frequency_rate = 0;     /* enable saves by t/t_relax cadence */
    config->save_frequency_rate = 0.0;       /* snapshots per relaxation time */
    config->has_save_frequency_timing = 0;   /* enable saves by fractional time change */
    config->save_frequency_timing = 0.0;     /* fractional time change per snapshot */
    config->has_save_frequency_density = 1;  /* enable saves by central density change */
    config->save_frequency_density = 0.1;    /* fractional central density change */
    config->t_epsilon = 1.0e-4;              /* heat conduction timestep factor */
    config->r_epsilon = 1.0e-14;             /* hydrostatic adjustment tolerance */

    config->use_openmp = 1;
    config->omp_threads = 0;
    config->resume_enabled = 0;
    config->resume_target[0] = '\0';
}

static int gsidm_parse_bool(const char *text, int *value) {
    if (gsidm_string_equals(text, "1") || gsidm_string_equals(text, "true")) {
        *value = 1;
        return 0;
    }
    if (gsidm_string_equals(text, "0") || gsidm_string_equals(text, "false")) {
        *value = 0;
        return 0;
    }
    return -1;
}

static int gsidm_parse_optional_double(
    const char *text,
    int *has_value,
    double *value) {
    if (gsidm_string_equals(text, "none")) {
        *has_value = 0;
        *value = 0.0;
        return 0;
    }
    if (gsidm_string_equals(text, "inf")) {
        *has_value = 1;
        *value = INFINITY;
        return 0;
    }
    *has_value = 1;
    *value = strtod(text, NULL);
    return 0;
}

static double gsidm_parse_double(const char *text) {
    if (gsidm_string_equals(text, "inf")) {
        return INFINITY;
    }
    return strtod(text, NULL);
}

void gsidm_print_help(const char *program_name) {
    printf("Usage: %s [options]\n", program_name);
    printf("  --run-name NAME\n");
    printf("  --output-dir DIR\n");
    printf("  --profile {NFW|Hernquist|Isothermal|tNFW|tHern}\n");
    printf("  --r-s VALUE\n");
    printf("  --rho-s VALUE\n");
    printf("  --r-cut VALUE      tNFW/tHern cutoff radius in units of r_s\n");
    printf("  --sigma-m-with-units VALUE\n");
    printf("  --w-units VALUE\n");
    printf("  --n-shells VALUE\n");
    printf("  --r-min VALUE\n");
    printf("  --r-max VALUE\n");
    printf("  --p-rmax-factor VALUE\n");
    printf("  --n-adjustment-fixed VALUE\n");
    printf("  --n-adjustment-max VALUE\n");
    printf("  --flag-timestep-use-relaxation {0|1}\n");
    printf("  --flag-timestep-use-energy {0|1}\n");
    printf("  --flag-hydrostatic-initial {0|1}\n");
    printf("  --lmfp-model NAME\n");
    printf("  --smfp-model NAME\n");
    printf("  --t-end VALUE\n");
    printf("  --rho-factor-end VALUE|inf\n");
    printf("  --kn-end VALUE\n");
    printf("  --save-frequency-rate VALUE|none\n");
    printf("  --save-frequency-timing VALUE|none\n");
    printf("  --save-frequency-density VALUE|none\n");
    printf("  --t-epsilon VALUE\n");
    printf("  --r-epsilon VALUE\n");
    printf("  --use-openmp {0|1}\n");
    printf("  --threads VALUE\n");
    printf("  --resume latest|SNAPSHOT_FILE\n");
    printf("  --help\n");
}

int gsidm_parse_args(gsidm_config_t *config, int argc, char **argv) {
    int i;
    char error[128];

    for (i = 1; i < argc; ++i) {
        const char *arg = argv[i];
        const char *value;

        if (gsidm_string_equals(arg, "--help")) {
            gsidm_print_help(argv[0]);
            return 1;
        }

        if (i + 1 >= argc) {
            fprintf(stderr, "Missing value for option %s\n", arg);
            return -1;
        }
        value = argv[++i];

        if (gsidm_string_equals(arg, "--run-name")) {
            gsidm_copy_string(config->run_name, sizeof(config->run_name), value);
        } else if (gsidm_string_equals(arg, "--output-dir")) {
            gsidm_copy_string(config->output_dir, sizeof(config->output_dir), value);
        } else if (gsidm_string_equals(arg, "--profile")) {
            if (gsidm_parse_profile(value, &config->profile) != 0) {
                fprintf(stderr, "Unknown profile: %s\n", value);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--r-s")) {
            config->r_s = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--rho-s")) {
            config->rho_s = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--r-cut")) {
            config->r_cut = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--sigma-m-with-units")) {
            config->sigma_m_with_units = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--w-units")) {
            config->w_units = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--n-shells")) {
            config->n_shells = atoi(value);
        } else if (gsidm_string_equals(arg, "--r-min")) {
            config->r_min = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--r-max")) {
            config->r_max = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--p-rmax-factor")) {
            config->p_rmax_factor = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--n-adjustment-fixed")) {
            config->n_adjustment_fixed = atoi(value);
        } else if (gsidm_string_equals(arg, "--n-adjustment-max")) {
            config->n_adjustment_max = atoi(value);
        } else if (gsidm_string_equals(arg, "--flag-timestep-use-relaxation")) {
            if (gsidm_parse_bool(value, &config->flag_timestep_use_relaxation) != 0) {
                fprintf(stderr, "Bad boolean for %s: %s\n", arg, value);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--flag-timestep-use-energy")) {
            if (gsidm_parse_bool(value, &config->flag_timestep_use_energy) != 0) {
                fprintf(stderr, "Bad boolean for %s: %s\n", arg, value);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--flag-hydrostatic-initial")) {
            if (gsidm_parse_bool(value, &config->flag_hydrostatic_initial) != 0) {
                fprintf(stderr, "Bad boolean for %s: %s\n", arg, value);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--lmfp-model")) {
            if (gsidm_parse_scattering_model(value, &config->lmfp_model, error, sizeof(error)) != 0) {
                fprintf(stderr, "%s\n", error);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--smfp-model")) {
            if (gsidm_parse_scattering_model(value, &config->smfp_model, error, sizeof(error)) != 0) {
                fprintf(stderr, "%s\n", error);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--t-end")) {
            config->t_end = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--rho-factor-end")) {
            config->rho_factor_end = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--kn-end")) {
            config->Kn_end = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--save-frequency-rate")) {
            gsidm_parse_optional_double(value, &config->has_save_frequency_rate, &config->save_frequency_rate);
        } else if (gsidm_string_equals(arg, "--save-frequency-timing")) {
            gsidm_parse_optional_double(value, &config->has_save_frequency_timing, &config->save_frequency_timing);
        } else if (gsidm_string_equals(arg, "--save-frequency-density")) {
            gsidm_parse_optional_double(value, &config->has_save_frequency_density, &config->save_frequency_density);
        } else if (gsidm_string_equals(arg, "--t-epsilon")) {
            config->t_epsilon = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--r-epsilon")) {
            config->r_epsilon = gsidm_parse_double(value);
        } else if (gsidm_string_equals(arg, "--use-openmp")) {
            if (gsidm_parse_bool(value, &config->use_openmp) != 0) {
                fprintf(stderr, "Bad boolean for %s: %s\n", arg, value);
                return -1;
            }
        } else if (gsidm_string_equals(arg, "--threads")) {
            config->omp_threads = atoi(value);
        } else if (gsidm_string_equals(arg, "--resume")) {
            config->resume_enabled = 1;
            gsidm_copy_string(config->resume_target, sizeof(config->resume_target), value);
        } else {
            fprintf(stderr, "Unknown option: %s\n", arg);
            return -1;
        }
    }

    return 0;
}

static int gsidm_mkdir_if_needed(const char *path) {
    if (GSIDM_MKDIR(path) == 0) {
        return 0;
    }
    if (errno == EEXIST) {
        return 0;
    }
    return -1;
}

static int gsidm_ensure_dir_recursive(const char *path) {
    char tmp[GSIDM_PATH_MAX];
    size_t len;
    size_t i;

    if (path == NULL || path[0] == '\0') {
        return -1;
    }

    gsidm_copy_string(tmp, sizeof(tmp), path);
    len = strlen(tmp);
    if (len == 0) {
        return -1;
    }

    for (i = 1; i < len; ++i) {
        if (tmp[i] == '/' || tmp[i] == '\\') {
            char saved = tmp[i];
            tmp[i] = '\0';
            if (strlen(tmp) > 0 && strchr(tmp, ':') != tmp + strlen(tmp) - 1) {
                if (gsidm_mkdir_if_needed(tmp) != 0) {
                    return -1;
                }
            }
            tmp[i] = saved;
        }
    }

    if (gsidm_mkdir_if_needed(tmp) != 0) {
        return -1;
    }
    return 0;
}

static void gsidm_join_path(
    char *dst,
    size_t dst_size,
    const char *left,
    const char *right) {
    size_t len = strlen(left);
    if (len > 0 && (left[len - 1] == '/' || left[len - 1] == '\\')) {
        snprintf(dst, dst_size, "%s%s", left, right);
    } else {
        snprintf(dst, dst_size, "%s/%s", left, right);
    }
}

static void gsidm_record_init(gsidm_record_t *record, const gsidm_config_t *config) {
    gsidm_join_path(record->dir_data, sizeof(record->dir_data), config->output_dir, config->run_name);
    gsidm_join_path(record->path_init, sizeof(record->path_init), record->dir_data, "halo_init.txt");
    gsidm_join_path(record->path_index, sizeof(record->path_index), record->dir_data, "snapshot_index.csv");
}

static int gsidm_resolve_resume_snapshot_path(
    const gsidm_record_t *record,
    const gsidm_config_t *config,
    char *path_snapshot,
    size_t path_snapshot_size) {
    if (!config->resume_enabled) {
        return -1;
    }

    if (gsidm_string_equals(config->resume_target, "latest")) {
        FILE *fp = fopen(record->path_index, "rb");
        char line[8192];
        char last_line[8192];
        int found = 0;
        if (fp == NULL) {
            fprintf(stderr, "Cannot open snapshot index for resume: %s\n", record->path_index);
            return -1;
        }
        last_line[0] = '\0';
        while (fgets(line, sizeof(line), fp) != NULL) {
            gsidm_trim(line);
            if (line[0] == '\0') {
                continue;
            }
            if (strncmp(line, "save_id,", 8) == 0) {
                continue;
            }
            gsidm_copy_string(last_line, sizeof(last_line), line);
            found = 1;
        }
        fclose(fp);

        if (!found) {
            fprintf(stderr, "No snapshot entries found in %s\n", record->path_index);
            return -1;
        }

        {
            char csv_copy[8192];
            char *token;
            char *filename = NULL;
            int field = 0;
            gsidm_copy_string(csv_copy, sizeof(csv_copy), last_line);
            token = strtok(csv_copy, ",");
            while (token != NULL) {
                field += 1;
                if (field == 8) {
                    filename = token;
                    break;
                }
                token = strtok(NULL, ",");
            }
            if (filename == NULL || filename[0] == '\0') {
                fprintf(stderr, "Failed to parse latest snapshot filename from %s\n", record->path_index);
                return -1;
            }
            gsidm_join_path(path_snapshot, path_snapshot_size, record->dir_data, filename);
            return 0;
        }
    }

    if (strchr(config->resume_target, '/') != NULL || strchr(config->resume_target, '\\') != NULL) {
        gsidm_copy_string(path_snapshot, path_snapshot_size, config->resume_target);
    } else {
        gsidm_join_path(path_snapshot, path_snapshot_size, record->dir_data, config->resume_target);
    }
    return 0;
}

static int gsidm_apply_saved_initialization(gsidm_config_t *config, const gsidm_record_t *record) {
    FILE *fp = fopen(record->path_init, "rb");
    char line[4096];
    char error[128];
    if (fp == NULL) {
        fprintf(stderr, "Cannot open initialization file for resume: %s\n", record->path_init);
        return -1;
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
        char *eq;
        char *key;
        char *value;
        gsidm_trim(line);
        if (line[0] == '\0' || line[0] == '#') {
            continue;
        }
        eq = strchr(line, '=');
        if (eq == NULL) {
            continue;
        }
        *eq = '\0';
        key = line;
        value = eq + 1;
        gsidm_trim(key);
        gsidm_trim(value);

        if (gsidm_string_equals(key, "profile")) {
            if (gsidm_parse_profile(value, &config->profile) != 0) {
                fclose(fp);
                fprintf(stderr, "Bad profile in %s: %s\n", record->path_init, value);
                return -1;
            }
        } else if (gsidm_string_equals(key, "r_s_kpc")) {
            config->r_s = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "rho_s_msun_per_pc3")) {
            config->rho_s = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "r_cut")) {
            config->r_cut = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "a")) {
            config->a = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "b")) {
            config->b = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "C")) {
            config->C = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "gamma")) {
            config->gamma = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "lmfp_model")) {
            if (gsidm_parse_scattering_model(value, &config->lmfp_model, error, sizeof(error)) != 0) {
                fclose(fp);
                fprintf(stderr, "%s\n", error);
                return -1;
            }
        } else if (gsidm_string_equals(key, "smfp_model")) {
            if (gsidm_parse_scattering_model(value, &config->smfp_model, error, sizeof(error)) != 0) {
                fclose(fp);
                fprintf(stderr, "%s\n", error);
                return -1;
            }
        } else if (gsidm_string_equals(key, "sigma_m_with_units_cm2_per_g")) {
            config->sigma_m_with_units = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "w_units_km_per_s")) {
            config->w_units = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "n_shells")) {
            config->n_shells = atoi(value);
        } else if (gsidm_string_equals(key, "r_min")) {
            config->r_min = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "r_max")) {
            config->r_max = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "p_rmax_factor")) {
            config->p_rmax_factor = gsidm_parse_double(value);
        } else if (gsidm_string_equals(key, "n_adjustment_fixed")) {
            config->n_adjustment_fixed = atoi(value);
        } else if (gsidm_string_equals(key, "n_adjustment_max")) {
            config->n_adjustment_max = atoi(value);
        } else if (gsidm_string_equals(key, "flag_timestep_use_relaxation")) {
            if (gsidm_parse_bool(value, &config->flag_timestep_use_relaxation) != 0) {
                fclose(fp);
                fprintf(stderr, "Bad boolean in %s: %s=%s\n", record->path_init, key, value);
                return -1;
            }
        } else if (gsidm_string_equals(key, "flag_timestep_use_energy")) {
            if (gsidm_parse_bool(value, &config->flag_timestep_use_energy) != 0) {
                fclose(fp);
                fprintf(stderr, "Bad boolean in %s: %s=%s\n", record->path_init, key, value);
                return -1;
            }
        } else if (gsidm_string_equals(key, "flag_hydrostatic_initial")) {
            if (gsidm_parse_bool(value, &config->flag_hydrostatic_initial) != 0) {
                fclose(fp);
                fprintf(stderr, "Bad boolean in %s: %s=%s\n", record->path_init, key, value);
                return -1;
            }
        }
    }

    fclose(fp);
    return 0;
}

static int gsidm_load_snapshot_into_state(gsidm_state_t *state, const char *path_snapshot) {
    FILE *fp = fopen(path_snapshot, "rb");
    uint64_t magic;
    uint64_t version;
    uint64_t shell_count;
    int64_t n_adjustment;
    int64_t n_conduction;
    uint64_t n_save;

    if (fp == NULL) {
        fprintf(stderr, "Cannot open snapshot for resume: %s\n", path_snapshot);
        return -1;
    }

    /* retrieve information from binary data */
    if (fread(&magic, sizeof(uint64_t), 1, fp) != 1 ||
        fread(&version, sizeof(uint64_t), 1, fp) != 1 ||
        fread(&shell_count, sizeof(uint64_t), 1, fp) != 1) {
        fclose(fp);
        fprintf(stderr, "Failed to read snapshot header: %s\n", path_snapshot);
        return -1;
    }
    if (magic != GSIDM_MAGIC || (version != GSIDM_VERSION_BASE && version != GSIDM_VERSION)) {
        fclose(fp);
        fprintf(stderr, "Snapshot format mismatch: %s\n", path_snapshot);
        return -1;
    }
    if ((size_t)shell_count != (size_t)state->config->n_shells) {
        fclose(fp);
        fprintf(stderr,
                "Snapshot shell count (%llu) does not match current configuration (%d): %s\n",
                (unsigned long long)shell_count,
                state->config->n_shells,
                path_snapshot);
        return -1;
    }

    if (fread(&state->rho_center, sizeof(double), 1, fp) != 1 ||
        fread(&n_adjustment, sizeof(int64_t), 1, fp) != 1 ||
        fread(&n_conduction, sizeof(int64_t), 1, fp) != 1 ||
        fread(&n_save, sizeof(uint64_t), 1, fp) != 1 ||
        fread(&state->t_epsilon, sizeof(double), 1, fp) != 1 ||
        fread(&state->r_epsilon, sizeof(double), 1, fp) != 1 ||
        fread(&state->t, sizeof(double), 1, fp) != 1 ||
        fread(&state->t_before, sizeof(double), 1, fp) != 1) {
        fclose(fp);
        fprintf(stderr, "Failed to read snapshot scalars: %s\n", path_snapshot);
        return -1;
    }

    state->n_adjustment = (int)n_adjustment;
    state->n_conduction = (int)n_conduction;
    state->n_save = (size_t)n_save;

    if (fread(state->m, sizeof(double), (size_t)shell_count, fp) != (size_t)shell_count ||
        fread(state->r, sizeof(double), (size_t)shell_count, fp) != (size_t)shell_count ||
        fread(state->rho, sizeof(double), (size_t)shell_count, fp) != (size_t)shell_count ||
        fread(state->p, sizeof(double), (size_t)shell_count, fp) != (size_t)shell_count) {
        fclose(fp);
        fprintf(stderr, "Failed to read snapshot arrays: %s\n", path_snapshot);
        return -1;
    }
    if (version >= GSIDM_VERSION &&
        fread(state->Mcy, sizeof(double), (size_t)shell_count, fp) != (size_t)shell_count) {
        fclose(fp);
        fprintf(stderr, "Failed to read snapshot Mcy array: %s\n", path_snapshot);
        return -1;
    }

    fclose(fp);

    /* recover derived parameters if halo save state is loaded */
    gsidm_update_derived_parameters(state);
    return 0;
}

static int gsidm_record_prepare(gsidm_record_t *record) {
    FILE *index;
    if (gsidm_ensure_dir_recursive(record->dir_data) != 0) {
        fprintf(stderr, "Failed to create output directory: %s\n", record->dir_data);
        return -1;
    }
    index = fopen(record->path_index, "rb");
    if (index != NULL) {
        fclose(index);
        return 0;
    }
    index = fopen(record->path_index, "wb");
    if (index == NULL) {
        fprintf(stderr, "Failed to create snapshot index: %s\n", record->path_index);
        return -1;
    }
    fprintf(index, "save_id,prefix,t,t_before,rho_center,n_adjustment,n_conduction,filename\n");
    fclose(index);
    return 0;
}

static int gsidm_write_initialization(
    const gsidm_state_t *state) {
    FILE *fp = fopen(state->record.path_init, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write initialization file: %s\n", state->record.path_init);
        return -1;
    }

    fprintf(fp, "# GravothermalSIDM simple C port\n");
    fprintf(fp, "# Default values are editable in gsidm_set_default_config() inside gravothermal_sidm.c\n");
    fprintf(fp, "run_name=%s\n", state->config->run_name);
    fprintf(fp, "output_dir=%s\n", state->config->output_dir);
    fprintf(fp, "profile=%s\n", gsidm_profile_name(state->config->profile));
    fprintf(fp, "r_s_kpc=%.17g\n", state->config->r_s);
    fprintf(fp, "rho_s_msun_per_pc3=%.17g\n", state->config->rho_s);
    fprintf(fp, "r_cut=%.17g\n", state->config->r_cut);
    fprintf(fp, "a=%.17g\n", state->config->a);
    fprintf(fp, "b=%.17g\n", state->config->b);
    fprintf(fp, "C=%.17g\n", state->config->C);
    fprintf(fp, "gamma=%.17g\n", state->config->gamma);
    fprintf(fp, "lmfp_model=%s\n", state->config->lmfp_model.raw_name);
    fprintf(fp, "smfp_model=%s\n", state->config->smfp_model.raw_name);
    fprintf(fp, "sigma_m_with_units_cm2_per_g=%.17g\n", state->config->sigma_m_with_units);
    fprintf(fp, "w_units_km_per_s=%.17g\n", state->config->w_units);
    fprintf(fp, "n_shells=%d\n", state->config->n_shells);
    fprintf(fp, "r_min=%.17g\n", state->config->r_min);
    fprintf(fp, "r_max=%.17g\n", state->config->r_max);
    fprintf(fp, "p_rmax_factor=%.17g\n", state->config->p_rmax_factor);
    fprintf(fp, "n_adjustment_fixed=%d\n", state->config->n_adjustment_fixed);
    fprintf(fp, "n_adjustment_max=%d\n", state->config->n_adjustment_max);
    fprintf(fp, "flag_timestep_use_relaxation=%d\n", state->config->flag_timestep_use_relaxation);
    fprintf(fp, "flag_timestep_use_energy=%d\n", state->config->flag_timestep_use_energy);
    fprintf(fp, "flag_hydrostatic_initial=%d\n", state->config->flag_hydrostatic_initial);

    /* save dimensionless quantities and unit scales */
    fprintf(fp, "sigma_m_dimensionless=%.17g\n", state->sigma_m);
    fprintf(fp, "w_dimensionless=%.17g\n", state->w);
    fprintf(fp, "t_relax_dimensionless=%.17g\n", state->t_relax);
    fprintf(fp, "scale_r_m=%.17g\n", state->scale_r);
    fprintf(fp, "scale_rho_kg_per_m3=%.17g\n", state->scale_rho);
    fprintf(fp, "scale_m_kg=%.17g\n", state->scale_m);
    fprintf(fp, "scale_u_m2_per_s2=%.17g\n", state->scale_u);
    fprintf(fp, "scale_p_pa=%.17g\n", state->scale_p);
    fprintf(fp, "scale_v_m_per_s=%.17g\n", state->scale_v);
    fprintf(fp, "scale_t_s=%.17g\n", state->scale_t);
    fprintf(fp, "scale_L_w=%.17g\n", state->scale_L);
    fprintf(fp, "scale_sigma_m_m2_per_kg=%.17g\n", state->scale_sigma_m);
    fprintf(fp, "use_openmp=%d\n", state->config->use_openmp);
    fprintf(fp, "omp_threads=%d\n", state->config->omp_threads);
    fprintf(fp, "openmp_min_shells=%d\n", GSIDM_OPENMP_MIN_SHELLS);
    fclose(fp);
    return 0;
}

static void gsidm_sanitize_time(double value, char *dst, size_t dst_size) {
    size_t i;
    snprintf(dst, dst_size, "%.9e", value);
    for (i = 0; dst[i] != '\0'; ++i) {
        if (dst[i] == '.') {
            dst[i] = 'd';
        } else if (dst[i] == '+') {
            dst[i] = 'p';
        } else if (dst[i] == '-') {
            dst[i] = 'm';
        }
    }
}

static int gsidm_write_binary_vector(FILE *fp, const double *values, size_t count) {
    if (count == 0) {
        return 0;
    }
    return fwrite(values, sizeof(double), count, fp) == count ? 0 : -1;
}

static double gsidm_projected_sphere_volume_factor(double radius, double projected_radius) {
    double radius2;                 /* spherical radius squared */
    double projected_radius2;       /* cylindrical radius squared */
    if (projected_radius >= radius) {
        return radius * radius * radius;
    }
    radius2 = radius * radius;
    projected_radius2 = projected_radius * projected_radius;
    return radius * radius * radius - pow(radius2 - projected_radius2, 1.5);
}

/* obtain cylindrical mass profile */
static void gsidm_update_cylindrical_mass(gsidm_state_t *state) {
    int i;
    int j;
    int n = state->config->n_shells;
    int enable_omp = gsidm_parallel_enabled(state->config);

#ifdef _OPENMP
#pragma omp parallel for private(j) if(enable_omp) schedule(static)
#endif
    for (i = 0; i < n; ++i) {
        double projected_radius = state->r[i]; /* cylindrical aperture radius R */
        double mass_cyl = 0.0;                 /* projected mass inside R */

        for (j = 0; j < n; ++j) {
            double r_inner = j == 0 ? 0.0 : state->r[j - 1]; /* shell inner radius */
            double r_outer = state->r[j];                    /* shell outer radius */
            double outer_factor;                             /* projected outer volume factor */
            double inner_factor;                             /* projected inner volume factor */

            if (projected_radius <= 0.0 || r_outer <= 0.0 || state->rho[j] <= 0.0) {
                continue;
            }

            outer_factor = gsidm_projected_sphere_volume_factor(r_outer, projected_radius);
            inner_factor = r_inner > 0.0 ? gsidm_projected_sphere_volume_factor(r_inner, projected_radius) : 0.0;
            mass_cyl += state->rho[j] * (outer_factor - inner_factor) / 3.0;
        }

        state->Mcy[i] = mass_cyl;
    }
}

static int gsidm_save_snapshot(
    const gsidm_record_t *record,
    const char *prefix,
    const gsidm_snapshot_t *snapshot) {
    char time_text[128];
    char filename[GSIDM_PATH_MAX];
    char path_snapshot[GSIDM_PATH_MAX];
    uint64_t magic = GSIDM_MAGIC;
    uint64_t version = GSIDM_VERSION;
    uint64_t shell_count = (uint64_t)snapshot->n_shells;
    int64_t n_adjustment = (int64_t)snapshot->n_adjustment;
    int64_t n_conduction = (int64_t)snapshot->n_conduction;
    uint64_t n_save = (uint64_t)snapshot->n_save;
    FILE *fp;
    FILE *index;

    gsidm_sanitize_time(snapshot->t, time_text, sizeof(time_text));
    snprintf(filename, sizeof(filename), "%s_%zu_t_%s.bin", prefix, snapshot->n_save, time_text);
    gsidm_join_path(path_snapshot, sizeof(path_snapshot), record->dir_data, filename);

    fp = fopen(path_snapshot, "wb");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write snapshot: %s\n", path_snapshot);
        return -1;
    }

    /* save data */
    fwrite(&magic, sizeof(uint64_t), 1, fp);
    fwrite(&version, sizeof(uint64_t), 1, fp);
    fwrite(&shell_count, sizeof(uint64_t), 1, fp);
    fwrite(&snapshot->rho_center, sizeof(double), 1, fp);
    fwrite(&n_adjustment, sizeof(int64_t), 1, fp);
    fwrite(&n_conduction, sizeof(int64_t), 1, fp);
    fwrite(&n_save, sizeof(uint64_t), 1, fp);
    fwrite(&snapshot->t_epsilon, sizeof(double), 1, fp);
    fwrite(&snapshot->r_epsilon, sizeof(double), 1, fp);
    fwrite(&snapshot->t, sizeof(double), 1, fp);
    fwrite(&snapshot->t_before, sizeof(double), 1, fp);
    if (gsidm_write_binary_vector(fp, snapshot->m, snapshot->n_shells) != 0 ||
        gsidm_write_binary_vector(fp, snapshot->r, snapshot->n_shells) != 0 ||
        gsidm_write_binary_vector(fp, snapshot->rho, snapshot->n_shells) != 0 ||
        gsidm_write_binary_vector(fp, snapshot->p, snapshot->n_shells) != 0 ||
        gsidm_write_binary_vector(fp, snapshot->Mcy, snapshot->n_shells) != 0) {
        fclose(fp);
        fprintf(stderr, "Failed while writing snapshot arrays: %s\n", path_snapshot);
        return -1;
    }
    fclose(fp);

    index = fopen(record->path_index, "a");
    if (index == NULL) {
        fprintf(stderr, "Failed to append snapshot index: %s\n", record->path_index);
        return -1;
    }
    fprintf(
        index,
        "%zu,%s,%.17g,%.17g,%.17g,%d,%d,%s\n",
        snapshot->n_save,
        prefix,
        snapshot->t,
        snapshot->t_before,
        snapshot->rho_center,
        snapshot->n_adjustment,
        snapshot->n_conduction,
        filename);
    fclose(index);
    return 0;
}

static int gsidm_state_alloc_arrays(gsidm_state_t *state, int n_shells) {
    size_t n = (size_t)n_shells;
    state->r = (double *)calloc(n, sizeof(double));
    state->m = (double *)calloc(n, sizeof(double));
    state->Mcy = (double *)calloc(n, sizeof(double));
    state->rho = (double *)calloc(n, sizeof(double));
    state->p = (double *)calloc(n, sizeof(double));
    state->u = (double *)calloc(n, sizeof(double));
    state->v = (double *)calloc(n, sizeof(double));
    state->L = (double *)calloc(n, sizeof(double));
    state->Kn = (double *)calloc(n, sizeof(double));
    state->Kinv_lmfp = (double *)calloc(n, sizeof(double));
    state->Kinv_smfp = (double *)calloc(n, sizeof(double));
    state->r_ext = (double *)calloc(n + 1, sizeof(double));
    state->m_ext = (double *)calloc(n + 1, sizeof(double));
    state->a_arr = (double *)calloc(n - 2, sizeof(double));
    state->b_arr = (double *)calloc(n - 1, sizeof(double));
    state->c_arr = (double *)calloc(n - 2, sizeof(double));
    state->d_arr = (double *)calloc(n - 1, sizeof(double));
    state->delta_r = (double *)calloc(n, sizeof(double));
    state->delta_uc = (double *)calloc(n, sizeof(double));
    state->delta_rho = (double *)calloc(n, sizeof(double));
    state->delta_ph = (double *)calloc(n, sizeof(double));
    state->keff = (double *)calloc(n, sizeof(double));
    state->tdma_ba = (double *)calloc(n > 0 ? n - 1 : 0, sizeof(double));
    state->tdma_dc = (double *)calloc(n > 0 ? n - 1 : 0, sizeof(double));
    state->r2_ext = (double *)calloc(n + 1, sizeof(double));
    state->r3_ext = (double *)calloc(n + 1, sizeof(double));
    state->r2_times_dr = (double *)calloc(n, sizeof(double));
    state->r3_shell = (double *)calloc(n, sizeof(double));

    if (state->r == NULL || state->m == NULL || state->Mcy == NULL || state->rho == NULL || state->p == NULL ||
        state->u == NULL || state->v == NULL || state->L == NULL || state->Kn == NULL ||
        state->Kinv_lmfp == NULL || state->Kinv_smfp == NULL || state->r_ext == NULL ||
        state->m_ext == NULL || state->a_arr == NULL || state->b_arr == NULL ||
        state->c_arr == NULL || state->d_arr == NULL || state->delta_r == NULL ||
        state->delta_uc == NULL || state->delta_rho == NULL || state->delta_ph == NULL ||
        state->keff == NULL || state->tdma_ba == NULL || state->tdma_dc == NULL ||
        state->r2_ext == NULL || state->r3_ext == NULL || state->r2_times_dr == NULL ||
        state->r3_shell == NULL) {
        return -1;
    }
    return 0;
}

static void gsidm_state_free_arrays(gsidm_state_t *state) {
    free(state->r);
    free(state->m);
    free(state->Mcy);
    free(state->rho);
    free(state->p);
    free(state->u);
    free(state->v);
    free(state->L);
    free(state->Kn);
    free(state->Kinv_lmfp);
    free(state->Kinv_smfp);
    free(state->r_ext);
    free(state->m_ext);
    free(state->a_arr);
    free(state->b_arr);
    free(state->c_arr);
    free(state->d_arr);
    free(state->delta_r);
    free(state->delta_uc);
    free(state->delta_rho);
    free(state->delta_ph);
    free(state->keff);
    free(state->tdma_ba);
    free(state->tdma_dc);
    free(state->r2_ext);
    free(state->r3_ext);
    free(state->r2_times_dr);
    free(state->r3_shell);
}

static double gsidm_truncated_mass_integrand(double radius, void *ctx_ptr) {
    const gsidm_config_t *config = (const gsidm_config_t *)ctx_ptr;
    if (radius <= 0.0) {
        return 0.0;
    }
    if (config->profile == GSIDM_PROFILE_TNFW) {
        return radius * exp(-radius / config->r_cut) / pow(1.0 + radius, 2.0);
    }
    if (config->profile == GSIDM_PROFILE_THERN) {
        return radius * exp(-radius / config->r_cut) / pow(1.0 + radius, 3.0);
    }
    return 0.0;
}

/* obtain halo mass profile */
static double gsidm_get_initial_mass_scalar(const gsidm_config_t *config, double x) {
    if (config->profile == GSIDM_PROFILE_NFW) {
        return -x / (1.0 + x) + log(1.0 + x);
    }
    if (config->profile == GSIDM_PROFILE_HERNQUIST) {
        return x * x / (2.0 * pow(1.0 + x, 2.0));
    }
    if (config->profile == GSIDM_PROFILE_ISOTHERMAL) {
        return x - atan(x);
    }
    if (config->profile == GSIDM_PROFILE_TNFW) {
        if (x <= 0.0) {
            return 0.0;
        }
        return gsidm_adaptive_simpson(gsidm_truncated_mass_integrand, (void *)config, 0.0, x, 1.0e-10);
    }
    if (config->profile == GSIDM_PROFILE_THERN) {
        if (x <= 0.0) {
            return 0.0;
        }
        return gsidm_adaptive_simpson(gsidm_truncated_mass_integrand, (void *)config, 0.0, x, 1.0e-10);
    }
    return 0.0;
}

/* obtain halo density profile */
static double gsidm_get_initial_rho_scalar(const gsidm_config_t *config, double x) {
    if (config->profile == GSIDM_PROFILE_NFW) {
        return 1.0 / (x * pow(1.0 + x, 2.0));
    }
    if (config->profile == GSIDM_PROFILE_HERNQUIST) {
        return 1.0 / (x * pow(1.0 + x, 3.0));
    }
    if (config->profile == GSIDM_PROFILE_ISOTHERMAL) {
        return 1.0 / (1.0 + x * x);
    }
    if (config->profile == GSIDM_PROFILE_TNFW) {
        return exp(-x / config->r_cut) / (x * pow(1.0 + x, 2.0));
    }
    if (config->profile == GSIDM_PROFILE_THERN) {
        return exp(-x / config->r_cut) / (x * pow(1.0 + x, 3.0));
    }
    return 0.0;
}

static double gsidm_pressure_integrand(double radius, void *ctx_ptr) {
    gsidm_pressure_context_t *ctx = (gsidm_pressure_context_t *)ctx_ptr;
    return gsidm_get_initial_mass_scalar(ctx->config, radius) *
           gsidm_get_initial_rho_scalar(ctx->config, radius) /
           (radius * radius);
}

static double gsidm_simpson(
    double (*fn)(double, void *),
    void *ctx,
    double a,
    double b) {
    double midpoint = 0.5 * (a + b);
    return (b - a) * (fn(a, ctx) + 4.0 * fn(midpoint, ctx) + fn(b, ctx)) / 6.0;
}

static double gsidm_adaptive_simpson_recursive(
    double (*fn)(double, void *),
    void *ctx,
    double a,
    double b,
    double tol,
    double whole,
    int depth) {
    double midpoint = 0.5 * (a + b);
    double left = gsidm_simpson(fn, ctx, a, midpoint);
    double right = gsidm_simpson(fn, ctx, midpoint, b);
    double delta = left + right - whole;
    if (depth <= 0 || fabs(delta) <= 15.0 * tol) {
        return left + right + delta / 15.0;
    }
    return gsidm_adaptive_simpson_recursive(fn, ctx, a, midpoint, tol * 0.5, left, depth - 1) +
           gsidm_adaptive_simpson_recursive(fn, ctx, midpoint, b, tol * 0.5, right, depth - 1);
}

static double gsidm_adaptive_simpson(
    double (*fn)(double, void *),
    void *ctx,
    double a,
    double b,
    double tol) {
    double whole = gsidm_simpson(fn, ctx, a, b);
    return gsidm_adaptive_simpson_recursive(fn, ctx, a, b, tol, whole, 20);
}

/* obtain halo pressure profile */
static double gsidm_get_initial_pressure_scalar(const gsidm_config_t *config, double x) {
    if (config->profile == GSIDM_PROFILE_HERNQUIST) {
        return (log(1.0 + 1.0 / x) -
                (25.0 + 2.0 * x * (26.0 + 3.0 * x * (7.0 + 2.0 * x))) /
                    (12.0 * pow(1.0 + x, 4.0))) /
               2.0;
    }
    if (config->profile == GSIDM_PROFILE_ISOTHERMAL) {
        return (M_PI * M_PI) / 8.0 - atan(x) * (2.0 + x * atan(x)) / (2.0 * x);
    }
    if (config->profile == GSIDM_PROFILE_NFW ||
        config->profile == GSIDM_PROFILE_TNFW ||
        config->profile == GSIDM_PROFILE_THERN) {
        gsidm_pressure_context_t ctx;
        double xmax = config->p_rmax_factor * config->r_max; /* pressure integral upper radius */
        ctx.config = config;
        return gsidm_adaptive_simpson(gsidm_pressure_integrand, &ctx, x, xmax, 1.0e-10);
    }
    return 0.0;
}

static int gsidm_parallel_enabled(const gsidm_config_t *config) {
    return config->use_openmp && config->n_shells >= GSIDM_OPENMP_MIN_SHELLS;
}

/* given r, mass, rho, and p, set derived halo quantities */
static void gsidm_update_derived_parameters(gsidm_state_t *state) {
    int i;
    int n = state->config->n_shells;
    int enable_omp = gsidm_parallel_enabled(state->config);

#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 0; i < n; ++i) {
        double scaled_velocity;     /* v/w for velocity-dependent scattering */

        /* specific energy */
        state->u[i] = 1.5 * state->p[i] / state->rho[i];

        /* 1d velocity dispersion */
        state->v[i] = sqrt(state->p[i] / state->rho[i]);
        scaled_velocity = state->v[i] / state->w;

        /* effective thermal conductivity */
        state->Kinv_smfp[i] =
            state->sigma_m * gsidm_eval_scattering_model(&state->config->smfp_model, scaled_velocity) /
            (state->config->b * state->v[i]);
        state->Kinv_lmfp[i] =
            1.0 /
            (state->config->a * state->config->C * state->v[i] * state->p[i] * state->sigma_m *
             gsidm_eval_scattering_model(&state->config->lmfp_model, scaled_velocity));
        state->keff[i] = 1.0 / (state->Kinv_smfp[i] + state->Kinv_lmfp[i]);

        /* Knudsen number */
        state->Kn[i] = 1.0 / (sqrt(state->p[i]) * state->sigma_m);
    }

    if (n > 1) {
        /* luminosity */
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
        for (i = 1; i < n - 1; ++i) {
            state->L[i] = -state->r[i] * state->r[i] * (state->keff[i] + state->keff[i + 1]) / 2.0 *
                          (state->u[i + 1] - state->u[i]) /
                          ((state->r[i + 1] - state->r[i - 1]) / 2.0);
        }
        state->L[0] = -state->r[0] * state->r[0] * (state->keff[0] + state->keff[1]) / 2.0 *
                      (state->u[1] - state->u[0]) / (state->r[1] / 2.0);
    } else {
        state->L[0] = 0.0;
    }
    state->L[n - 1] = 0.0;
}

static double gsidm_get_central_quantity(const double *values) {
    return values[GSIDM_CENTRAL_INDEX];
}

static void gsidm_update_timestep_min(double candidate, double *current_min) {
    if (isfinite(candidate) && candidate > 0.0 && candidate < *current_min) {
        *current_min = candidate;
    }
}

/* determine heat conduction time step */
static double gsidm_get_timestep(const gsidm_state_t *state) {
    int i;
    double delta_t1 = DBL_MAX;      /* energy-change timestep */
    double delta_t2 = DBL_MAX;      /* relaxation timestep */
    double delta_t = 0.0;           /* chosen heat conduction timestep */

    if (state->m[0] > 0.0) {
        double rate = state->L[0] / state->m[0]; /* central du/dt */
        if (rate != 0.0) {
            gsidm_update_timestep_min(fabs(state->u[0] / rate), &delta_t1);
        }
    }

    for (i = 1; i < state->config->n_shells; ++i) {
        double dm = state->m[i] - state->m[i - 1]; /* shell mass */
        if (dm > 0.0) {
            double rate = (state->L[i] - state->L[i - 1]) / dm; /* shell du/dt */
            if (rate != 0.0) {
                gsidm_update_timestep_min(fabs(state->u[i] / rate), &delta_t1);
            }
        }
    }

    for (i = 0; i < state->config->n_shells; ++i) {
        if (state->rho[i] > 0.0 && state->v[i] > 0.0) {
            gsidm_update_timestep_min(1.0 / (state->rho[i] * state->v[i]), &delta_t2);
        }
    }

    if (state->config->flag_timestep_use_relaxation && state->config->flag_timestep_use_energy) {
        delta_t = delta_t1 < delta_t2 ? delta_t1 : delta_t2;
    } else if (state->config->flag_timestep_use_relaxation) {
        delta_t = delta_t2;
    } else if (state->config->flag_timestep_use_energy) {
        delta_t = delta_t1;
    }

    if (!isfinite(delta_t) || delta_t <= 0.0 || delta_t == DBL_MAX) {
        fprintf(
            stderr,
            "Failed to compute a finite timestep: energy=%g, relaxation=%g. "
            "Check the initial pressure/density profile and timestep flags.\n",
            delta_t1,
            delta_t2);
        exit(1);
    }

    return state->t_epsilon * delta_t;
}

/* perform heat conduction step */
static void gsidm_conduct_heat(gsidm_state_t *state, double max_delta_t) {
    int i;
    int enable_omp = gsidm_parallel_enabled(state->config);
    double delta_t;                 /* heat conduction timestep */

    /* update conduction counter */
    state->n_conduction += 1;

    /* determine minimum time step */
    delta_t = gsidm_get_timestep(state);
    if (max_delta_t > 0.0 && isfinite(max_delta_t) && delta_t > max_delta_t) {
        delta_t = max_delta_t;
    }
    if (!isfinite(delta_t) || delta_t <= 0.0) {
        fprintf(stderr, "Invalid heat-conduction timestep: %g\n", delta_t);
        exit(1);
    }

    /* calculate delta_uc */
    state->delta_uc[0] = -delta_t * (state->L[0] / state->m[0]);
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 1; i < state->config->n_shells; ++i) {
        state->delta_uc[i] =
            -delta_t * ((state->L[i] - state->L[i - 1]) / (state->m[i] - state->m[i - 1]));
    }

#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 0; i < state->config->n_shells; ++i) {
        double delta_pc = state->p[i] * state->delta_uc[i] / state->u[i]; /* pressure change by conduction */

        /* update variables */
        state->p[i] += delta_pc;
        state->u[i] += state->delta_uc[i];
    }

    /* update the current dimensionless time */
    state->t_before = state->t;
    state->t += delta_t;
}

static void gsidm_tdma_solver(
    const double *aa,
    const double *ba_in,
    const double *ca,
    const double *da_in,
    int nf,
    double *ba,
    double *dc,
    double *x_out) {
    int i;
    for (i = 0; i < nf; ++i) {
        ba[i] = ba_in[i];
        dc[i] = -da_in[i];
    }
    for (i = 1; i < nf; ++i) {
        double wa = aa[i - 1] / ba[i - 1];
        ba[i] = ba[i] - wa * ca[i - 1];
        dc[i] = dc[i] - wa * dc[i - 1];
    }
    x_out[nf - 1] = dc[nf - 1] / ba[nf - 1];
    for (i = nf - 2; i >= 0; --i) {
        x_out[i] = (dc[i] - ca[i] * x_out[i + 1]) / ba[i];
    }
}

static void gsidm_set_hydrostatic_coefficients(gsidm_state_t *state) {
    int i;
    int j;
    int n = state->config->n_shells;
    int enable_omp = gsidm_parallel_enabled(state->config);

    /* fill extended r and m arrays */
    state->r_ext[0] = 0.0;
    state->m_ext[0] = 0.0;
    for (i = 0; i < n; ++i) {
        state->r_ext[i + 1] = state->r[i];
        state->m_ext[i + 1] = state->m[i];
    }

    for (i = 0; i < n + 1; ++i) {
        state->r2_ext[i] = state->r_ext[i] * state->r_ext[i];
        state->r3_ext[i] = state->r2_ext[i] * state->r_ext[i];
    }

    /* set tridiagonal matrix coefficients */
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (j = 0; j < n - 2; ++j) {
        state->a_arr[j] =
            (12.0 * state->config->gamma * state->p[j + 1] * state->r2_ext[j + 1] * state->r2_ext[j + 2] +
             state->m_ext[j + 2] *
                 (-3.0 * state->r2_ext[j + 1] * state->r_ext[j + 3] * state->rho[j + 1] +
                  state->r3_ext[j + 1] * (2.0 * state->rho[j + 1] - state->rho[j + 2]) +
                  state->r3_ext[j + 2] * (state->rho[j + 1] + state->rho[j + 2]))) /
            (4.0 * (state->r3_ext[j + 1] - state->r3_ext[j + 2]));
    }

#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (j = 0; j < n - 2; ++j) {
        state->c_arr[j] =
            (12.0 * state->config->gamma * state->p[j + 1] * state->r2_ext[j + 1] * state->r2_ext[j + 2] +
             state->m_ext[j + 1] *
                 (state->r3_ext[j + 1] * (state->rho[j] + state->rho[j + 1]) -
                  state->r2_ext[j + 2] *
                      (state->r_ext[j + 2] * state->rho[j] +
                       3.0 * state->r_ext[j] * state->rho[j + 1] -
                       2.0 * state->r_ext[j + 2] * state->rho[j + 1]))) /
            (4.0 * (state->r3_ext[j + 1] - state->r3_ext[j + 2]));
    }

#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (j = 0; j < n - 1; ++j) {
        state->b_arr[j] =
            state->r_ext[j + 1] *
            (-4.0 * state->p[j] *
                 (2.0 * state->r3_ext[j] + (-2.0 + 3.0 * state->config->gamma) * state->r3_ext[j + 1]) *
                 (state->r3_ext[j + 1] - state->r3_ext[j + 2]) -
             4.0 * state->p[j + 1] * (state->r3_ext[j] - state->r3_ext[j + 1]) *
                 ((-2.0 + 3.0 * state->config->gamma) * state->r3_ext[j + 1] + 2.0 * state->r3_ext[j + 2]) +
             3.0 * state->m_ext[j + 1] * state->r_ext[j + 1] *
                 (state->r_ext[j] - state->r_ext[j + 2]) *
                 (state->r3_ext[j + 2] * state->rho[j] + state->r3_ext[j] * state->rho[j + 1] -
                  state->r3_ext[j + 1] * (state->rho[j] + state->rho[j + 1]))) /
            (4.0 * (state->r3_ext[j] - state->r3_ext[j + 1]) *
             (state->r3_ext[j + 1] - state->r3_ext[j + 2]));

        state->d_arr[j] = -state->p[j] * state->r2_ext[j + 1] + state->p[j + 1] * state->r2_ext[j + 1] +
                          0.25 * state->m_ext[j + 1] *
                              (-state->r_ext[j] + state->r_ext[j + 2]) *
                              (state->rho[j] + state->rho[j + 1]);
    }
}

static void gsidm_hydrostatic_adjustment_step(gsidm_state_t *state) {
    int i;
    int n = state->config->n_shells;
    int enable_omp = gsidm_parallel_enabled(state->config);

    /* update the counter */
    state->n_adjustment += 1;

    gsidm_set_hydrostatic_coefficients(state);
    gsidm_tdma_solver(
        state->a_arr,
        state->b_arr,
        state->c_arr,
        state->d_arr,
        n - 1,
        state->tdma_ba,
        state->tdma_dc,
        state->delta_r);

    state->delta_r[n - 1] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 0; i < n; ++i) {
        state->r2_times_dr[i] = state->r[i] * state->r[i] * state->delta_r[i];
        state->r3_shell[i] = state->r[i] * state->r[i] * state->r[i];
    }

    /* update delta_rho */
    state->delta_rho[0] = -state->rho[0] * state->r2_times_dr[0] / (state->r3_shell[0] / 3.0);

    /* update delta_ph */
    state->delta_ph[0] =
        -state->config->gamma * state->p[0] * state->r2_times_dr[0] / (state->r3_shell[0] / 3.0);

#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 1; i < n; ++i) {
        double shell_volume_factor = (state->r3_shell[i] - state->r3_shell[i - 1]) / 3.0; /* shell volume / 4pi */
        state->delta_rho[i] =
            -state->rho[i] * (state->r2_times_dr[i] - state->r2_times_dr[i - 1]) / shell_volume_factor;
        state->delta_ph[i] =
            -state->config->gamma * state->p[i] *
            (state->r2_times_dr[i] - state->r2_times_dr[i - 1]) / shell_volume_factor;
    }

    /* update r, rho, and p */
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 0; i < n; ++i) {
        state->r[i] += state->delta_r[i];
        state->rho[i] += state->delta_rho[i];
        state->p[i] += state->delta_ph[i];
    }
}

static void gsidm_hydrostatic_adjustment(gsidm_state_t *state) {
    int i;

    /* initialize number of adjustment steps taken */
    state->n_adjustment = 0;

    if (state->config->n_adjustment_fixed == 0 || state->config->n_adjustment_max == 0) {
        /* no-op */
    } else if (state->config->n_adjustment_fixed > 0) {
        /* perform fixed number of adjustments */
        for (i = 0; i < state->config->n_adjustment_fixed; ++i) {
            gsidm_hydrostatic_adjustment_step(state);
        }
    } else {
        gsidm_hydrostatic_adjustment_step(state);
        while (1) {
            double max_relative_change = 0.0;

            /* adjust until convergence is met */
            for (i = 0; i < state->config->n_shells; ++i) {
                double value = fabs(state->delta_r[i] / state->r[i]);
                if (value > max_relative_change) {
                    max_relative_change = value;
                }
            }
            if (max_relative_change <= state->r_epsilon) {
                break;
            }
            if (state->config->n_adjustment_max >= 0 &&
                state->n_adjustment == state->config->n_adjustment_max) {
                /* break if max adjustments is reached */
                break;
            }
            gsidm_hydrostatic_adjustment_step(state);
        }
    }

    gsidm_update_derived_parameters(state);
}

static int gsidm_save_halo(gsidm_state_t *state, const char *prefix) {
    gsidm_snapshot_t snapshot;

    /* update the counter if prefix is default */
    if (gsidm_string_equals(prefix, GSIDM_PREFIX_DEFAULT)) {
        state->n_save += 1;
    }

    /* calculate derived output quantities */
    gsidm_update_cylindrical_mass(state);

    /* save data */
    snapshot.rho_center = state->rho_center;
    snapshot.n_adjustment = state->n_adjustment;
    snapshot.n_conduction = state->n_conduction;
    snapshot.n_save = state->n_save;
    snapshot.t_epsilon = state->t_epsilon;
    snapshot.r_epsilon = state->r_epsilon;
    snapshot.t = state->t;
    snapshot.t_before = state->t_before;
    snapshot.m = state->m;
    snapshot.r = state->r;
    snapshot.rho = state->rho;
    snapshot.p = state->p;
    snapshot.Mcy = state->Mcy;
    snapshot.n_shells = (size_t)state->config->n_shells;
    return gsidm_save_snapshot(&state->record, prefix, &snapshot);
}

static int gsidm_initialize_state(gsidm_state_t *state, const gsidm_config_t *config) {
    int i;
    int n = config->n_shells;
    int enable_omp = gsidm_parallel_enabled(config);
    double log_r_min;               /* log10 of innermost shell radius */
    double log_r_max;               /* log10 of outer halo radius */
    double *r_mid;                  /* shell midpoint radius */

    memset(state, 0, sizeof(*state));
    state->config = config;

    gsidm_record_init(&state->record, config);
    if (gsidm_record_prepare(&state->record) != 0) {
        return -1;
    }
    if (gsidm_state_alloc_arrays(state, config->n_shells) != 0) {
        fprintf(stderr, "Failed to allocate simulation arrays\n");
        return -1;
    }

    state->scale_r = config->r_s * GSIDM_KPC_METERS; /* radius unit, m */
    state->scale_rho = config->rho_s * GSIDM_SOLAR_MASS_KG / pow(GSIDM_PARSEC_METERS, 3.0); /* density unit */
    state->scale_m = 4.0 * M_PI * state->scale_rho * pow(state->scale_r, 3.0); /* mass unit */
    state->scale_u = GSIDM_G * state->scale_m / state->scale_r; /* specific energy unit */
    state->scale_p = state->scale_u * state->scale_rho; /* pressure unit */
    state->scale_v = sqrt(state->scale_u); /* velocity unit */
    state->scale_t = 1.0 / sqrt(4.0 * M_PI * state->scale_rho * GSIDM_G); /* time unit */
    state->scale_L = GSIDM_G * state->scale_m * state->scale_m / (state->scale_r * state->scale_t); /* luminosity unit */
    state->scale_sigma_m = 1.0 / (state->scale_r * state->scale_rho); /* sigma/m unit */

    /* convert dimensionful quantities to be dimensionless */
    state->sigma_m =
        config->sigma_m_with_units * GSIDM_CM2_PER_G_TO_M2_PER_KG / state->scale_sigma_m;
    state->w = config->w_units * GSIDM_KM_PER_S_TO_M_PER_S / state->scale_v; /* velocity scale */
    state->t_relax =
        (2.0 /
         (3.0 * config->a * config->C *
          config->sigma_m_with_units * GSIDM_CM2_PER_G_TO_M2_PER_KG *
          gsidm_eval_scattering_model(&config->lmfp_model, 1.0 / state->w) *
          state->scale_v * state->scale_rho)) /
        state->scale_t;

    /* set initial halo radius, mass, density, and pressure profiles */

    /* initial radius, location of outermost edge of each shell */
    log_r_min = log10(config->r_min);
    log_r_max = log10(config->r_max);
    for (i = 0; i < n; ++i) {
        double fraction = n == 1 ? 0.0 : (double)i / (double)(n - 1);
        state->r[i] = pow(10.0, log_r_min + fraction * (log_r_max - log_r_min));
    }

    /* location of midpoints of shells */
    r_mid = (double *)calloc((size_t)n, sizeof(double));
    if (r_mid == NULL) {
        fprintf(stderr, "Out of memory creating shell midpoints\n");
        return -1;
    }
    r_mid[0] = state->r[0] / 2.0;
    for (i = 1; i < n; ++i) {
        r_mid[i] = 0.5 * (state->r[i - 1] + state->r[i]);
    }

    /* set mass */
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 0; i < n; ++i) {
        state->m[i] = gsidm_get_initial_mass_scalar(config, state->r[i]);
    }

    /* set rho */
    state->rho[0] = 3.0 * state->m[0] / pow(state->r[0], 3.0);
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 1; i < n; ++i) {
        state->rho[i] = gsidm_get_initial_rho_scalar(config, r_mid[i]);
    }

    /* set pressure */
#ifdef _OPENMP
#pragma omp parallel for if(enable_omp) schedule(static)
#endif
    for (i = 0; i < n; ++i) {
        state->p[i] = gsidm_get_initial_pressure_scalar(config, r_mid[i]);
    }
    state->p[0] = gsidm_get_initial_pressure_scalar(config, state->r[0]);

    /* store initial density of innermost shell */
    state->rho_center = gsidm_get_central_quantity(state->rho);
    state->t = 0.0;
    state->t_before = 0.0;
    state->n_save = 0;
    state->n_adjustment = 0;
    state->n_conduction = 0;
    state->t_epsilon = 0.0;
    state->r_epsilon = 0.0;

    /* with necessary quantities initialized, set derived quantities */
    gsidm_update_derived_parameters(state);
    free(r_mid);
    return 0;
}

static void gsidm_destroy_state(gsidm_state_t *state) {
    gsidm_state_free_arrays(state);
}

int gsidm_run(const gsidm_config_t *config) {
    gsidm_config_t runtime_config;
    gsidm_state_t state;
    double last_save_time;          /* time of previous snapshot */
    double last_save_rho;           /* central density at previous snapshot */
    int rc = 0;
    int resumed = 0;
    char resume_snapshot_path[GSIDM_PATH_MAX];

    runtime_config = *config;

    /* load saved initialization if resuming */
    if (runtime_config.resume_enabled) {
        gsidm_record_t record;
        gsidm_record_init(&record, &runtime_config);
        if (gsidm_apply_saved_initialization(&runtime_config, &record) != 0) {
            return 1;
        }
    }

    if (runtime_config.n_shells <= GSIDM_CENTRAL_INDEX) {
        fprintf(stderr, "n_shells must be greater than %d\n", GSIDM_CENTRAL_INDEX);
        return 1;
    }
    if ((runtime_config.profile == GSIDM_PROFILE_TNFW ||
         runtime_config.profile == GSIDM_PROFILE_THERN) &&
        runtime_config.r_cut <= 0.0) {
        fprintf(stderr, "r_cut must be positive for truncated profiles\n");
        return 1;
    }
    if (!runtime_config.flag_timestep_use_relaxation && !runtime_config.flag_timestep_use_energy) {
        fprintf(stderr, "At least one time-step criterion must be enabled\n");
        return 1;
    }
    if (runtime_config.n_adjustment_fixed > -1 && runtime_config.n_adjustment_max > -1) {
        fprintf(stderr, "Cannot set both n_adjustment_fixed and n_adjustment_max\n");
        return 1;
    }

#ifdef _OPENMP
    if (runtime_config.omp_threads > 0) {
        omp_set_num_threads(runtime_config.omp_threads);
    }
#else
    if (runtime_config.omp_threads > 0 || runtime_config.use_openmp) {
        fprintf(stderr, "Note: compiled without OpenMP support; running single-threaded.\n");
    }
#endif

    if (gsidm_initialize_state(&state, &runtime_config) != 0) {
        gsidm_destroy_state(&state);
        return 1;
    }

    if (runtime_config.resume_enabled) {
        if (gsidm_resolve_resume_snapshot_path(
                &state.record,
                &runtime_config,
                resume_snapshot_path,
                sizeof(resume_snapshot_path)) != 0) {
            gsidm_destroy_state(&state);
            return 1;
        }
        if (gsidm_load_snapshot_into_state(&state, resume_snapshot_path) != 0) {
            gsidm_destroy_state(&state);
            return 1;
        }
        resumed = 1;
        printf("Resumed from snapshot: %s\n", resume_snapshot_path);
    }

    state.t_epsilon = runtime_config.t_epsilon;
    state.r_epsilon = runtime_config.r_epsilon;

    if (!runtime_config.has_save_frequency_rate &&
        !runtime_config.has_save_frequency_timing &&
        !runtime_config.has_save_frequency_density) {
        printf("WARNING: no intermediate halo files will be saved\n");
    }

    if (!resumed) {
        /* initialization for new evolution calculation */

        /* save halo initialization file */
        if (gsidm_write_initialization(&state) != 0) {
            gsidm_destroy_state(&state);
            return 1;
        }

        if (runtime_config.flag_hydrostatic_initial) {
            /* enforce hydrostatic condition initially */
            gsidm_hydrostatic_adjustment(&state);
        }
        /* save initial state */
        if (gsidm_save_halo(&state, GSIDM_PREFIX_DEFAULT) != 0) {
            gsidm_destroy_state(&state);
            return 1;
        }
    }

    last_save_time = state.t;
    last_save_rho = gsidm_get_central_quantity(state.rho);

    /* evolve the halo */
    while (1) {
        double current_rho;         /* current central density */
        double smallest_kn;         /* minimum Knudsen number in halo */
        int i;
        int save_snapshot_now = 0;

        if (state.t >= runtime_config.t_end) {
            if (gsidm_save_halo(&state, GSIDM_PREFIX_DEFAULT) != 0) {
                rc = 1;
            }
            printf("Success: evolution reached dimensionless time %.17g\n", state.t);
            break;
        }

        /* conduct heat */
        gsidm_conduct_heat(&state, runtime_config.t_end - state.t);

        /* make hydrostatic adjustment */
        gsidm_hydrostatic_adjustment(&state);

        /* determine if halo has met conditions for terminating evolution */
        current_rho = gsidm_get_central_quantity(state.rho);
        if (state.t >= runtime_config.t_end) {
            gsidm_save_halo(&state, GSIDM_PREFIX_DEFAULT);
            printf("Success: evolution reached dimensionless time %.17g\n", state.t);
            break;
        }
        if (current_rho > runtime_config.rho_factor_end * state.rho_center) {
            gsidm_save_halo(&state, GSIDM_PREFIX_DEFAULT);
            printf("Success: central density reached %.17g times its initial value\n", runtime_config.rho_factor_end);
            break;
        }
        if (runtime_config.Kn_end > 0.0) {
            smallest_kn = state.Kn[0];
            for (i = 1; i < runtime_config.n_shells; ++i) {
                if (state.Kn[i] < smallest_kn) {
                    smallest_kn = state.Kn[i];
                }
            }
            if (smallest_kn < runtime_config.Kn_end) {
                gsidm_save_halo(&state, GSIDM_PREFIX_DEFAULT);
                printf("Success: smallest Knudsen number in halo is %.17g\n", smallest_kn);
                break;
            }
        }

        /* save halo state at intermediate stages of evolution */
        if (runtime_config.has_save_frequency_rate) {
            /* save based on relaxation time scale */
            save_snapshot_now =
                (state.t / state.t_relax) >= ((double)state.n_save / runtime_config.save_frequency_rate);
        } else if (runtime_config.has_save_frequency_timing) {
            /* save based on relative time increase */
            double delta_t_fraction; /* fractional time change since previous snapshot */
            delta_t_fraction = last_save_time == 0.0 ? INFINITY : (state.t - last_save_time) / last_save_time;
            save_snapshot_now = delta_t_fraction >= runtime_config.save_frequency_timing;
        } else if (runtime_config.has_save_frequency_density) {
            /* save based on central density change */
            double delta_rho_fraction = fabs((current_rho - last_save_rho) / last_save_rho); /* central density change */
            save_snapshot_now = delta_rho_fraction >= runtime_config.save_frequency_density;
        }

        if (save_snapshot_now) {
            if (gsidm_save_halo(&state, GSIDM_PREFIX_DEFAULT) != 0) {
                rc = 1;
                break;
            }
            last_save_time = state.t;
            last_save_rho = current_rho;
        }
    }

    gsidm_destroy_state(&state);
    return rc;
}

int main(int argc, char **argv) {
    gsidm_config_t config;
    int parse_rc;
    gsidm_set_default_config(&config);

    parse_rc = gsidm_parse_args(&config, argc, argv);
    if (parse_rc == 1) {
        return 0;
    }
    if (parse_rc != 0) {
        return 1;
    }

    return gsidm_run(&config);
}
