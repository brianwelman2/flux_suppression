# Script with functions to generate params

from kalcal.calibration import vanilla
import numpy as np
import optuna
from optuna import visualization as vis
import contextlib, os
from casacore.tables import table
from datetime import datetime


def rmse(A, B):
    res = A - B
    return np.sqrt(np.vdot(res, res).real/res.size)


def optuna_run(lb, ub, percent, n_trials, config, antenna, src):
    msname = f"ms/{antenna}.ms"
    gains_path = f"gains/true/{antenna}/{antenna}-src-{src}-gains.npy"
    true_gains = np.load(gains_path)[..., (0, 3)]
    results = np.zeros((2, n_trials), dtype=np.float64)

    def objective(trial):
        options = config.copy()
        sigma_f = trial.suggest_float("sigma_f", lb, ub)
        options["sigma_f"] = sigma_f

        with open(os.devnull, 'w') as devnull:
            with contextlib.redirect_stdout(devnull):
                vanilla.calibrate(msname, **options)
        
        smooth_gains = np.load(options["out_smoother"])[..., (0, 3), 0]
        rmse_result = rmse(true_gains, smooth_gains)
        
        results[0, trial.number] = options["sigma_f"]
        results[1, trial.number] = rmse_result

        return rmse_result

    study = optuna.create_study(study_name=f"ANT={antenna} with SRC={src} at FP={percent}", sampler=optuna.samplers.TPESampler())

    study.optimize(objective, n_trials=n_trials, catch=(ZeroDivisionError,))
    best_params = study.best_params
    indices = results[1].argsort()

    top_25 = results[:, indices[0:25]]
    rlb = top_25[0].min()
    rub = top_25[0].max()
    
    print("Best Sigma_f:", best_params)
    print(f"Top 25: Min={rlb}, Max={rub}")

    return best_params["sigma_f"], results[1].min(), rlb, rub 


if __name__ == "__main__":

    # Calibration Config
    config = {
        "filter"        : 1,
        "smoother"      : 1,
        "algorithm"     : "NUMBA",
        "sigma_f"       : 0.01,
        "sigma_n"       : 1.4142135623730951,
        "step_control"  : 0.5,
        "model_column"  : f"",
        "vis_column"    : f"",
        "weight_column" : "WEIGHT",
        "out_filter"    : f"gains/kalcal/other/filter.npy",
        "out_smoother"  : f"gains/kalcal/other/smoother.npy",
        "out_data"      : None,
        "out_weight"    : None,
        "ncpu"          : 24,
        "yaml"          : None
    }

    setup = {
        "kat7": (1, 2),
        "meerkat": (1, 2, 100),
        "vlab": (1, 2, 100)
    }

    percents = [30, 50, 70, 100]
    
    with open("logs/optuna_grid.csv", "w") as file:
        file.write("time, antenna, sources, flux_percentage, sigma_f, rmse, lowerbound, upperbound\n")

    for antenna, srcs in setup.items():
        for src in srcs:
            if src == 1:
                config["model_column"] = "SRC1_MODEL"
                config["vis_column"] = "SRC1_DATA"
                best_sigma_f, best_rmse, lowerbound, upperbound = optuna_run(0.00001, 0.03, 100, 50, config, antenna, src)
                with open("logs/optuna_grid.csv", "a") as file:
                    file.write(f"{datetime.now()}, {antenna}, {src}, 100, {best_sigma_f}, {best_rmse}, {lowerbound}, {upperbound}\n")

            elif src == 2:
                config["model_column"] = "SRC2_MODEL"
                config["vis_column"] = "SRC2_DATA"
                best_sigma_f, best_rmse, lowerbound, upperbound = optuna_run(0.00001, 0.03, 100, 50, config, antenna, src)
                with open("logs/optuna_grid.csv", "a") as file:
                    file.write(f"{datetime.now()}, {antenna}, {src}, 100, {best_sigma_f}, {best_rmse}, {lowerbound}, {upperbound}\n")

                config["model_column"] = "SRC1_MODEL"
                config["vis_column"] = "SRC2_DATA"
                best_sigma_f, best_rmse, lowerbound, upperbound = optuna_run(0.00001, 0.03, 70, 50, config, antenna, src)
                with open("logs/optuna_grid.csv", "a") as file:
                    file.write(f"{datetime.now()}, {antenna}, {src}, 70, {best_sigma_f}, {best_rmse}, {lowerbound}, {upperbound}\n")
            
            elif src == 100:
                for percent in percents:
                    config["model_column"] = f"FP{percent}_SRC100_MODEL"
                    config["vis_column"] = "FP100_SRC100_DATA"
                    best_sigma_f, best_rmse, lowerbound, upperbound = optuna_run(0.00001, 0.03, percent, 50, config, antenna, src)
                    with open("logs/optuna_grid.csv", "a") as file:
                        file.write(f"{datetime.now()}, {antenna}, {src}, {percent}, {best_sigma_f}, {best_rmse}, {lowerbound}, {upperbound}\n")