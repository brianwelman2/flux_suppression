import os
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import Tigger


large_skymodels = [
    "model-fp-100-src-100.txt",
    "model-fp-70-src-100.txt",
    "model-fp-50-src-100.txt",
    "model-fp-30-src-100.txt"
]

arrays = ["kat7", "vlab", "meerkat"]

for ant in arrays:
    # Check small sky-models for correct percentage
    src_1_model = Tigger.load(f"skymodels/{ant}/model-src-1.txt", verbose=False)
    src_2_model = Tigger.load(f"skymodels/{ant}/model-src-2.txt", verbose=False)

    # Total flux for each model
    total_src_1 = sum(src.flux.I for src in src_1_model.sources)
    total_src_2 = sum(src.flux.I for src in src_2_model.sources)

    # Percentage for each model
    LHS, RHS = 100, round(total_src_1/total_src_1 * 100)
    print(f"[Name = '{ant}/model-src-1.txt', Sources = 1, Flux-Percentage = {LHS}%] => Sky-model = {RHS}%", u"\u2713" if LHS == RHS else u"\u2715")

    LHS, RHS = 100, round(total_src_2/total_src_2 * 100)
    print(f"[Name = '{ant}/model-src-2.txt', Sources = 2, Flux-Percentage = {LHS}%] => Sky-model = {RHS}%", u"\u2713" if LHS == RHS else u"\u2715")

    LHS, RHS = 70, round(total_src_1/total_src_2 * 100)
    print(f"[Name = '{ant}/model-src-2.txt', Sources = 2, Flux-Percentage = {LHS}%] => Sky-model = {RHS}%", u"\u2713" if LHS == RHS else u"\u2715")

    if ant == "kat7":
        continue
    
    src_100_model = Tigger.load(f"skymodels/{ant}/" + large_skymodels[0], verbose=False)
    total_src_100 = sum(src.flux.I for src in src_100_model.sources)

    for lsm in large_skymodels:        
        split_filename = lsm.split(".")[0].split("-")
        percent, src = int(split_filename[2]), int(split_filename[4])

        fp_src_100_model = Tigger.load(f"skymodels/{ant}/" + lsm, verbose=False)
        total_fp_src_100_model = sum(src.flux.I for src in fp_src_100_model.sources)

        LHS, RHS = percent, (total_fp_src_100_model/total_src_100 * 100)
        print(f"[Name = '{ant}/{lsm}', Sources = 100, Flux-Percentage = {LHS}%] => Sky-model = {RHS}%", u"\u2713" if LHS == RHS else u"\u2715")
