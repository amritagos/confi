import pandas as pd
from pathlib import Path


def main(
    output_csv: Path,
    molality: list,
    total_cation: list,
    total_anion: list,
    total_monomer: list,
    total_dimer: list,
    total_water: list,
    total_oct_water_shell: list,
    run_number: int,
):

    rows = []

    seed = 0  # The seed needed for PACKMOL. Each configuration should have a unique one
    for (
        i_molal,
        n_cations_tot,
        n_anions_tot,
        n_wat_tot,
        n_monomer,
        n_dimer,
        n_oct_water_shell,
    ) in zip(
        molality,
        total_cation,
        total_anion,
        total_water,
        total_monomer,
        total_dimer,
        total_oct_water_shell,
    ):
        if n_oct_water_shell is None:
            n_oct_water_shell = 0

        # Calculate the number of free cations and anions, and the number of free water molecules
        num_free_cation = n_cations_tot - n_oct_water_shell - n_monomer - n_dimer
        num_free_wat = n_wat_tot - n_oct_water_shell * 6
        num_free_anion = n_anions_tot - n_monomer - 2 * n_dimer

        for i_run in range(run_number):
            sample_id = f"{i_molal}_m_{n_monomer}_mono_{n_dimer}_dimer_{i_run}_run_num"  # Sample name
            seed += 1  # The seed for PACKMOL. Each run should have a unique one
            rows.append(
                {
                    "sample_name": sample_id,
                    "molality": i_molal,
                    "n_monomer": n_monomer,
                    "n_cations": num_free_cation,
                    "n_anions": num_free_anion,
                    "n_wat": num_free_wat,
                    "n_cation_water_oct_shell": n_oct_water_shell,
                    "n_dimer": n_dimer,
                    "seed": seed,
                    "run_number": i_run,
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    print(f"Written samples to file {output_csv.resolve()}")


if __name__ == "__main__":
    # Create the concentration
    molality = [4]
    total_cation = [320]
    total_anion = [int(2 * cat) for cat in total_cation]
    total_monomer = [0]
    total_dimer = [320]
    total_water = [4440 for m in molality]
    total_oct_water_shell = [0 for m in molality]
    run_number = 5  # Number of runs to perform per system
    output_csv = Path("config/samples.csv")  # samples will be written in tabular form

    main(
        output_csv,
        molality,
        total_cation,
        total_anion,
        total_monomer,
        total_dimer,
        total_water,
        total_oct_water_shell,
        run_number,
    )
