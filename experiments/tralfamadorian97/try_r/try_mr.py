import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

token="insert_token"
def go():
    # 1. Verify R is pointing to the Pixi environment
    # This should print a path inside your .pixi folder
    print(f"Using R at: {robjects.r['R.home']()[0]}")

    # 2. Import TwoSampleMR
    # This suppresses the R startup messages in your Python terminal
    try:
        twosamplemr = importr('TwoSampleMR')
        print("✅ TwoSampleMR loaded successfully!")
    except Exception as e:
        print(f"❌ Failed to load TwoSampleMR: {e}")

    print("getting outcomes")
    ao = twosamplemr.available_outcomes(
        opengwas_jwt=token
    )
    print(ao)
    print("getting exposure")
    exposure_dat = twosamplemr.extract_instruments("ieu-a-2",
                                                   opengwas_jwt=token
                                                   )
    print("getting outcome")
    # outcome_dat = twosamplemr.rextract_outcome_data(snps=exposure_dat$SNP, outcomes = "ieu-a-7",
    # opengwas_jwt = token
    # )


def tutorial_example_gemini():
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.robjects.packages import importr

    # Import the R packages
    tsmr = importr('TwoSampleMR')
    base = importr('base')

    print("--- Starting TwoSampleMR Analysis (Modern Context Manager) ---")

    # Define the converter we want to use (Default R rules + Pandas rules)
    # We will reuse this context in every step where we need conversion.
    conv = ro.default_converter + pandas2ri.converter

    # --- Step 1: List available GWASs ---
    # We use the context manager to automatically convert the R return value to a Pandas DataFrame
    with localconverter(conv):
        print("Fetching available outcomes...")
        ao = tsmr.available_outcomes(

            opengwas_jwt=token
        )
        # 'ao' is immediately a pandas DataFrame here
        print(ao.head(3))

    # --- Step 2: Get instruments ---
    with localconverter(conv):
        print("\nExtracting instruments for 'ieu-a-2' (BMI)...")
        exposure_dat = tsmr.extract_instruments(outcomes="ieu-a-2",

                                                opengwas_jwt=token
                                                )

    # --- Step 3: Get effects of instruments on outcome ---
    print("Extracting outcome data for 'ieu-a-7' (Coronary Heart Disease)...")

    # CRITICAL: We need to pass the SNP column (Pandas Series) to R.
    # While the converter *can* handle Series, explicit conversion is often safer
    # for function arguments to avoid ambiguity.
    snps_vector = ro.vectors.StrVector(exposure_dat['SNP'])

    with localconverter(conv):
        outcome_dat = tsmr.extract_outcome_data(
            snps=snps_vector,
            outcomes="ieu-a-7",
            opengwas_jwt=token
        )

    # --- Step 4: Harmonise the exposure and outcome data ---
    # Note: Since 'exposure_dat' and 'outcome_dat' are already Pandas DataFrames
    # (because we converted them on output previously), we need the converter active
    # so rpy2 knows how to turn them BACK into R DataFrames for the function input.
    with localconverter(conv):
        print("Harmonising data...")
        dat = tsmr.harmonise_data(exposure_dat, outcome_dat)

    # --- Step 5: Perform MR ---
    with localconverter(conv):
        print("Performing MR analysis...")
        res = tsmr.mr(dat)

        print("\n--- MR Results ---")
        print(res)

    # Optional: Save results to CSV using standard Pandas
    res.to_csv("mr_results.csv", index=False)
    print("\nResults saved to mr_results.csv")


if __name__ == "__main__":
    tutorial_example_gemini()