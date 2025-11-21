from corems.transient.input.brukerSolarix import ReadBrukerSolarix
from matplotlib import pyplot
from pathlib import Path

# file_path= 'tests/tests_data/ftms/ESI_NEG_SRFA.d'
# file_path= '/Users/anya/Coding/CoreMS/tests/tests_data/ftms/ESI_NEG_SRFA.d'
# When running in Docker with "-v /Users/anya/Coding/CoreMS/tests/tests_data:/testdata",
# the CoreMS test data will be available under /testdata inside the container.

# file_path = "/testdata/ftms/ESI_NEG_SRFA.d" # works_
# file_path = '/app/data/Measurement_d-files/noncalibrated_d-files/ESI_neg_G015548-0_100_200sc_000001.d/' # doesnt work


DATA_FOLDER = Path('/app/data')
RESULTS_FOLDER = Path('/app/results/')

# file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "ESI_neg_G015548-0_100_200sc_000001.d/" # doesnt work

file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "ESI_neg_G017736-0_100_200sc_000001.d" #works
# file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "Esi_neg_G000427-0_100_200scTol_000001.d"# doesnt work
# file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "Esi_neg_G000420-0_100_200scTol_000001.d" # doesnt work
# file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "ESI_neg_G015548-0_100_200sc_000001.d"# doesnt work
# file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "ESI_neg_G017736-0_100_200sc_000001.d" # works
# file_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "Esi_neg_G000427-0_100_200scTol_000001.d" # doesnt work
# ffile_path = DATA_FOLDER / "Measurement_d-files" / "calibrated_d-files" / "Esi_neg_G000420-0_100_200sc_000001.d" # doesnt work
# Instatiate the Bruker Solarix reader with the filepath
file_path = str(file_path)  # if some library requires string
bruker_reader = ReadBrukerSolarix(file_path)

# Use the reader to instatiate a transient object
bruker_transient_obj = bruker_reader.get_transient()

# Calculate the transient duration time
T =  bruker_transient_obj.transient_time

# Use the transient object to instatitate a mass spectrum object
mass_spectrum_obj = bruker_transient_obj.get_mass_spectrum(plot_result=False, auto_process=True)

# NOTE: Molecular formula search is disabled in this test script to avoid
# database configuration (PostgreSQL / SQLite) issues inside Docker.
# If you want to enable it later, uncomment the block below and make sure
# the CoreMS molecular formula database is configured.
#
# from corems.molecular_id.search.molecularFormulaSearch import SearchMolecularFormulas
# SearchMolecularFormulas(mass_spectrum_obj, first_hit=False).run_worker_mass_spectrum()
#
# # Iterate over mass spectral peaks objs within the mass_spectrum_obj
# for mspeak in mass_spectrum_obj.sort_by_abundance():
#
#     # If there is at least one molecular formula associated, mspeak returns True
#     if mspeak:
#
#         # Get the molecular formula with the highest mass accuracy
#         molecular_formula = mspeak.molecular_formula_lowest_error
#
#         # Plot mz and peak height
#         pyplot.plot(mspeak.mz_exp, mspeak.abundance, 'o', c='g')
#
#         # Iterate over all molecular formulas associated with the ms peaks obj
#         for molecular_formula in mspeak:
#
#             # Check if the molecular formula is a isotopologue
#             if molecular_formula.is_isotopologue:
#
#                 # Access the molecular formula text representation and print
#                 print(molecular_formula.string)
#
#                 # Get 13C atoms count
#                 print(molecular_formula['13C'])
#     else:
#         # Get mz and peak height
#         print(mspeak.mz_exp, mspeak.abundance)

# Save data
## to a csv file
mass_spectrum_obj.to_csv(str(RESULTS_FOLDER / "filename"))
mass_spectrum_obj.to_hdf(str(RESULTS_FOLDER / "filename"))
# to pandas Datarame pickle
mass_spectrum_obj.to_pandas(str(RESULTS_FOLDER / "filename"))

# Extract data as a pandas Dataframe
df = mass_spectrum_obj.to_dataframe()