import sys
from pathlib import Path
import matplotlib.pyplot as plt

from schon.read_bruker_fid_corems import load_mass_spectrum_from_d_folder


def plot_spectrum_from_d(d_folder: str, out_png: str | None = None):
    d_path = Path(d_folder)
    ms_obj = load_mass_spectrum_from_d_folder(str(d_path))

    # m/z and abundances from the MassSpectrum (consistent with CoreMS)
    mz = [p.mz_exp for p in ms_obj]          # experimental m/z
    inten = [p.abundance for p in ms_obj]    # intensity

    plt.figure(figsize=(10, 6))
    plt.vlines(mz, [0] * len(inten), inten, linewidth=0.5)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    plt.title(d_path.name)
    plt.tight_layout()

    if out_png:
        plt.savefig(out_png, dpi=200)
        print(f"✅ Spectrum saved to {out_png}")
    else:
        plt.show()


def main():
    """
    Usage:
        python -m SCHON.corems_pipeline.visualize_spectrum              # uses default .d folder
        python -m SCHON.corems_pipeline.visualize_spectrum path/to.d    # custom .d folder
        python -m SCHON.corems_pipeline.visualize_spectrum path/to.d out.png
    """
    # Default .d folder inside the Docker container
    default_d_folder = (
        "/app/Measurement_d-files/calibrated_d-files/"
        "ESI_neg_G017736-0_100_200sc_000001.d"
    )

    if len(sys.argv) == 1:
        # No arguments: use the default .d folder, no PNG output
        d_folder = default_d_folder
        out_png = None
        print(f"Using default .d folder: {d_folder}")
    elif len(sys.argv) == 2:
        # One argument: custom .d folder, no PNG output
        d_folder = sys.argv[1]
        out_png = None
    else:
        # Two arguments: custom .d folder and output PNG path
        d_folder = sys.argv[1]
        out_png = sys.argv[2]

    if out_png is None:
        out_png = f"/app/SCHON/output_plots/{Path(d_folder).name}.png"

    plot_spectrum_from_d(d_folder, out_png)


if __name__ == "__main__":
    main()