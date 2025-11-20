import pandas as pd
from itertools import product
from main import SAMPLE_TYPE

sample_type=SAMPLE_TYPE

if sample_type == 'crude_oil': 
    # (1) Crude Oil

    # Hydrogen-to-Carbon Ratio Constraints
    H_C_MIN = 0.2
    H_C_MAX = 2.2

    # Elemental Ratios
    O_C_MIN = 0.0
    O_C_MAX = 0.6
    N_C_MIN = 0.0
    N_C_MAX = 0.2
    S_C_MIN = 0.0
    S_C_MAX = 0.2

    # CHNOS Limits including *C (13C) and Na
    CHNOS_LIMITS = {
        'C': (0, 100),
        'H': (0, 200),
        'N': (0, 4),
        'O': (0, 8),
        'S': (0, 4),
        '*C': (0, 2),
        'Na': (0, 1),
        '34S': (0, 1)  # Sulfur isotope
    }
elif sample_type == 'sedimentary_rock':
    # (2) Sedimentary Rock (Solvent Extract/Bitumen)

    # Hydrogen-to-Carbon Ratio Constraints
    H_C_MIN = 0.2
    H_C_MAX = 2.2

    # Elemental Ratios
    O_C_MIN = 0.0
    O_C_MAX = 0.6
    N_C_MIN = 0.0
    N_C_MAX = 0.3
    S_C_MIN = 0.0
    S_C_MAX = 0.3

    # CHNOS Limits including *C (13C), Na, Cl, and Fe
    CHNOS_LIMITS = {
        'C': (0, 100),
        'H': (0, 200),
        'N': (0, 4),
        'O': (0, 6),
        'S': (0, 4),
        '*C': (0, 2),
        'Na': (0, 1),
        '34S': (0, 1),
        '35Cl': (0, 1),
        '37Cl': (0, 1),
        '54Fe': (0, 1),
        '56Fe': (0, 1),
        '57Fe': (0, 1),
        '58Fe': (0, 1)
    }
elif sample_type == 'coal':

    # (3) Coal (Solvent Extract/Bitumen)

    # Hydrogen-to-Carbon Ratio Constraints
    H_C_MIN = 0.2
    H_C_MAX = 2.2

    # Elemental Ratios
    O_C_MIN = 0.0
    O_C_MAX = 0.6
    N_C_MIN = 0.0
    N_C_MAX = 0.3
    S_C_MIN = 0.0
    S_C_MAX = 0.3

    # CHNOS Limits including *C (13C) and Na
    CHNOS_LIMITS = {
        'C': (0, 100),
        'H': (0, 200),
        'N': (0, 2),
        'O': (0, 12),
        'S': (0, 2),
        '*C': (0, 2),
        'Na': (0, 1)
    }
elif sample_type == 'natural_water':

    # (4) Natural Water DOM Sample (After SPE)

    # Hydrogen-to-Carbon Ratio Constraints
    H_C_MIN = 0.2
    H_C_MAX = 2.2

    # Elemental Ratios
    O_C_MIN = 0.1
    O_C_MAX = 1.0
    N_C_MIN = 0.0
    N_C_MAX = 0.1
    S_C_MIN = 0.0
    S_C_MAX = 0.1

    # CHNOS Limits including *C (13C), Na, and P
    CHNOS_LIMITS = {
        'C': (0, 80),
        'H': (0, 160),
        'N': (0, 3),
        'O': (0, 35),
        'S': (0, 2),
        '*C': (0, 2),
        'Na': (0, 1),
        'P': (0, 1)
    }
else:
    # Default sample type (e.g., ESI-FT-ICR-MS) with [M-H]- ionization

    # Hydrogen-to-Carbon Ratio Constraints (adjusted for [M-H]- ions)
    H_C_MIN = 0.4
    H_C_MAX = 2.8

    # CHNOS limits including *C (13C) and Na
    CHNOS_LIMITS = {
        'C': (7, 120),
        'H': (0, 220),  # Hydrogen total before subtracting 1 for [M-H]-
        'N': (0, 3),
        'O': (0, 8),
        'S': (0, 2),
        '*C': (0, 2),
        'Na': (0, 2)  # Consider Na for adducts
    }

    # ------------ FILTER CHECK FUNCTION (H/C ratio + Rule 5 + DBE) ---------------
    def formula_filter(c, h, o, n, s, na):
        if c == 0:  # Avoid division by zero
            return False

        # Hydrogen for neutral molecule (adjusted for [M-H]- ion)
        h_neutral = h + 1  
        if h_neutral < 0:
            return False

        # H/C ratio for neutral molecule
        h_c_ratio = h_neutral / c

        # Other ratios for Rule 5
        o_c = o / c
        n_c = n / c
        s_c = s / c

        # Apply Rule 5 constraints
        if not ((H_C_MIN <= h_c_ratio <= H_C_MAX) and
                (0 <= o_c <= 3) and
                (0 <= n_c <= 1.3) and
                (0 <= s_c <= 0.8)):
            return False

        # ---------- DBE calculation ----------
        if na == 0:
            # Regular [M-H]- case
            DBE = 1 + (2 * c + 2 + n - h) / 2
        else:
            # [M + Na - 2H]- case
            DBE = 1 + (2 * c + 2 + n - h) / 2

        if DBE < 0 or DBE % 1 != 0:  # DBE must be a non-negative integer
            return False

        return True  # Passes all filters


# ------------ FORMULA GENERATION -----------
def generate_formulas(sample_type):
    
    
    formulas = []
    for c, h, n, o, s, c13, na in product(
        range(*CHNOS_LIMITS['C']),
        range(*CHNOS_LIMITS['H']),
        range(*CHNOS_LIMITS['N']),
        range(*CHNOS_LIMITS['O']),
        range(*CHNOS_LIMITS['S']),
        range(*CHNOS_LIMITS['*C']),
        range(*CHNOS_LIMITS['Na'])  # Include Na adducts
        ):
        # Validity checks
        if c == 0 or (c13 > c):  # At least one 12C and valid *C
            continue

        if not formula_filter(c, h, o, n, s, na):  # Apply filters
            continue

        # ---------- Mass for [M-H]- or [M + Na - 2H]- ion ----------
        neutral_mass = (
            c * 12.000000000 + 
            c13 * 13.0033548378 + 
            h * 1.007825032 + 
            n * 14.003074005 + 
            o * 15.994914620 + 
            s * 31.972071174 +
            na * 22.98977
        )

        # Ion mass based on presence of Na
        if na == 0:
            ion_mass = neutral_mass - 1.007825032 + 0.000548579909 # [M-H]-
        else:
            ion_mass = neutral_mass - 2 * 1.007825032 + 0.000548579909 # [M + Na - 2H]-

        # ---------- Assemble formula string for ion ----------
        h_deprotonated = h - 1  # For [M-H]- or [M + Na - 2H]-
        formula = f"C{c}" + (f"*C{c13}" if c13 > 0 else "") + \
                  (f"H{h_deprotonated}" if h_deprotonated > 0 else "") + \
                  (f"N{n}" if n > 0 else "") + \
                  (f"O{o}" if o > 0 else "") + \
                  (f"S{s}" if s > 0 else "") + \
                  (f"Na{na}" if na > 0 else "")

        # Collect valid formula entry
        formulas.append((ion_mass, formula, c, h, n, o, s, c13, na))

    # Return as DataFrame
    return pd.DataFrame(formulas, columns=['Mass', 'Formula', 'C', 'H', 'N', 'O', 'S', '*C', 'Na'])