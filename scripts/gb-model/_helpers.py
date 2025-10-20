# SPDX-FileCopyrightText: gb-open-market-model contributors
#
# SPDX-License-Identifier: MIT


import pandas as pd


def _strip_str(series: pd.Series) -> pd.Series:
    """Strip whitespace from strings in a pandas Series."""
    return series.str.strip() if series.dtype == "object" else series
