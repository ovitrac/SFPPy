"""
survey — Production Survey Module for SFPPy
============================================

Survey-scale exposure estimation with persistent caching,
parallel master curve computation, and dynamic PDF generation.

@project: SFPPy/INSERM — Survey-scale exposure estimation
@author: Olivier Vitrac, PhD, HDR
@email: olivier.vitrac@gmail.com
@license: MIT
"""

from survey.survey import Survey
from survey.models import (
    LayerSpec,
    PackagingSpec,
    SubstanceSpec,
    PriorSpec,
    SurveyConfig,
)
from survey.visualization import (
    plot_pdf_cdf,
    save_figure,
    generate_survey_plots,
    HAS_MATPLOTLIB,
)
from survey.tables import (
    SurveyResult,
    extract_survey_result,
    format_inputs_table,
    format_outputs_table,
    generate_summary_report,
    save_summary_report,
)
from survey.aggregation import (
    combine_tensors,
    combine_multiple_tensors,
    aggregate_components,
    aggregate_family_weighted,
    compute_pdf_from_samples,
    compute_statistics,
    compute_percentiles,
    compute_risk_percentiles,
    RISK_PERCENTILES,
    aggregate_packaging_components,
    aggregate_food_exposure,
    quantile_from_samples,
    quantile_from_cdf,
)

__version__ = "0.4.0"
__all__ = [
    "Survey",
    "LayerSpec",
    "PackagingSpec",
    "SubstanceSpec",
    "PriorSpec",
    "SurveyConfig",
    # Visualization
    "plot_pdf_cdf",
    "save_figure",
    "generate_survey_plots",
    "HAS_MATPLOTLIB",
    # Tables
    "SurveyResult",
    "extract_survey_result",
    "format_inputs_table",
    "format_outputs_table",
    "generate_summary_report",
    "save_summary_report",
    # Aggregation
    "combine_tensors",
    "combine_multiple_tensors",
    "aggregate_components",
    "aggregate_family_weighted",
    "compute_pdf_from_samples",
    "compute_statistics",
    "compute_percentiles",
    "compute_risk_percentiles",
    "RISK_PERCENTILES",
    "aggregate_packaging_components",
    "aggregate_food_exposure",
    "quantile_from_samples",
    "quantile_from_cdf",
]
