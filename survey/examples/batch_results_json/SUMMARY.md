# SFPPy Batch Simulation Report

**Generated:** 2025-12-19T14:22:07.459850
**Duration:** 5.9s

## Summary

| Metric | Value |
|--------|-------|
| Total scenarios | 4 |
| Successful | 4 |
| Failed | 0 |

## Results

| Scenario | Polymer | Simulant | Thickness (µm) | Substances | Q95 (mg/kg) | Status |
|----------|---------|----------|----------------|------------|-------------|--------|
| cheese_PP_tray | PP | oliveoil | 400 | 3 | 1.6774 | ✅ |
| water_HDPE_cap | HDPE | water | 800 | 3 | 125.9409 | ✅ |
| water_PET_1500mL | gPET | water | 180 | 3 | 0.0513 | ✅ |
| yogurt_PET_125g | gPET | ethanol50 | 200 | 3 | 0.0177 | ✅ |

## Quantile Details

| Scenario | Q50 | Q75 | Q90 | Q95 | Q99 |
|----------|-----|-----|-----|-----|-----|
| cheese_PP_tray | 0.7695 | 1.1238 | 1.4781 | 1.6774 | 2.0427 |
| water_HDPE_cap | 59.8494 | 87.0204 | 112.7226 | 125.9409 | 146.5026 |
| water_PET_1500mL | 0.0212 | 0.0330 | 0.0444 | 0.0513 | 0.0631 |
| yogurt_PET_125g | 0.0074 | 0.0113 | 0.0153 | 0.0177 | 0.0217 |