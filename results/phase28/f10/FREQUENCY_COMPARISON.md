# Phase 28 item 2: 10.0 GHz against 10.54 GHz

Generated from `frequency_comparison.json` only. Both sweeps use the same code (`src/hologram/coded_feasibility.py`) and the same cases; they differ in the frequency, the supercell table (rebuilt at the new frequency with the same method)
and, as an ASSUMPTION carried over unchanged, the launch fraction (the 10.54 GHz Phase 26 value; it shifts every absolute gain equally and no loss). `difference` = 10.0 GHz value minus 10.54 GHz value.

| quantity | 10.54 GHz | 10.0 GHz | difference |
|---|---|---|---|
| frequency_ghz | 10.54 | 10 | -0.54 |
| table_alpha_max_np_per_m | 4.259 | 4.123 | -0.1359 |
| table_alpha_min_np_per_m | 0.4216 | 0.1127 | -0.3089 |
| supported_scan_theta_deg | 20 | 50 | 30 |
| design_point_eta_conv | 0.5384 | 0.504 | -0.03438 |
| design_point_outside_table_fraction | 0.223 | 0 | -0.223 |
| design_point_cp_continuous_gain_dbi | 23.6 | 22.85 | -0.7463 |
| design_point_cp_1bit_loss_db | 3.819 | 3.896 | 0.07754 |
| design_point_cp_2bit_loss_db | 0.9456 | 0.9009 | -0.04468 |
| design_point_cp_3bit_loss_db | 0.224 | 0.231 | 0.007056 |
| alpha_cap_for_1db_np_per_m | 2.308 | 2.089 | -0.2194 |
| tensor_delta_ms | 2.296 | 2.444 | 0.148 |
| tensor_dx_aniso_ohm | 92.33 | 82.6 | -9.731 |
| mean_beta_over_k0 | 1.062 | 1.053 | -0.00943 |
| mean_loading_ms | 7.559 | 7.483 | -0.07535 |
| combined_1db_alpha_cap_np_per_m | 4.01 | 3.661 | -0.3485 |
| combined_1db_dx_aniso_ohm | 132 | 118.4 | -13.66 |
| combined_1db_delta_ms | 3.027 | 3.236 | 0.2095 |
| combined_1.5db_alpha_cap_np_per_m | 2.737 | 2.443 | -0.2946 |
| combined_1.5db_dx_aniso_ohm | 102.6 | 90.88 | -11.69 |
| combined_1.5db_delta_ms | 2.5 | 2.643 | 0.1426 |
| combined_2db_alpha_cap_np_per_m | 2.263 | 2.017 | -0.246 |
| combined_2db_dx_aniso_ohm | 91.23 | 80.88 | -10.35 |
| combined_2db_delta_ms | 2.273 | 2.402 | 0.1282 |
| cp_2bit_loss_db_pitch_4 | 0.9102 | 0.911 | 0.0007914 |
| cp_2bit_loss_db_pitch_5 | 0.9206 | 0.916 | -0.0046 |
| cp_2bit_loss_db_pitch_6.4 | 0.9456 | 0.9009 | -0.04468 |
| cp_2bit_loss_db_pitch_8 | 0.9036 | 0.9085 | 0.004936 |
| outside_table_fraction_theta_0 | 0 | 0 | 0 |
| cp_2bit_gain_dbi_theta_0 | 22.67 | 22.43 | -0.2428 |
| outside_table_fraction_theta_10 | 0 | 0 | 0 |
| cp_2bit_gain_dbi_theta_10 | 22.75 | 22.12 | -0.6247 |
| outside_table_fraction_theta_20 | 0 | 0 | 0 |
| cp_2bit_gain_dbi_theta_20 | 22.95 | 22.06 | -0.8929 |
| outside_table_fraction_theta_30 | 0.223 | 0 | -0.223 |
| cp_2bit_gain_dbi_theta_30 | 22.65 | 21.95 | -0.7016 |
| outside_table_fraction_theta_40 | 0.4809 | 0 | -0.4809 |
| cp_2bit_gain_dbi_theta_40 | 22.08 | 21.35 | -0.7338 |
| outside_table_fraction_theta_50 | 0.5719 | 0.01156 | -0.5604 |
| cp_2bit_gain_dbi_theta_50 | 21.07 | 20.37 | -0.705 |
