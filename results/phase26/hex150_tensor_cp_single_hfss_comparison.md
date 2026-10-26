# hex150_tensor_cp_single: HFSS vs prediction (pin model: [user to state: scripted pin or connector])

| Quantity | HFSS | Predicted | Target / falsifier | State |
|---|---|---|---|---|
| Peak RHCP realized gain [dBi] | 23.03 | 23.5 (bare) .. 25.2 (loaded) | >= 24; falsified < 19 | GAIN_WITHIN_PREDICTED_RANGE |
| Peak direction (theta, phi) [deg] | 30, 225 | 30, 224 | within 2 deg of 30, 225 | OK (0.0 deg) |
| XPD at peak [dB] | 23.0 | 16.9 | >= 25 | delta +6.1 dB |
| AR at peak / worst in -3 dB beam [dB] | 1.23 / 1.58 | 2.50 / 3.07 | <= 3 | |
| RHCP / total power within 8 deg | 0.994 | 0.979 | falsified < 0.8 | OK |
| RHCP SLL [dB] | -18.5 | -22.8 | <= -15 | |

## Power budget (fractions of the radiated power; HFSS radiation efficiency ~ 1 so ~ share of accepted)
| Region | HFSS | Model |
|---|---|---|
| Main beam <= 8 deg (RHCP) | 0.271 (0.272 total) | 0.302 |
| Near sidelobes 8-20 deg | 0.024 | in aperture pattern only |
| Far sidelobes (theta <= 70) | 0.171 | in aperture pattern only |
| Grazing 70-90 deg | 0.264 | 0.681 (upper bound incl. back, not modelled) |
| Back hemisphere | 0.269 | (in the above) |
| LHCP, all angles | 0.310 | 0.008 (aperture only) |

- s_parameter: `{"s_db_at_design": -15.002165023156, "f_nearest_ghz": 10.54, "min_s_db": -27.2692214220148, "f_at_min_ghz": 10.08, "band_below_minus10_db_ghz": [8.81, 11.59], "design_inside_minus10_db": true, "prediction": "NO_PREDICTION (match not modelled; the pin is the scripted/connector feed)", "brief_target": "S11 <= -10 dB at 10.54 GHz (launcher redesign constraint)"}`
- ar_sweep: `{"phi_cut_deg": 224.0, "frequencies_ghz": [8.0, 8.125, 8.25, 8.375, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0, 10.25, 10.5, 11.0, 11.25, 11.5, 11.625, 11.75, 12.0], "ar_at_beam_theta_db": [18.2156305400679, 9.65541294149227, 1.2250334292437, 9.45647034138737, 17.0656584917165, 22.3711341874021, 12.8159619852049, 5.98265632236716, 1.54544599119052, 1.51705259964405, 1.92351385034272, 1.6962766563518, 1.30453051322988, 1.27541688740696, 1.33098355254461, 1.10756734492841, 3.08129157609375, 2.46710394160797, 2.49748377694982], "frequency_used_ghz": 10.5, "ar_at_beam_theta_at_used_freq_db": 1.30453051322988, "min_ar_in_cut_db": 0.248517525594402, "min_ar_theta_deg": 26.0, "freqs_with_ar_le_3db_at_beam_theta_ghz": [8.25, 9.5, 9.75, 10.0, 10.25, 10.5, 11.0, 11.25, 11.5, 11.75, 12.0], "contiguous_3db_band_around_design_ghz": [9.5, 11.5]}`
- ar_cut: `"SEE ar_sweep"`
- surface_wave_decay: `"NOT_PROVIDED"`
- antenna_params_selected: `{"Total Efficiency": 0.975038, "Radiation Efficiency": 1.003544, "Accepted Power": 971.595098, "Radiated Power": 975.038402, "Incident Power": 1.0}`

The model gives no pattern to the pin's direct space wave (0.41) or the unconverted remainder (0.27): they are lumped into the grazing+back upper bound. HFSS LHCP includes the pin's direct radiation, the model's LHCP does not.
