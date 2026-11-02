# Rogers RO3010 - datasheet values used (checked 2026-11-02 against the Rogers product page, rogerscorp.com RO3010 laminates)

| Quantity | Datasheet | Used in the model |
|---|---|---|
| Dk, process, 10 GHz | 10.20 +/- 0.30 (typical) | eps_r = 10.2 (isotropic) |
| Dk, design | 11.2 (typical) | NOT used; sensitivity run requested (env `P30_EPS=11.2`) |
| Dissipation factor, 10 GHz | 0.0022 (typical) | tan d = 0.0022 |
| Thickness | 0.635 mm (25 mil) standard; tolerance NOT stated on the product page | 0.635 mm |
| CTE x/y/z | 13 / 11 / 16 ppm/C | not used |
| Copper | 18 um (0.5 oz) per the request; the page only mentions 1 oz (35 um) peel-strength tests | sigma 5.8e7 S/m, smooth, Rs = 25.6 mohm/sq at 10.54 GHz |

Open points (NOT_AVAILABLE_FROM_CURRENT_DATA): the page does not say whether 10.2 is the z-axis value, nor why design Dk is
1 higher; the Dk test method is not stated on the page. The surface-wave cells are very close to a resonance, so +-0.3 in Dk
and the 10.2-vs-11.2 ambiguity are a first-order uncertainty of every predicted X_bar/beta/alpha in this repository.
