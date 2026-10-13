import numpy as np, src.fields.reference_field_grid as R
def _cell_centers(vmin,vmax,period_m):
    span=vmax-vmin; n=max(1,int(np.floor(span/period_m+1e-6))); c=0.5*(vmin+vmax)
    return c+(np.arange(n)-(n-1)/2.0)*period_m, n
R._cell_centers=_cell_centers
