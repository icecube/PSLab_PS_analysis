#ifndef FLUXUSDEV_PLOTFLUX_READ_FLUX_CURVES_H_
#define FLUXUSDEV_PLOTFLUX_READ_FLUX_CURVES_H_

// Forward Declarations (when feasible, more efficient than including headers)
class FluxPlotManager;


void read_flux_curves(char* filename, FluxPlotManager& fpm,
		      FluxPlotManager& fpmEvents,
		      int color=1, int width=1, int style=1,
		      bool optConnect = false, int connectPoint = 0);


#endif // FLUXUSDEV_PLOTFLUX_READ_FLUX_CURVES_H_
