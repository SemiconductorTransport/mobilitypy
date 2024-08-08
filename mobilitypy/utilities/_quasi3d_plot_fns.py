from matplotlib import colors, cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri 
import matplotlib.ticker as ticker
import numpy as np
from ._general_plot_functions import _GeneratePlots

## ============================================================================
class PlotQuasi3DFuns(_GeneratePlots):
    def __init__(self, save_figure_dir='.'):
        _GeneratePlots.__init__(self, save_figure_dir=save_figure_dir)    
    
    @classmethod
    def _triangulation(cls, x_values, y_values):
        return tri.Triangulation(x_values, y_values)
    
    @classmethod
    def _meshgrid(cls, x_values, y_values, npoints:int=20): 
        xi, yi = np.meshgrid(np.linspace(x_values.min(), x_values.max(), npoints), 
                             np.linspace(y_values.min(), y_values.max(), npoints))
        return xi, yi
    
    @classmethod                    
    def _linear_interpolation(cls, triangles_, x_values, y_values, z_values):
        interp_lin = tri.LinearTriInterpolator(triangles_, z_values)
        return interp_lin(x_values, y_values)
    
    def InterPolation(self, x_values, y_values, z_values, method:str='linear',
                      interpolation_points:int=20):
        assert method in ['linear'], 'Requested interpolation method not implemented yet.'
        triang = self._triangulation(x_values, y_values)
        xi, yi = self._meshgrid(x_values, y_values, npoints=interpolation_points)
        zi = None
        if method == 'linear':
            zi = self._linear_interpolation(triang, xi, yi, z_values)
        return xi, yi, zi
    
    def CreateColorbarMapableObject(self, vmin=None, vmax=None, colorbar_scale='normal',
                                    color_map='viridis'):
        if colorbar_scale=='normal':
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
        elif colorbar_scale=='log':
            norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        mappable = cm.ScalarMappable(norm=norm, cmap=color_map)
        return norm, mappable
    
    def _PlotContour(self, xi, yi, zi, fig=None, axs=None, x_label:str='', y_label:str='',
                    title_label:str=None, z_label:str='', 
                    tick_multiplicator:list=[None, None, None, None],
                    vmin=None, vmax=None, cbar_mappable=None, norm=None,
                    color_map='viridis', show_contour_lines:bool=False, 
                    cbar_text:str=None, show_colorbar:bool=False):
        self.fig, axs = self._set_figure(fig=fig, axs=axs)

        CS = axs.contourf(xi, yi, zi, cmap=color_map, norm=norm)
        
        if show_contour_lines: 
            CS2 = axs.contour(CS, levels=CS.levels, colors='k')
        if show_colorbar:
            cbar=self._set_colorbar(axs, self.fig, CS=CS, cbar_mappable=cbar_mappable, 
                                    cbar_text=cbar_text)        
        self._set_labels(axs, x_label=x_label, y_label=y_label, title_label=title_label)
        self._set_tickers(axs, tick_multiplicator=tick_multiplicator)
        return self.fig, axs

    def _PlotScatter(self, xi, yi, zi, fig=None, axs=None, x_label:str='', y_label:str='',
                    title_label:str=None, z_label:str='', show_colorbar:bool=False,
                    tick_multiplicator:list=[None, None, None, None],
                    marker_size=None, vmin=None, vmax=None, cbar_mappable=None, norm=None,
                    color_map='viridis', cbar_text:str=None, marker='o'):

        self.fig, axs = self._set_figure(fig=fig, axs=axs)
            
        CS = axs.scatter(xi, yi, c=zi, cmap= color_map, norm=norm, 
                         marker=marker, s=marker_size, edgecolor='none')
        
        if show_colorbar:
            cbar = self._set_colorbar(axs, self.fig, CS=CS, cbar_mappable=cbar_mappable, 
                                      cbar_text=cbar_text)        
        self._set_labels(axs, x_label=x_label, y_label=y_label, title_label=title_label)
        self._set_tickers(axs, tick_multiplicator=tick_multiplicator)
        return self.fig, axs

    def _set_figure(self, fig=None, axs=None):
        if axs is None: 
            self.fig, axs = plt.subplots(constrained_layout=True)
        else:
            self.fig = fig 
        return self.fig, axs

    @classmethod
    def _set_colorbar(cls, axs, fig, CS=None, cbar_mappable=None, cbar_text:str=None):
        if cbar_mappable is None:
                cbar = fig.colorbar(CS, ax=axs)
        else:
            cbar = fig.colorbar(cbar_mappable, ax=axs)         
        cbar.ax.set_ylabel(cbar_text)
        return cbar

    @classmethod
    def _set_labels(cls, axs, x_label:str='', y_label:str='', title_label:str=None):
        axs.set_ylabel(y_label)
        axs.set_xlabel(x_label)
        if title_label is not None: axs.set_title(title_label)

    @classmethod
    def _set_tickers(cls, axs, tick_multiplicator=[None]*4):
        if all(tick_multiplicator[:2]):
            axs.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            axs.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            axs.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            axs.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))

    def Plotq3D(self, xi, yi, zi, fig=None, ax=None, x_label:str='', y_label:str='',
                title_label:str=None, z_label:str='', show_colorbar:bool=False,
                tick_multiplicator:list=[None, None, None, None],
                xmin=None, xmax=None, ymin=None, ymax=None,
                vmin=None, vmax=None, cbar_mappable=None, norm=None,
                color_map='viridis', show_contour_lines:bool=False, 
                cbar_text:str=None, marker='o', marker_size=None,
                plot_controur:bool=False, plot_scatter:bool=True, 
                interpolation_method='linear', interpolation_points:int = 20, 
                colorbar_scale:bool='log', CountFig=None, show_plot:bool=False,
                save_file_name:str='test', savefigure:bool=False, **kwargs_savefig):
        
        if plot_scatter:
            fig, ax = self._PlotScatter(xi, yi, zi, fig=fig, axs=ax, x_label=x_label, y_label=y_label,
                                         title_label=title_label, z_label=z_label, show_colorbar=show_colorbar,
                                         tick_multiplicator=tick_multiplicator, marker_size=marker_size,
                                         vmin=vmin, vmax=vmax, cbar_mappable=cbar_mappable, norm=norm,
                                         color_map=color_map, cbar_text=cbar_text, marker=marker)
        elif plot_controur:
            ##### Generate data with interpolation
            xii, yii, zii = self.InterPolation(xi,yi,zi, method=interpolation_method,
                                               interpolation_points=interpolation_points)
            fig, ax = self._PlotContour(xii, yii, zii, fig=fig, axs=ax, x_label=x_label, y_label=y_label,
                                         title_label=title_label, z_label=z_label, show_colorbar=show_colorbar,
                                         tick_multiplicator=tick_multiplicator, cbar_text=cbar_text,
                                         vmin=vmin, vmax=vmax, cbar_mappable=cbar_mappable, norm=norm,
                                         color_map=color_map, show_contour_lines=show_contour_lines)
        else:
            raise AttributeError('Not implemented')
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        self.save_figure(save_file_name, savefig=savefigure, show_plot=show_plot,
                         fig=fig, CountFig=CountFig, **kwargs_savefig)
        return fig, ax