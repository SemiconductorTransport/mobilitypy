from matplotlib import colors, cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri 
import matplotlib.ticker as ticker
from pathlib import Path
import os
import numpy as np

## ============================================================================
class general_plot_functions:
    def __init__(self):
        pass
    
    @classmethod
    def save_figs(cls, fig, filename:str='test', figs_path='.', savefigure:bool=False, 
                   FigFormat='.png', FigDpi:int=75):
        if savefigure and fig is not None:
            Path(figs_path).mkdir(parents=True, exist_ok=True)
            filename_ = f'{filename}{FigFormat}'
            fig.savefig(os.path.join(figs_path, filename_), bbox_inches='tight', dpi=FigDpi) 
            plt.close()
        return
        
class PlotQuasi3DFuns(general_plot_functions):
    def __init__(self):
        pass    
    
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
                    FigDpi:int=75, FigFormat='.png', figs_path='.', 
                    vmin=None, vmax=None, cbar_mappable=None, norm=None,
                    color_map='viridis', show_contour_lines:bool=False, 
                    cbar_text:str=None, show_colorbar:bool=False,
                    filename:str='test', savefigure:bool=False):
        # Set up the figure
        if axs is None: 
            self.fig, axs = plt.subplots(constrained_layout=True)
        else:
            self.fig = fig 

        # Plot linear interpolation to quad grid.
        CS = axs.contourf(xi, yi, zi, cmap=color_map)
        
        if norm is not None: CS.set_norm(norm)
        
        if show_contour_lines: 
            CS2 = axs.contour(CS, levels=CS.levels, colors='k')

        if show_colorbar:
            if cbar_mappable is None:
                cbar = fig.colorbar(CS, ax=axs)
            else:
                cbar = fig.colorbar(cbar_mappable, ax=axs)
            
            cbar.ax.set_ylabel(cbar_text)

        # Plot the triangulation.
        #cs = axs.tricontourf(triang, ZZ, vmin = vmin, vmax = vmax)
        #axs.tricontour(cs, colors='k', linewidths=1)
        axs.set_ylabel(y_label)
        axs.set_xlabel(x_label)
        if title_label is not None: axs.set_title(title_label)

        if all(tick_multiplicator[:2]):
            axs.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            axs.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            axs.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            axs.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))
            
        self.save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, axs

    def _PlotScatter(self, xi, yi, zi, fig=None, axs=None, x_label:str='', y_label:str='',
                    title_label:str=None, z_label:str='', show_colorbar:bool=False,
                    tick_multiplicator:list=[None, None, None, None],
                    FigDpi:int=75, FigFormat='.png', figs_path='.', 
                    vmin=None, vmax=None, cbar_mappable=None, norm=None,
                    color_map='viridis', cbar_text:str=None, marker='o',
                    filename:str='test', savefigure:bool=False):
        # Set up the figure
        if axs is None: 
            self.fig, axs = plt.subplots(constrained_layout=True)
        else:
            self.fig = fig 

        CS = axs.scatter(xi, yi, c=zi, cmap= color_map, vmin=vmin, vmax=vmax,marker=marker)

        if show_colorbar:
            if cbar_mappable is None:
                cbar = fig.colorbar(CS, ax=axs)
            else:
                cbar = fig.colorbar(cbar_mappable, ax=axs)
            
            cbar.ax.set_ylabel(cbar_text)

        axs.set_ylabel(y_label)
        axs.set_xlabel(x_label)
        if title_label is not None: axs.set_title(title_label)

        if all(tick_multiplicator[:2]):
            axs.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            axs.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            axs.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            axs.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))
            
        self.save_figs(fig, filename=filename, figs_path=figs_path, 
                        savefigure=savefigure, FigFormat=FigFormat, FigDpi=FigDpi)
        return fig, axs

    def Plot2D(self, xi, yi, zi, fig=None, ax=None, x_label:str='', y_label:str='',
               title_label:str=None, z_label:str='', show_colorbar:bool=False,
               tick_multiplicator:list=[None, None, None, None],
               FigDpi:int=75, FigFormat='.png', figs_path='.', 
               vmin=None, vmax=None, cbar_mappable=None, norm=None,
               color_map='viridis', show_contour_lines:bool=False, 
               cbar_text:str=None, marker='o', plot_controur:bool=False,
               plot_scatter:bool=True, interpolation_method='linear',
               interpolation_points:int = 20, colorbar_scale:bool='log',
               filename:str='test', savefigure:bool=False):
        
        if plot_scatter:
            return self._PlotScatter(xi, yi, zi, fig=fig, axs=ax, x_label=x_label, y_label=y_label,
                                     title_label=title_label, z_label=z_label, show_colorbar=show_colorbar,
                                     tick_multiplicator=tick_multiplicator,
                                     FigDpi=FigDpi, FigFormat=FigFormat, figs_path=figs_path, 
                                     vmin=vmin, vmax=vmax, cbar_mappable=cbar_mappable, norm=norm,
                                     color_map=color_map, cbar_text=cbar_text, marker=marker, 
                                     filename=filename, savefigure=savefigure)
        elif plot_controur:
            ##### Generate data with interpolation
            xii, yii, zii = self.InterPolation(xi,yi,zi, method=interpolation_method, 
                                                  interpolation_points=interpolation_points)
            return self._PlotContour(xii, yii, zii, fig=fig, axs=ax, x_label=x_label, y_label=y_label,
                                     title_label=title_label, z_label=z_label, show_colorbar=show_colorbar,
                                     tick_multiplicator=tick_multiplicator,
                                     FigDpi=FigDpi, FigFormat=FigFormat, figs_path=figs_path, 
                                     vmin=vmin, vmax=vmax, cbar_mappable=cbar_mappable, norm=norm,
                                     color_map=color_map, show_contour_lines=show_contour_lines,
                                     cbar_text=cbar_text, filename=filename, savefigure=savefigure)
        else:
            raise AttributeError('Not implemented')