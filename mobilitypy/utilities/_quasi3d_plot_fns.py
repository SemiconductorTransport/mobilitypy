from matplotlib import colors, cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri 
import matplotlib.ticker as ticker
import numpy as np
from ._general_plot_functions import _GeneratePlots

## ============================================================================
class PlotQuasi3DFuns(_GeneratePlots):
    """
    The function in this class are used to plot 2d heatmaps.
    """
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
        """
        This function perform interpolation using matplotlib tri interpolation
        library.

        Parameters
        ----------
        x_values : 1d numpy array
            x coordinates. Will be used to generate mesh grid and 
            passed to matplotlib.tri.Triangulation() function.
        y_values : 1d numpy array
            y coordinates. Will be used to generate mesh grid and 
            passed to matplotlib.tri.Triangulation() function.
        z_values : 1d numpy array
            z coordinates. Will be used to generate mesh grid and 
            passed to matplotlib.tri.Triangulation() function.
        method : str, optional => ['linear']
            Method of interpolation. The default is 'linear'.
        interpolation_points : int, optional
            Numper of interpolation points. The default is 20.

        Returns
        -------
        xi : nd numpy array
            x coordinate after interpolation.
        yi : nd numpy array
            y coordinate after interpolation.
        zi : nd numpy array
            z coordinate after interpolation.

        """
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
        """
        This function plots 2d heatmaps.

        Parameters
        ----------
        xi : 1d numpy array
            x coordinates.
        yi : 1d numpy array
            y coordinates.
        zi : 1d numpy array
            z coordinates.
        fig : matplotlib.pyplot figure instance, optional
            Figure instance to plot on. The default is None.
        ax : matplotlib.pyplot axis, optional
            Figure axis to plot on. If None, new figure will be created.
            The default is None.
        x_label : str, optional
            x-axis label text. The default is ''.
        y_label : str, optional
            y-axis label. The default is ''.
        title_label : str, optional
            Figure title. The default is None.
        z_label : str, optional
            Label for colorbar. The default is ''.
        show_colorbar : bool, optional
            Show colorbar. The default is False.
        tick_multiplicator : list, optional
            Controls ticks in axes. 
            Format: [x_major_tick_multiplicator, x_minor_tick_multiplicator,
                     y_major_tick_multiplicator, y_minor_tick_multiplicator]
            The default is [None, None, None, None].
        ymin : float, optional
            Minimum in y. The default is None.
        ymax : float, optional
            Maximum in y. The default is None.
        xmin : float, optional
            Minimum in x. The default is None.
        xmax : float, optional
            Maximum in x. The default is None.
        vmin : float, optional
            Minimum value in colorbar. The default is None.
        vmax : float, optional
            Maximum value in colorbar. The default is None.
            When using scalar data and no explicit norm, vmin and vmax define 
            the data range that the colormap covers. By default, the colormap 
            covers the complete value range of the supplied data. It is an error 
            to use vmin/vmax when a norm instance is given (but using a str norm 
            name together with vmin/vmax is acceptable).
            This parameter is ignored if c is RGB(A).
        cbar_mappable : Matplotlib colormap mappable instance, optional
            Matplotlib colormap mappable. The default is None.
        norm : normstr or Normalize, optional
            The normalization method used to scale scalar data to the [0, 1] range 
            before mapping to colors using cmap. By default, a linear scaling is used, 
            mapping the lowest value to 0 and the highest to 1.
            If given, this can be one of the following:
                An instance of Normalize or one of its subclasses (see Colormap normalization).
                A scale name, i.e. one of "linear", "log", "symlog", "logit", etc. For a list of 
                available scales, call matplotlib.scale.get_scale_names(). In that case,
                a suitable Normalize subclass is dynamically generated and instantiated.
            This parameter is ignored if c is RGB(A). The default is None.
        color_map : cmapstr or Colormap, optional
            cmapstr or Colormap, default: rcParams["image.cmap"]. The default is 'viridis'.
        show_contour_lines : bool, optional
            The the contour line in contour plot. The default is False.
        cbar_text : str, optional
            Colorbar label. The default is None.
        marker : matplotlib marker, optional
            Marker for matplotlib plot. The default is 'o'.
        marker_size : matplotlib marker size, optional
            Marker size of matplotlib plot. The default is None.
        plot_controur : bool, optional
            Plot contour plot. The default is False.
        plot_scatter : bool, optional
            Plot scatter plot. The default is True.
        interpolation_method : str, optional
            Interpolation method. The default is 'linear'.
        interpolation_points : int, optional
            Numper of interpolation points. The default is 20.
        colorbar_scale : bool, optional
            A scale name, i.e. one of "linear", "log", "symlog", "logit", etc.
            The default is 'log'.
        CountFig : int, optional
            Figure count/index. The default is None.
        show_plot : bool, optional
            Show plot. The default is False.
        save_file_name : str, optional
            File name to save the figure. The default is 'test'.
        savefigure : bool, optional
            Save the figure. The default is False.
        **kwargs_savefig : TYPE
            DESCRIPTION.

        Raises
        ------
        AttributeError
            DESCRIPTION.

        Returns
        -------
        fig : matplotlib figure instance
        ax : matplotlib axis instance

        """
        
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