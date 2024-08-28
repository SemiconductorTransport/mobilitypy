import numpy as np
import matplotlib.pyplot as plt
from ._general_plot_functions import _GeneratePlots
import matplotlib.ticker as ticker

### ===========================================================================
class _plot_mobilities(_GeneratePlots):
    """
    Plotting mobility figures.

    """
    def __init__(self, save_figure_dir='.'):
        """
        Initialize the plotting class.

        Parameters
        ----------
        save_figure_dir : str/path, optional
            Directory where to save the figure. The default is current directory.
       """
        _GeneratePlots.__init__(self, save_figure_dir=save_figure_dir)
    
    def _plot(self, results, fig=None, ax=None, save_file_name=None, CountFig=None, ymin=None, 
              ymax=None, xmax=None, xmin=None, y_scale_log:bool=True, mode:str= '2d_carrier_mobility',
              mobility_model:str='Bassaler', annotate_pos=(0,0), show_right_ticks:bool=False,
              title_text:str=None, xaxis_label:str='Composition', ls_2d='-', 
              yaxis_label:str=r'Electron mobility ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)',   
              color='gray', color_map='viridis', show_legend:bool=False, 
              show_colorbar:bool=False, colorbar_label:str=None, savefig:bool=False,
              vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        """
        

        Parameters
        ----------
        fig : matplotlib.pyplot figure instance, optional
            Figure instance to plot on. The default is None.
        ax : matplotlib.pyplot axis, optional
            Figure axis to plot on. If None, new figure will be created.
            The default is None.
        save_file_name : str, optional
            Name of the figure file. If None, figure will be not saved. 
            The default is None.
        CountFig: int, optional
            Figure count. The default is None.
        ymin : float, optional
            Minimum in y. The default is None.
        ymax : float, optional
            Maximum in y. The default is None.
        yaxis_label : str, optional
            Y-axis label text. The default is 'E (eV)'.
        color : str/color, optional
            Color of plot of unfolded band structure. The color of supercell
            band structures is gray. The default is 'gray'.
        color_map: str/ matplotlib colormap
            Colormap for density plot. The default is viridis.
        show_legend : bool, optional
            If show legend or not. The default is True.
        show_colorbar : bool, optional
            Plot the colorbar in the figure or not. If fig=None, this is ignored.
            The default is False.
        colorbar_label : str, optional
            Colorbar label. The default is None. If None, ignored.
        vmin, vmax : float, optional
            vmin and vmax define the data range that the colormap covers. 
            By default, the colormap covers the complete value range of the supplied data.
        show_plot : bool, optional
            To show the plot when not saved. The default is True.
        **kwargs_savefig : dict
            The matplotlib keywords for savefig function.
        
        Raises
        ------
        ValueError
            If plot mode is unknown.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            Figure instance. If ax is not None previously generated/passed fig instance
            will be returned. Return None, if no fig instance is inputed along with ax.
        ax : Axis instance
            Figure axis instance.
        CountFig: int or None
            Figure count.

        """
        if ax is None: 
            self.fig, ax = plt.subplots(constrained_layout=True)
        else:
            self.fig = fig 
            
        if yaxis_label is None: yaxis_label=''
        if xaxis_label is None: xaxis_label=''
 
        if mode == '2d_carrier_mobility':
            if mobility_model=='Bassaler':
                ax, return_plot = self._plot_2d_carrier_mobilities(results, ax, annotate_pos=annotate_pos, color=color)
        elif mode == 'plane_2d':
                ax, return_plot = self._plot_2d_plane(results, ax, color=color, ls=ls_2d)
        else:
            raise ValueError("Unknownplot mode: '{}'".format(mode))
            
        if show_colorbar and (self.fig is not None):
            cbar = self.fig.colorbar(return_plot, ax=ax)
            if colorbar_label is not None:
                cbar.set_label(colorbar_label)
        
        if y_scale_log: ax.set_yscale('log')
        ax.set_ylabel(yaxis_label)
        ax.set_xlabel(xaxis_label)
        ax.set_ylim([ymin, ymax])
        ax.set_xlim([xmin, xmax])

        if title_text is not None: ax.set_title(title_text)
        
        if show_right_ticks:
            ax.yaxis.set_ticks_position('both')

        if save_file_name is None:
            if show_plot: plt.show()
        else:
            CountFig = self._save_figure(save_file_name, savefig=savefig, show_plot=show_plot, 
                                         fig=self.fig, CountFig=CountFig, **kwargs_savefig)
        return self.fig, ax, CountFig

    @classmethod          
    def _plot_2d_carrier_mobilities(cls, mobility_df, ax, annotate_pos=(0,0), color=None):
        """

        Parameters
        ----------
        mobility_df : pandas dataframe
            Data to plot.
        ax : matplotlib.pyplot axis, optional
            Figure axis to plot on.
        color : matplotlib color/str, optional
            Color of the plots. The default is None.
        plot_colormap : bool, optional
        color_map : str/ matplotlib colormap
            Colormap for density plot.The default is 'viridis'.
        min_z : float, optional
            Minimum in the color scale. 
            The default is None. If None, determined from the data array supplied.
        max_z : float, optional
            Maximum in the color scale.
            The default is None. If None, determined from the data array supplied.

        Returns
        -------
        ax : matplotlib.pyplot axis
            Figure axis to plot on. 
        cmap_mappable :
            Figure instance or colormap instance for colorbar.

        """
        comp_ = np.array(mobility_df['comp'], dtype=float)
        for mu in list(mobility_df.keys())[1:]:
            ls='--' if 'TOT' in mu else '-'
            YY = np.array(mobility_df[mu], dtype=float)
            pp, = ax.plot(comp_, YY, ls=ls, color=color)
            color_pp = pp.get_color()
            ax.annotate(mu, (comp_[annotate_pos[0]], YY[annotate_pos[1]]), xycoords='data',
                        color=color_pp, xytext=(0, -20),  # -20 points vertical offset.
                        textcoords='offset points', ha='center', va='bottom', size=18)
        return ax, None

    @classmethod          
    def _plot_2d_plane(cls, plot_data, ax, color=None, ls='-'):
        pp, = ax.plot(plot_data[:, 0], plot_data[:, 1], ls=ls, color=color)
        return ax, None

    @classmethod
    def set_tickers(cls, axs, tick_multiplicator=[None]*4):
        if all(tick_multiplicator[:2]):
            axs.xaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[0]))
            axs.xaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[1]))
        if all(tick_multiplicator[2:]):
            axs.yaxis.set_major_locator(ticker.MultipleLocator(base=tick_multiplicator[2]))
            axs.yaxis.set_minor_locator(ticker.MultipleLocator(base=tick_multiplicator[3]))