import numpy as np
import matplotlib.pyplot as plt
from ._general_plot_functions import _GeneratePlots

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
              ymax=None, xmax=None, xmin=None, y_scale_log:bool=True, mode:str= '2deg_mobility',
              mobility_model:str='Bassaler',
              yaxis_label:str=r'Electron mobility ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)',   
              xaxis_label:str='Composition', color='gray', color_map='viridis', show_legend:bool=False, 
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
 
        if mode == '2deg_mobility':
            if mobility_model=='Bassaler':
                ax, return_plot = self._plot_2deg_mobilities(results, ax, color=color)
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

        

        if save_file_name is None:
            if show_plot: plt.show()
        else:
            CountFig = self._save_figure(save_file_name, savefig=savefig, show_plot=show_plot, 
                                         fig=self.fig, CountFig=CountFig, **kwargs_savefig)
            plt.close()
        return self.fig, ax, CountFig

    @classmethod          
    def _plot_2deg_mobilities(cls, mobility_df, ax, color=None):
        """
        Plot the band centers and band width.

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
        comp_ = np.array(mobility_df.index, dtype=float)
        for mu in mobility_df:
            ls='--' if 'TOT' in mu else '-'
            pp, = ax.plot(comp_, mobility_df[mu], ls=ls, color=color)
            color_pp = pp.get_color()
            ax.annotate(mu[2:], (comp_[11], mobility_df[mu].iloc[11]), color=color_pp, xytext=(0, -20),  # -13 points vertical offset.
                                textcoords='offset points', ha='center', va='bottom', size=18)
        return ax, None