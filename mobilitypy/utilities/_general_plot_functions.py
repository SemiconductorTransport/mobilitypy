import matplotlib.pyplot as plt
from pathlib import Path
### ===========================================================================

class _GeneratePlots:
    def __init__(self, save_figure_dir='.'):
        """
        Initialize the plotting class.

        Parameters
        ----------
        save_figure_dir : str/path, optional
            Directory where to save the figure. The default is current directory.

        Returns
        -------
        None.

        """
        self.save_figure_folder = save_figure_dir
        params = {'figure.figsize': (8, 6),
                  'legend.fontsize': 18,
                  'legend.title_fontsize': 18,
                  'axes.labelsize': 24,
                  'axes.titlesize': 24,
                  'xtick.labelsize':24,
                  'xtick.major.width':2,
                  'xtick.major.size':5,
                  'xtick.minor.width':2,
                  'xtick.minor.size':3,
                  'ytick.labelsize': 24,
                  'ytick.major.width':2,
                  'ytick.major.size':5,
                  'ytick.minor.width':2,
                  'ytick.minor.size':3,
                  'errorbar.capsize':2}
        plt.rcParams.update(params)
        plt.rc('font', size=24)

    def _save_figure(self,fig_name, savefig:bool=True, show_plot:bool=True,
                     fig=None, CountFig=None, **kwargs_savefig):
        """
        The function saves the matplotlib figure.

        Parameters
        ----------
        fig_name : str
            name of the figure to use including figure extension.
        savefig : bool, optional
            Save the figure or not. The default is True.
        show_plot : bool, optional
            Show the plot or not. The default is True.
        fig : matplotlib figure instance, optional
            matplotlib figure instance. The default is None.
        CountFig : int, optional
            Figure count/index. The default is None.
        **kwargs_savefig : matplotlib savefig kwargs
            Any other matplotlib savefig keywords arguments.

        Returns
        -------
        CountFig : int
            Figure count/index.

        """
        if not savefig: 
            if show_plot: plt.show()
            return CountFig

        Path(self.save_figure_folder).mkdir(parents=True, exist_ok=True)
        if fig is not None:
            fig.savefig(f'{self.save_figure_folder}/{fig_name}', 
                        bbox_inches='tight', **kwargs_savefig)
        else:
            plt.savefig(f'{self.save_figure_folder}/{fig_name}', 
                        bbox_inches='tight', **kwargs_savefig)
        plt.close()
        if CountFig is not None: CountFig += 1
        return CountFig

    def save_figure(self,fig_name, savefig:bool=True, show_plot:bool=True,
                     fig=None, CountFig=None, **kwargs_savefig):
        """
        The function saves the matplotlib figure.

        Parameters
        ----------
        fig_name : str
            name of the figure to use including figure extension.
        savefig : bool, optional
            Save the figure or not. The default is True.
        show_plot : bool, optional
            Show the plot or not. The default is True.
        fig : matplotlib figure instance, optional
            matplotlib figure instance. The default is None.
        CountFig : int, optional
            Figure count/index. The default is None.
        **kwargs_savefig : matplotlib savefig kwargs
            Any other matplotlib savefig keywords arguments.

        Returns
        -------
        CountFig : int
            Figure count/index.

        """
        return self._save_figure(fig_name, savefig=savefig, show_plot=show_plot,
                                 fig=fig, CountFig=CountFig, **kwargs_savefig)
