from .src import _Mobility2DEG, _AlloyParams
from .utilities import _plot_mobilities

## ==============================================================================
class AlloyParams(_AlloyParams):
    def __init__(self):
        pass
            
    def get_alloy_params(self, system='ternary', compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN'):
        _AlloyParams.__init__(self, compositions=compositions, binaries=binaries, alloy=alloy)
        return self._get_alloy_params(system=system)

class Mobility2DEG(_Mobility2DEG):
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN', 
                 system='ternary', print_log=None):
        _Mobility2DEG.__init__(self, compositions=compositions, binaries=binaries, 
                               alloy=alloy, system=system, print_log=print_log)
        
    def calculate_mobility(self, n_2d=0.1, rms_roughness=0.1, corr_len=1, n_dis=1, f_dis=0.1, 
                           T=300, alloy_disordered_effect:bool=False,
                           interface_roughness_effect:bool=False,
                           dislocation_effect:bool=False,
                           deformation_potential_effect:bool=False, 
                           piezoelectric_effect:bool=False,
                           acoustic_phonon_effect:bool=False,
                           polar_optical_phonon_effect:bool=False,
                           total_mobility:bool=True,
                           mobility_model='Bassaler'):
        self.alloy_disordered_effect_=alloy_disordered_effect,
        self.interface_roughness_effect_=interface_roughness_effect,
        self.dislocation_effect_=dislocation_effect,
        self.deformation_potential_effect_=deformation_potential_effect,
        self.piezoelectric_effect_=piezoelectric_effect,
        self.acoustic_phonon_effect_=acoustic_phonon_effect,
        self.polar_optical_phonon_effect_=polar_optical_phonon_effect,
        self.total_mobility_=total_mobility,
        self.mobility_model_=mobility_model
        return self._calculate_mobility(n_2d=n_2d, rms_roughness=rms_roughness, corr_len=corr_len, 
                                        n_dis=n_dis, f_dis=f_dis, T=T)

    def calculate_figure_of_merit(self, n_2d, mobility, mode:str='LFOM', 
                                   direct_bandgap:bool=True, indirect_bandgap:bool=False):
        return self._calculate_figure_of_merit(n_2d, mobility, mode=mode, 
                                               direct_bandgap=direct_bandgap, 
                                               indirect_bandgap=indirect_bandgap)


class Plottings(_plot_mobilities):   
    def __init__(self, save_figure_dir='.'):
        """
        Intializing BandUPpy Plotting class.

        Parameters
        ----------
        save_figure_dir : str, optional
            Directory where to save the figure. The default is current directory.

        """
        self.save_figure_directory = save_figure_dir

    def plot_2d(self, data2plot, fig=None, ax=None, save_file_name=None, CountFig=None, ymin=None,
                ymax=None, xmax=None, xmin=None, y_scale_log:bool=True, show_right_ticks:bool=False,
                yaxis_label:str='', xaxis_label:str='', color=None, color_map='viridis', 
                show_legend:bool=False, show_colorbar:bool=False, colorbar_label:str=None, 
                savefig:bool=True, vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        
        _plot_mobilities.__init__(self, save_figure_dir=self.save_figure_directory)
        return self._plot(data2plot, fig=fig, ax=ax, save_file_name=save_file_name, 
                          CountFig=CountFig, ymin=ymin, ymax=ymax, xmax=xmax, xmin=xmin, 
                          y_scale_log=y_scale_log, mode='plane_2d', yaxis_label=yaxis_label, 
                          xaxis_label=xaxis_label, color=color, show_right_ticks=show_right_ticks,
                          color_map=color_map, show_legend=show_legend, show_colorbar=show_colorbar, 
                          colorbar_label=colorbar_label, savefig=savefig,
                          vmin=vmin, vmax=vmax, show_plot=show_plot, **kwargs_savefig)
    
    def plot_2deg_mobilities(self, mobility_dataframe, fig=None, ax=None, save_file_name=None, CountFig=None, ymin=None, 
                             ymax=None, xmax=None, xmin=None, y_scale_log:bool=True, mode:str= '2deg_mobility',
                             mobility_model:str='Bassaler',annotate_pos=(0,0), show_right_ticks:bool=False,
                             yaxis_label:str=r'Electron mobility ($\mathrm{cm}^2\mathrm{V}^{-1}\mathrm{s}^{-1}$)',
                             xaxis_label:str='Composition', color=None, color_map='viridis', show_legend:bool=False, 
                             show_colorbar:bool=False, colorbar_label:str=None, savefig:bool=True,
                             vmin=None, vmax=None, show_plot:bool=True, **kwargs_savefig):
        
        _plot_mobilities.__init__(self, save_figure_dir=self.save_figure_directory)
        return self._plot(mobility_dataframe, fig=fig, ax=ax, save_file_name=save_file_name, 
                          CountFig=CountFig, ymin=ymin, ymax=ymax, xmax=xmax, xmin=xmin, 
                          annotate_pos=annotate_pos, show_right_ticks=show_right_ticks,
                          y_scale_log=y_scale_log, mode= mode, yaxis_label=yaxis_label, 
                          xaxis_label=xaxis_label, color=color, mobility_model=mobility_model,
                          color_map=color_map, show_legend=show_legend, show_colorbar=show_colorbar, 
                          colorbar_label=colorbar_label, savefig=savefig,
                          vmin=vmin, vmax=vmax, show_plot=show_plot, **kwargs_savefig)