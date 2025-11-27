# Package Documentation

## 1. Import modules
```
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import colors, cm
import matplotlib.ticker as ticker
from mobilitypy import AlloyParams, Mobility2DCarrier
from mobilitypy import Plottings, PlotQuasi3DFuns
```

## 2. Calculate mobilities
```
## Intialize the mobility class
mu2deg = Mobility2DCarrier(compositions=0.1, binaries=['AlN', 'GaN'], alloy='AlGaN', system='ternary', 
                           psedomorphic_strain=False, substrate=None, alloy_type='WZ',
                           print_log=None, eps_n_2d=1e-10)
                           
## Set parameters for mobility calculations
alloy_disordered_effect = True
interface_roughness_effect = True
dislocation_effect = True
deformation_potential_effect = True 
piezoelectric_effect = True
acoustic_phonon_effect = True
polar_optical_phonon_effect = True
total_mobility = True
mobility_model = 'Bassaler'
density_2deg = 0.1 # nm^-2
irf_rms_roughness = 0.3 # nm
irf_corr_length = 3.0 # nm
dislocation_density = 1e-4 # nm^-2
occup_dislocation = 0.3
t = 300 # K
T_corect_bandgap_in_LFOM = False                          
                           
# Calculate mobilities
mobility_ref = mu2deg.calculate_sheet_mobility(n_2d=density_2deg_comp_ref, rms_roughness=irf_rms_roughness,
                                               corr_len=irf_corr_length, n_dis=dislocation_density,
                                               f_dis=occup_dislocation, T=T,
                                               alloy_disordered_effect=alloy_disordered_effect,
                                               interface_roughness_effect=interface_roughness_effect,
                                               dislocation_effect=dislocation_effect,
                                               deformation_potential_effect=deformation_potential_effect,
                                               piezoelectric_effect=piezoelectric_effect,
                                               acoustic_phonon_effect=acoustic_phonon_effect,
                                               polar_optical_phonon_effect=polar_optical_phonon_effect,
                                               total_mobility=total_mobility,
                                               calculate_total_mobility_only=False,
                                               mobility_model=mobility_model,
                                               return_sc_rates=False)

```

## 3. Plot mobilities
### 3.1 Mobility plots 2d
```
plt2deg = Plottings(save_figure_dir='.')
fig, ax,_ = plt2deg.plot_2d_carrier_mobilities(mobility_ref, save_file_name='save_file_name.png',
                                               ymin=5e1, ymax=2e5, xmax=0.9, xmin=0.5, y_scale_log=True, 
                                               annotate_pos=(2,2), annotatetextoffset=(0,-20),show_right_ticks=True,
                                               mode='2d_carrier_mobility', yaxis_label=y_label, xaxis_label=x_label,
                                               color=None, color_map='viridis', savefig=True, show_plot=False)
#plt2deg.save_figure(save_file_name_, fig=fig, savefig=true, dpi=300, show_plot=False)
```
### 3.2 Heatmap plots
```
lpltq3d = PlotQuasi3DFuns(save_figure_dir='.')
##### Plot composition map 
cbar_mapable = cm.ScalarMappable(norm=norm, cmap='viridis')
fig, axt = lpltq3d.Plotq3D(xx, yy, zz, 
                          xmin=0.475, xmax=0.975, ymin=0.525, ymax=1.025,
                          x_label='x-label', y_label='y-label',
                          interpolation_method='linear',
                          interpolation_points = 20,
                          tick_multiplicator=[1,0.5,2,0.1],
                          title_label=None, 
                          cbar_mappable=cbar_mapable, norm=norm,
                          show_contour_lines=False, marker='s', marker_size=24**2,
                          cbar_text=z_label_,show_colorbar=False,
                          plot_controur=0, plot_scatter=True,
                          savefigure=False, show_plot=False)
cbar = fig.colorbar(cbar_mapable, ax=axt)         
cbar.ax.set_ylabel('z-label', size=18)
lpltq3d.save_figure('test.png', savefig=True, show_plot=True,
                    fig=fig, CountFig=None, dpi=300)
```
<!-- =========================================================== -->



##
__If you have new suggestions, please feel free to reach out to us. We are committed to providing the best experience for our users and greatly value your feedback.__

__Have fun with mobilitypy!__

__Best wishes,__  
__The mobilitypy Team__
