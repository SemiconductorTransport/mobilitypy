# `mobilitypy`: Python package for mobility calculations in semiconductor heterostructures

<!-- =========================================================== -->

<!-- =========================================================== -->
<!-- Here absolute path is added since pypi documentation could not access package relative path-->
![](https://raw.githubusercontent.com/SemiconductorTransport/mobilitypy/refs/heads/main/imgs/mobilities_300K.png) | ![](https://raw.githubusercontent.com/SemiconductorTransport/mobilitypy/refs/heads/main/imgs/LFOM_norm_300K.png) | ![](https://raw.githubusercontent.com/SemiconductorTransport/mobilitypy/refs/heads/main/imgs/LFOM_norm_T_dependence.png)
:------------------------------:| :------------------------------:| :------------------------------:
AlN/AlGaN HEMT 2DEG mobilities | AlN/AlGaN HEMT lateral figure of merit | AlN/AlGaN HEMT lateral figure of merit (Temperature effect) 
<!-- =========================================================== -->

<!-- =========================================================== -->
## Systems
<!-- =========================================================== -->
### 1. High electron mobility transistors (HEMTs)
Details of the mobility models are described here: [2DEG mobility models](https://raw.githubusercontent.com/SemiconductorTransport/mobilitypy/refs/heads/main/docs/Mobility2DEGAnalyticalModels.pdf)

### 2. Polarization field effect transistors (PolFETs)
Details of the mobility models are described here: [3DEG mobility models](https://raw.githubusercontent.com/SemiconductorTransport/mobilitypy/refs/heads/main/docs/Mobility3DEGAnalyticalModels.pdf)

<!-- =========================================================== -->
## Developers and contributors
<!-- =========================================================== -->

__Developer of mobilitypy :__
[Badal Mondal](https://github.com/bmondal94) 

__Maintainer of mobilitypy :__
[Badal Mondal](https://github.com/bmondal94)

__mobilitypy Contributors:__  [Contributors](https://github.com/SemiconductorTransport/mobilitypy/graphs/contributors)

* We sincerely thank each and every contributors for their valuable input and support.

__Contact us:__ [Email developer/maintainer team](mailto:badalmondal.chembgc@gmail.com)

* If you would like to contribute to the development of `mobilitypy` or request new functionality, please get in touch with [us](mailto:badalmondal.chembgc@gmail.com) or open a pull request or discuss [here](https://github.com/SemiconductorTransport/mobilitypy/discussions). We appreciate and respect our users' views and are committed to providing the best experience possible. Your feedback is highly valued. We will be happy to support your request ASAP.

<!-- =========================================================== -->

<!-- =========================================================== -->
## Installation

### 1. Requirements
```
    1. python>=3.12
    2. numpy
    3. scipy
    4. matplotlib
    5. pathlib
    6. pandas
```

### 2. Installation using `pip` [*recommended]

```
    pip install mobilitypy
```

### 3. Installation from github repository

```
    git clone https://github.com/SemiconductorTransport/mobilitypy.git
    cd mobilitypy
    pip install .  
```
Or, 
```
    pip install git+https://github.com/SemiconductorTransport/mobilitypy.git@specific_branch
```

<!-- =========================================================== -->


<!-- =========================================================== -->
## Usage
__Wiki page__: [Welcome to mobilitypy](https://github.com/SemiconductorTransport/mobilitypy/wiki)

__Documentation__: [Package documentation](https://github.com/SemiconductorTransport/mobilitypy/wiki/01.-Package-documentation)

__Discussions__: [Discuss more about the package here](https://github.com/SemiconductorTransport/mobilitypy/discussions)

__Material database__: [here](https://github.com/SemiconductorTransport/mobilitypy/blob/main/mobilitypy/src/database.py)

__Tutorials__: [tutorial](https://github.com/SemiconductorTransport/mobilitypy/tree/main/tutorials)

<!-- =========================================================== -->
## Tips and tricks:

__FAQs__: [here](https://github.com/SemiconductorTransport/mobilitypy/wiki/02.-Frequently-asked-questions-(FAQs))

You can find a list of common user issues encountered while using this software [here](https://github.com/SemiconductorTransport/mobilitypy/wiki/02.-Frequently-asked-questions-(FAQs)). We appreciate and respect our users' views and are committed to providing the best experience possible. Your feedback is highly valued.

<!-- =========================================================== -->

<!-- =========================================================== -->
## Citations and references:

If you use `mobilitypy` in your work, please:

  * **State EXPLICITLY that you have used the mobilitypy code** (or a modified version of it, if this is the case), for instance, adding a sentence like:
  
    "The mobility calculation has been performed using the mobilitypy code."

  * **How to cite the package:** (use appropriate version number and doi corresponding to your installed mobilitypy)
>> Badal Mondal, "SemiconductorTransport/mobilitypy: version-V.V.V (vV.V.V))". Zenodo, 2026. doi: XXXXX

  * **Read and cite the following papers where applicable** (and the appropriate references therein):
    *  The analytical mobility models are implemented based on the following publications and references therein:
>> 1. J. Bassaler, J. Mehta, I. Abid, L. Konczewicz, S. Juillaguet, S. Contreras, S. Rennesson, S. Tamariz, M. Nemoz, F. Semond, J. Pernot, F. Medjdoub, Y. Cordier, P. Ferrandis, Al-Rich AlGaN Channel High Electron Mobility Transistors on Silicon: A Relevant Approach for High Temperature Stability of Electron Mobility. [Adv. Electron. Mater. 11, 2400069 (2025).](https://doi.org/10.1002/aelm.202400069)
>> 2. J. Zhang, Y. Hao, J. Zhang, J. Ni, The mobility of two-dimensional electron gas in AlGaN/GaN heterostructures with varied Al content. [Sci. China Ser. F-Inf. Sci. 51, 780–789 (2008).](https://doi.org/10.1007/s11432-008-0056-7)
>> 3. B. Mondal, P. Pampili, J. Mukherjee, D. Moran, P.J. Parbrook, S. Schulz, Interplay of carrier density and mobility in Al-rich (Al,Ga)N-channel HEMTs: Impact on high-power device performance potential. [APL Electronic Devices 1, 026117 (2025).](https://doi.org/10.1063/5.0277051) [preprint arXiv:2502.13809](https://doi.org/10.48550/arXiv.2502.13809)

  * **Other publications that used mobilitypy**:
>> TBA

__Bibliography file:__ Here is the [bibliography file](https://github.com/SemiconductorTransport/mobilitypy/wiki/03.-References-(bibliography-style)) for your convenience.

<!-- =========================================================== -->

<!-- =========================================================== -->
## Version release

Chekout out [version release history here](https://github.com/SemiconductorTransport/mobilitypy/wiki/04.-Release-history) for the full list of updates and upgrades.

<!-- =========================================================== -->

<!-- =========================================================== -->
## License
* [GNU General Public License v3.0](https://raw.githubusercontent.com/SemiconductorTransport/mobilitypy/refs/heads/main/LICENSE)
<!-- =========================================================== -->

<!-- =========================================================== -->
## Upcoming (TBD)
1. TBA
<!-- =========================================================== -->

