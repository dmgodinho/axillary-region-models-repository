# Assigning Dielectric Properties to MRI-based models

This is a brief guide explaining how users can assign dielectric properties to the MRI-based axillary region models of this repository. We first detail how a Debye model is defined, and then, we present two alternatives to assign dielectric properties to the MRI-segmented tissues using Debye models. 

A Debye model can be used to simulate the frequency-dependent dielectric properties of biological tissues and can be written as:

![Eq1](/images/eq1.PNG)

where ![\omega](https://latex.codecogs.com/svg.latex?\omega) is the angular frequency, ![\Delta\epsilon=\epsilon_S-\epsilon_\infty](https://latex.codecogs.com/svg.latex?\Delta\epsilon=\epsilon_S-\epsilon_\infty), ![\epsilon_\infty](https://latex.codecogs.com/svg.latex?\epsilon_\infty) is the permittivity at ![\omega=\infty](https://latex.codecogs.com/svg.latex?\omega=\infty), ![\epsilon_S](https://latex.codecogs.com/svg.latex?\epsilon_S) and ![\sigma_S](https://latex.codecogs.com/svg.latex?\sigma_S) are the static permittivity and conductivity at ![\omega=0](https://latex.codecogs.com/svg.latex?\omega=0), and ![\tau](https://latex.codecogs.com/svg.latex?\tau) is the relaxation time constant.

The ![\epsilon^*](https://latex.codecogs.com/svg.latex?\epsilon^*) is a complex number which can be written as ![\epsilon'-j\epsilon''](https://latex.codecogs.com/svg.latex?\epsilon'-j\epsilon''). The real part ![\epsilon'](https://latex.codecogs.com/svg.latex?\epsilon') corresponds to the relative permittivity and the conductivity can be calculated from the imaginary part as ![\sigma=\omega\epsilon_0\epsilon''](https://latex.codecogs.com/svg.latex?\sigma=\omega\epsilon_0\epsilon'').

## 1. Recommended Debye Models of Dielectric Properties

The dielectric properties of the tissues included in the axillary region models can be inserted in the electromagnetic simulation software inputting the Debye models parameters [1-3], as shown in Table I.

![Table1](/images/table1.PNG  =150x)

## 2. Creation of Dielectric Property Maps

Some software packages for electromagnetic simulations may allow the input of dielectric property maps for each frequency. In order to apply this method, each model is represented by two files which can be used to create these maps: the _label\_map_ and _weight\_map_ files.

The _label\_map_ file has a matrix where each voxel is represented by one of 7 labels (background and 6 groups of tissues). Each label corresponds to a tissue type or group of tissues, and defines which dielectric properties curves should be used to interpolate the voxel intensities into the corresponding dielectric properties, as shown in Table II.

The _weight\_map_ file includes the weight assigned to each voxel which was calculated from the linear interpolation using the voxel intensities:

![Eq2](/images/eq2.PNG =20x)

where _v_ is the voxel intensity, _v\_min_ and _v\_max_ and are the minimum and maximum voxel intensity of the corresponding cluster, respectively. _w_ is a value between 0 and 1.

The corresponding dielectric property value for each voxel can be calculated using the following equation:

![Eq2](/images/eq3.PNG)

where _c\_upper_ and _c\_lower_ are the dielectric property values at the upper and lower curves considered for interpolation, respectively, at a given frequency _f_.

![Table2](/images/table2.PNG)

## 3. Example of a Python script for permittivity map calculation

```python
import numpy as np

def calculate_properties_from_debye(debye_parameters, curve, frequency):
    """
    Function to calculate dielectric properties curve based on parameters of Debye model.
    Inputs:
        - "debye_parameters": parameters of Debye model for each
        curve
        - "curve": string of the curve name
        - "frequency": frequency value(s), in Hz
    Outputs:
        - "e_r": relative permittivity values
        - "sigma_eff": effective conductivity values
    """
    e0 = 8.85e-12
    w = 2*np.pi*frequency
    e = debye_parameters[curve]['e_inf'] + \
            debye_parameters[curve]['delta_e']/ \
            (1+(1j*w*debye_parameters[curve]['tau'])**(1-debye_parameters[curve]['alpha'])) + \
            debye_parameters[curve]['sigma']/(1j*w*e0)

    e_r = np.real(e)
    sigma_eff = w*e0*np.absolute(np.imag(e))
    return e_r, sigma_eff 

```
```python
# Load the files

# model = LOAD MAT OR RAW label_map FILE
# w = LOAT MAT OR RAW weight_map FILE


# Define parameters of Debye models
debye_parameters = {}
debye_parameters['Skin'] = {'e_inf':15.93, 'delta_e':23.83, 'tau':13e-12, 'alpha':0, 'sigma':0.831}
# [...] apply same rationale for remaining curves


# Calculation of relative permittivity and conductivity of Debye model curve for the specified frequency
# e_r is the relative permittivity
# sigma_eff is the effective conductivity
# f is the frequency of interest in Hz 

f = 5
e_r = {}
sigma_eff = {}
for tissue_name in debye_parameters.keys():
    e_r[tissue_name], sigma_eff[tissue_name] = calculate_properties_from_debye(debye_parameters, tissue_name, f)


# Create matrix of zeros with the same dimensions of the model
permittivity_map = np.zeros(model.shape)

# For loop scanning each voxel 
for x in range(model.shape[0]):
    for y in range(model.shape[1]):
        for z in range(model.shape[2]):
            if model[x,y,z] == 0:
                c_upper = 0
                c_lower = 0
            elif model[x,y,z] == -2:                    
                c_upper = e_r['Skin']
                c_lower = e_r['Skin']
            elif model[x,y,z] == -1:                    
                c_upper = e_r['Lung']
                c_lower = e_r['Lung']
            elif model[x,y,z] == 1:                    
                c_upper = e_r['Adipose Min']    
                c_lower = e_r['Adipose Q1']
            elif model[x,y,z] == 2:                    
                c_upper = e_r['Adipose Q1']
                c_lower = e_r['Adipose Q3']
            elif model[x,y,z] == 3:                    
                c_upper = e_r['Fibroglandular Q1']
                c_lower = e_r['Fibroglandular Q3']
            elif model[x,y,z] == 4:                    
                c_upper = e_r['Fibroglandular Q3']
                c_lower = e_r['Fibroglandular Max']
            
            permittivity_map[x,y,z] = w[x,y,z]*c_upper + (1-w[x,y,z])*c_lower

```

## 4. References

[1] D. M. Godinho, J. M. Felício, T. Castela, N. A. Silva, M. L. Orvalho, C. A. Fernandes, R. C. Conceição, "MRI-based Axillary Numerical Models and Extraction of Axillary Lymph Nodes Dielectric Properties for Microwave Imaging" (under review in Medical Physics, 2021).

[2] S. Gabriel, R. W. Lau, and C. Gabriel, "The dielectric properties of biological tissues: III Parametric models for the dielectric spectrum of tissues," Phys. Med. Biol., vol. 41, no. 11, pp. 2271–2293, 1996.

[3] M. Lazebnik et al., "A large-scale study of the ultrawideband microwave dielectric properties of normal breast tissue obtained from reduction surgeries," Phys. Med. Biol., vol. 52, no. 10, pp. 2637–3656, 2007.


## Cite Us

To cite this work, please cite as: 	

D. M. Godinho, J. M. Felício, T. Castela, N. A. Silva, M. L. Orvalho, C. A. Fernandes, R. C. Conceição, "MRI-based Axillary Numerical Models and Extraction of Axillary Lymph Nodes Dielectric Properties for Microwave Imaging" (under review in Medical Physics, 2021).