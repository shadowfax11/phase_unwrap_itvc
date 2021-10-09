===ABOUT===

Contains code for performing ADMM-based optimization for unwrapping 2D phase maps. The methods by which one can unwrap phase maps are -
	- ITV: 
		Using first-order Itoh information alongwith TV-based and Irrotationality/Integrability-based regularizers
	- IHTV: 
		Using first-order Itoh information alongwith HTV-based and Irrotationality/Integrability-based regularizers
	- ITVC: 
		Using first and second order Itoh information alongwith TV-based and Irrotationality/Integrability-based regularizers

===CODE FROM OTHER SOURCES===
Terrain generation code from: 
Tucker McClure (2021). Automatic Terrain Generation (https://www.mathworks.com/matlabcentral/fileexchange/39559-automatic-terrain-generation), MATLAB Central File Exchange. Retrieved September 28, 2021. 

FSIM implementation code from:
http://www4.comp.polyu.edu.hk/~cslzhang/IQA/FSIM/FSIM.htm

MLV implementation code from:
Khosro Bahrami (2021). Maximum Local Variation (MLV) Code for sharpness assessment of images (https://www.mathworks.com/matlabcentral/fileexchange/49991-maximum-local-variation-mlv-code-for-sharpness-assessment-of-images), MATLAB Central File Exchange. Retrieved September 28, 2021. 

===USAGE===

Given noisy phase map 'Pn' only: 
To perform unwrapping, call the following MATLAB functions:
	For ITV:  wrapper_unwrap_itv(Pn,Pn,'ITV',save_str)
	For IHTV: wrapper_unwrap_itv3(Pn,Pn,'IHTV',save_str) 
	For ITVC: wrapper_unwrap_itvc(Pn,Pn,'ITVC',save_str)
where save_str is the string used for the .mat file used to store the results of the unwrapping

Given noisy phase map 'Pn' and ground truth 'P': 
To perform unwrapping, call the following MATLAB functions:
	For ITV:  wrapper_unwrap_itv(P,Pn,'ITV',save_str)
	For IHTV: wrapper_unwrap_itv3(P,Pn,'IHTV',save_str) 
	For ITVC: wrapper_unwrap_itvc(P,Pn,'ITVC',save_str)
where save_str is the string used for the .mat file used to store the results of the unwrapping

 