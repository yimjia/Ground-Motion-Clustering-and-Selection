# Ground Motion Clustering and Selection 
We developed a clustering-based ground motion (GM) selection method that can properly identify representative GM records to match the conditional spectra that a structure will probabilistically experience. The developed GM selection method includes four main steps: 1) leveraging domain-specific knowledge to pre-select candidate GMs; 2) using a convolutional autoencoder to learn low-dimensional underlying characteristics of candidate GMs’ response spectra – i.e., latent features; 3) performing k-means clustering to classify the learned latent features, equivalent to cluster the response spectra of candidate GMs; and 4) embedding the clusters in the conditional spectra-based GM selection. The selected GMs can match the conditional spectra’s mean and variability well while fully describing the candidate GMs.

For more information, please refer to the following:

Jia, Y., and Sasani, M. (2024). "Convolutional Autoencoder-Based Ground Motion Clustering and Selection", In Review (will be available after this paper is published). 
<br/><br/>

As prerequisites, users need to have knowledge of earthquake engineering to perform step one GM pre-selection and have access to run Python and Matlab codes. No training in machine learning is required to perform the developed clustering-based GM selection. 

To download the codes, pleaase navigate to the main page of this repository. Above the list of files, click the green **Code** button, and in the menu that appears, click **Download ZIP**. 

The "Model" folder includes the codes and a step-by-step instruction to use the developed clustering-based GM selection. 

The "Example" folder includes the code, data, and a instruction to reproduce the GM selection presented in Section 4.1 of the above-mentioned paper.

Note that the codes are run on Python 3.10.14 (Tensorflow 2.13.1 required), and Matlab R2023b. The authors recommend running the Python codes on Jupyter Notebook. Running the codes on different versions of Python, Tensorflow, or Matlab may result in compatibility issues. However, most potential compatibility issues can be easily resolved by following the suggestions in the error messages.
