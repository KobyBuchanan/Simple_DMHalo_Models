# Simple_DMHalo_Models

## What is This?
In short, this program generates a galaxy from random realizations of model parameters. Currently Implemented models are the Isothermal Sphere and the NFW sphere. 

## Why did I make this?
This project originated as an assignment for a computational physics course I took during my sophomore year at UMass. While the project document that outlines the assignment is the property of UMass and cannot be shared here, the goal was to create toy galaxy models, plot their phase space, and use numerical integration techniques to simulate their evolution over time. The course introduced me to Python and its applications in astronomy.

I was never fully satisfied with how my original code turned out. Now, three years later, with more experience under, I have decided to revisit the project and completely overhaul it. The original project guidelines included numerous simplifications and assumptions to make these "toy models" manageable. While many of these simplifications remain in the current version, I am gradually working toward more complex models that incorporate scientifically justified assumptions 

## Known Issues
For the NFW Halo I was provided a table of 1-D velocity dispersions, I no longer have that file or a means to access it again. Currently I imploy the Isotropic Jeans equation to calculate the dispersions at my accepted radius distribution. This technique should work, however the galxies generated from the NFW model don't reach equalibrium. 

## To Do
* Make better plots and figures
* Finish adding NFW sphere
* Add remaining models from old project
* Add galaxy evolution
* Make more robust and accurate models 

## Acknowledgements 
I would like to acknowledge UMass Amherst and Dr. Neal Katz for providing the original project guidelines.  
The book ["Dynamics and Astrophysics of Galaxies"](https://galaxiesbook.org/index.html) by Dr. Jo Bovy has been immensely helpful in deepening my understanding of the equations used in this project. 

