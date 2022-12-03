# indentationTools.jl

IndentationTools.jl is a package of help functions to process indentation experiments. It is particularly adapt at handling AFM-based indentation data. A switch is included to allow the import of force-indentation curves directly. Currently, this option is mostly used to compare against commercial software.

This repository is written with a specific work-flow in mind. Before diving into your analysis, take a day to read the code and check that used logic applies to your problem. 

## Essentials
The effect of a non-rigid indenter can be accounted for by defining a measured modulus $E_r$ via

$$\frac{1}{E_r} = \frac{1-\nu^2}{E} + \frac{1-\nu_i^2}{E_i}$$

where $E$ and $\nu$ are the properties of the tested material and $E_i$, $\nu_i$ are the properties of the indenter. The key underpinning of all the equations in the package is that the measured modulus $E_r$ is related to the unloading stiffness $S$ and the projected area of elastic contact $A$ by the equation

$$ S = \frac{dP}{dh} = \frac{2}{\sqrt{\pi}} E_r \sqrt{A}$$

## Unloading curve fitting equations

### Oliver-Pharr [1]
The classical method from reference [1] below. Assumes that the unloading curve can be described by the equation 

$$P(h) = B(h-h_f)^m $$

where $B$, $h_f$ and $m$ are fitting parameters. The initial unloading slope $S$ is found by differentiating this expression and evaluating it at the onset of unloading, 

$$ S = mB(h_{\textrm{max}}-h_f)^{m-1}$$

Note that the value of $m$ predicted by theory is 1 / 1.5 / 2.0 for a flat/paraboloid/conical indenter, respectively (p. 1575, [1]).

## Where can I read more about the equations implemented?
The basics of indentation testing is the theory of Hertzian contact. Beyond a general understanding of Hertz contact theory, I recommend the following articles to understand this repository:

    [1] Oliver, W. C., & Pharr, G. M. (1992). 
    An improved technique for determining hardness and elastic modulus using load and displacement sensing indentation experiments. 
    Journal of Materials Research, 7(6), 1564–1583. [https://doi.org/10.1557/JMR.1992.1564](https://doi.org/10.1557/JMR.1992.1564)
    
    [2] Oliver, W. C., & Pharr, G. M. (2004). 
    Measurement of hardness and elastic modulus by instrumented indentation: Advances in understanding and refinements to methodology. 
    Journal of Materials Research, 19(1), 3–20. [https://doi.org/10.1557/jmr.2004.19.1.3](https://doi.org/10.1557/jmr.2004.19.1.3) 
    
    [3] Feng, G., & Ngan, A. H. W. (2002). 
    Effects of creep and thermal drift on modulus measurement using depth-sensing indentation. 
    Journal of Materials Research, 17(3), 660–668. [https://doi.org/10.1557/JMR.2002.0094](https://doi.org/10.1557/JMR.2002.0094) 
    
    [4] Tang, B., & Ngan, A. H. W. (2003). 
    Accurate measurement of tip - Sample contact size during nanoindentation of viscoelastic materials.
    Journal of Materials Research, 18(5), 1141–1148. [https://doi.org/10.1557/JMR.2003.0156](https://doi.org/10.1557/JMR.2003.0156) 
    
    [5] Cheng, Y. T., & Cheng, C. M. (2005). 
    Relationships between initial unloading slope, contact depth, and mechanical properties for conical indentation in linear viscoelastic solids. 
    Journal of Materials Research, 20(4), 1046–1053. [https://doi.org/10.1557/JMR.2005.0141](https://doi.org/10.1557/JMR.2005.0141) 
    
    [6] Cheng, Y. T., & Cheng, C. M. (2005). 
    Relationships between initial unloading slope, contact depth, and mechanical properties for spherical indentation in linear viscoelastic solids. 
    Materials Science and Engineering A, 409(1–2), 93–99. [https://doi.org/10.1016/j.msea.2005.05.118](https://doi.org/10.1016/j.msea.2005.05.118) 
