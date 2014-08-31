## Error Analysis Notes

# Abstract
These are pedagogical error analysis notes written for 1st/2nd year physics
undergraduate students as a response to feedback from students that there was
not sufficient coverage of these ideas in their coursework. The notes contain
Basic materials, roughly commensurate with what a typical undergraduate physics
major requires (basic concepts, error estimation in measurement and from
graphs, and linear error propagation), and Advanced Topics for keen students
wishing to learn about probability theory, the two major schools of statistics
(Bayesian and frequentist) and Monte Carlo techniques for non-linear error
propagation. A list of materials for further reading and a summary 'cheat sheet'
are also provided at the end. This open source document is licensed under a [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).
The 'official' version is maintained in a public Git repository at [Github](https://github.com/mikejbrown/error-analysis).
Please feel free to contribute by submitting a pull request!

# Building

This document was built using LyX 2.1.1 and tested on Fedora 20 and Ubuntu 14.04.
It is built using the [dexy](http://dexy.it/) template and build system for reproducible
document generation. The latest tested version is dexy 1.0.13.

You must also have python installed. The build process has been tested with python versions 2.7.x,
though this should not be a strict version requirement. :eyes:

To run `fit_lines.py` requires the [triangle plotting package](https://pypi.python.org/pypi/triangle_plot/0.0.6). Note that this is `pip install triangle_plot` _not_ `pip install triangle` :exclamation:

To set up:

1. Install dexy if need be: `pip install dexy` (see [here](http://dexy.it/install/) for more details on different systems, building from source etc.).
2. Get the triangle plotting package if need be: `pip install triangle_plot` (see [here](https://github.com/dfm/triangle.py))
2. Clone the repo somewhere: `git clone https://github.com/mikejbrown/error-analysis.git`
3. cd into it: `cd error-analysis`
4. Set up the dexy build environment: `dexy setup`
5. Run dexy: `dexy`
6. If all went well you should have the fresh build document at `output/error_analysis.pdf`!


The files `accuracy-vs-precision.nb` and `uncertain-function.nb` requires Wolfram Mathematica or Mathematica
Player. It was tested on Mathematica 9.

# License
:copyright: 2014 Michael J. Brown
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License. To view a copy of this license, visit
http://creativecommons.org/licenses/by-sa/4.0/.

The human readable summary follows. This is *not* the license itself. You are free to:

* Share: copy and redistribute the material in any medium or format
* Adapt: remix, transform, and build upon the material for any purpose, even commercially.

The licensor cannot revoke these freedoms as long as you follow the license terms.

Under the following terms:

* Attribution: You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
* ShareAlike: If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.
* No additional restrictions: You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

