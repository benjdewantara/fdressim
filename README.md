<h1 align='center'>Single-phase Fluid Finite-difference Simulator using Python </h1>

# Contents

- [About](#about)
- [Using the simulator](#how-to-use-the-simulator)
    - [Install Python and other dependencies](#install-python-and-other-dependencies)
    - [Other specification](#other-specification)
    - [Downloading **fdressim**](#downloading-fdressim)
    - [Example script](#example-script)

# About
This project was submitted as my bachelor thesis in [Teknik Perminyakan ITB][]. It attempts to produce an implementation of finite-difference method to solve diffusivity equation for *single-phase* fluid flow in porous media in the form of [Python](https://www.python.org/) code. My cause behind this project was because many *undergraduate* students majoring in [Petroleum Engineering][] seem to fail to completely understand how a reservoir simulator works behind the monitor (how to numerically solve the diffusivity equation). Thus, I aim to show how the equation is *put into code* and how to generate the solution for each time level.

<!--Note: My pure and naive motivation was I just wanted to code ;).-->



# How to use the simulator

## Install Python and other dependencies

This simulator is written in Python 3.x

Normally, one may install Python directly using an available version of installer on <https://www.python.org/downloads/>. But since this simulator makes use of the [SciPy stack][], it is recommended to use a Python distribution instead. A [Python distribution][] is Python with more features installed in it. Features that are heavily used in this simulator are [NumPy][] and [Scipy library][]. There are many types of Python distributions such as [Anaconda] and [Canopy]. I myself coded this using Anaconda (Python 3.5 version).







## Other specification





## Downloading fdressim

There is a *Clone or download* button in fdressim Github's page. You can download the [zip file][fdressim-master-zip] and extract it.



## Example script










[Teknik Perminyakan ITB]: http://tm.itb.ac.id/
[Petroleum Engineering]: https://en.wikipedia.org/wiki/Petroleum_engineering
[Scipy stack]: http://scipy.org/
[Python distribution]: https://wiki.python.org/moin/PythonDistributions
[NumPy]: http://www.numpy.org/
[Scipy library]: http://scipy.org/scipylib/index.html
[Anaconda]: https://www.continuum.io/anaconda
[Canopy]: https://www.enthought.com/products/canopy/
[fdressim-master-zip]: https://github.com/benjdewantara/fdressim/archive/master.zip