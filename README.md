# Frame-codes

This repository contains the code used to create and show the results in
Royee Yosibash's master thesis, *Irregular Polynomial Codes for Coded
Computation with Numerical Stability and Graceful Degradation*.

## MATLAB Project

The active MATLAB project is in [`Matlab/`](Matlab/). It was tested with
MATLAB 2018a.

Initialize the active MATLAB paths with:

```matlab
run('Matlab/setup.m')
```

The main workflows are:

- [`Matlab/GUI/Run.m`](Matlab/GUI/Run.m): launches the Frame Analyzer GUI.
- [`Matlab/Coding Scheme/compareCodes.m`](Matlab/Coding%20Scheme/compareCodes.m): runs coded-computation comparisons.
- [`Matlab/Results/CreateFiguresForArticle.m`](Matlab/Results/CreateFiguresForArticle.m): creates distribution and article figures.

Shared code-construction, linear-algebra, eigenvalue, distribution, and
statistics utilities are located in [`Matlab/FramesTOOLBOX/`](Matlab/FramesTOOLBOX/).
Older exploratory code and examples are kept in [`Matlab/Archive/`](Matlab/Archive/).

See [`Matlab/README.md`](Matlab/README.md) for the detailed MATLAB project
structure and workflow notes.
