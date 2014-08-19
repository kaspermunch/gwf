---
layout: page
title: About
permalink: /about/
---

GWF (Grid WorkFlow) is a small utility for specifying workflows to run on a computer grid. It was developed to help me keep track on my workflows on our local computer grids and mostly driven by the needs I have in my own projects, but the goal is to make a generally useful workflow system.

GWF provides a Make-like system for specifying dependencies between tasks to be run, and it then produces scripts for submission to the grid system, with dependencies between tasks so tasks than can run in parallel will run in parallel while tasks that have to wait for other tasks to complete will be submitted but put on hold until those tasks are completed.

Dependencies are based on files, but supports tasks that takes several input files and produces several output files. Tasks then have dependencies on selected files, either produced by one or many other tasks or assumed to be present on the system.


Bug reports and feature requests
--------------------------------

To report bugs or make feature requests, please use the [github issues page](https://github.com/mailund/gwf/issues) for the project.


Contributing
------------

To contribute you should:

1. Fork the repository.
2. Make a new branch in your fork and name it feature/DESCRIPTIVE_NAME.
3. Implement your new feature.
4. Test it.
5. Open a pull request to the main repository.

