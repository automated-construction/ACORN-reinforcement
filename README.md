# acorn_reinforcement

Reinforcement filament layout optimisation and winding.

This package tackles automated design-to-production of reinforcement layouts using robotic filament winding.
The source code includes the key functions and the example folder suggests an end-to-end workflow to run outside or inside Rhino/Grasshopper.

**Publications**

* [Oval, R., Costa, E., Thomas-McEwen, D., Spadea, S., Orr, J. and Shepherd, P., 2020, October. Automated framework for the optimisation of spatial layouts for concrete structures reinforced with robotic filament winding. In ISARC 2020: The 37th International Symposium on Automation and Robotics in Construction (pp. 1541-1548).](http://automated.construction/publications-blog/2020/10/28/isarc-2020-37th-international-symposium-for-automation-and-robotics-in-construction)

* Automated construction of reinforcement layouts for concrete structures from constrained layout optimisation with robotic filament winding. In preparation

**Dependencies**

* [compas](https://compas.dev/) (for geometrical operations)
* [numpy](https://numpy.org/), [scipy](https://www.scipy.org/), [cvxpy](https://www.cvxpy.org/) (for numerical optimisation)
* [networkx](https://networkx.org/) (for graph data structure)

**Installation**

Via your Terminal or Prompt, go to the root of your local repo copy. Then, for intallation, run:

```
pip install -e .
```

and for access with Rhino/Grasshopper, run:

```
python -m compas_rhino.install -p compas acorn_reinforcement
```

**Contact**

Robin Oval - rpho2@cam.ac.uk
