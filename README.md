<p align="center">
  <img src="images/showcase1.gif" alt="Showcase falling into black hole" width="900">
</p>

The red sphere is the true 3D black hole position with a non-spinning Schwarzschild radius.

The grey distant sphere is the target object, other side of the black hole.

The primary and secondary images show the fastest and slowest paths for light to reach the observer (without looping).

The game engine must use a ring buffer to render in the target object `travelTime` units in the past (for their orientation etc). It must also use `timeScale` to slow inputs and physics around a relativistic player to keep players in sync (if multiplayer).

*Note: The black hole shadow rendering logic is a work in progress.*

