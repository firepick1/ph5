# ph5
C++ library for Pythagorean Hodograph Quintic curves

### Overview
Pythogarean Hodograph curves provide a clean and elegant representation for
machine traversable curves without the control point clutter that plagues cubic splines.
Quintic PH curves provide C2 continuity ideal for smooth machine acceleration.

This library is targeted for use in FirePick Delta machines, but may find use
elsewhere.

### What is a PH5Curve?
A Pythogarean Hodograph curve is a polynomial approximation to 2D curve. As such,
it is similar to a Bezier or NURBS curve. What makes a PH curve special is that
the curve passes through its control points. The on-curve control point greatly
reduces the human labor of specifying a PH curve (In reality, a PH curve does in
fact have off-curve control points, but these are secondary and entirely derived
from the on-curve control points.)

PH curves, in particular, quintic PH curves (PH5Curve) are well suited to generating
smooth curves for CNC robots. Although higher order PH curves (7,9, etc.) are also
acceptable, the PH5Curve exhibits C2 continuity for the least computational burden.
C2 continuity is important because C2 continuity asserts that the acceleration profile for a
point traversing the curve is continuous. If acceleration is discontinuous,
machinery jerks.

### Computation and Memory
Generating PH curves from a set of control points is computationally and memory 
intensive. However, traversing a PH curve is less demanding. And, finally,
traversing a piece-wise linear approximation to a a sampled PH curve is 
even less demanding.

This library is designed to run on processors as small as the ArduinoMega2560,
which only has 8KB of free RAM. On sucha a small processor, the support for
PH curves is necessarily limited--perhaps even to 3 control points.
Given more processing power and memory (e.g., BeagleBone Black) this library 
can traverse more complex PH curves having more control points.

Generally speaking, you want to generate PH curves on a fast, desktop/laptop,
but you can traverse them on more modest hardware.

### Installation: ArduinoMega2560
1. Download and install the library:

	`https://github.com/firepick1/ph5/releases`

1. Unzip and copy the **ph5** folder under the **libraries** folder of your
Arduino workspace folder

	*your-personal-Arduino-workspace-folder*`/libraries/ph5`

1. Restart your Arduino IDE

1. To use the library, include the following line in your sketch:

	`#include <ph5.h>`
