# Curves

- Open the file Curve/CurveFunctions.cpp where function implementation is present

- Use the file Curve.html to view the results

- the compilation was done as directed by the homework specifications using emscripten compiler using emscripten cmd prompt using below command:

emcc -Wall -Werror --bind -I. -O2 --memory-init-file 0 -o curvelib.html CurveBridge.cpp Curve.cpp CurveFunctions.cpp


- the part 2.c is giving results but is not as expected, the curve is not perfect, I have attached two snapshots for C2 continuity function(with_C2.jpeg) implemented and other one without C2 continuity(without_C2.jpeg) function. The function CalculateHermiteSplineDerivativesForC2Continuity() is commented in the attached submitted file, but can be run by removing /* */ as directed in the function itself

- I have developed all the functions first hand.

Results:

![catmullrom](https://user-images.githubusercontent.com/22354463/41386152-f8188716-6f4d-11e8-9fa0-72005b5a0d0f.PNG)

![cubichermit](https://user-images.githubusercontent.com/22354463/41386155-fb3b613e-6f4d-11e8-85e0-c1f96e21e860.PNG)

![cubicbeizerbernstein](https://user-images.githubusercontent.com/22354463/41386159-ff1b7b9a-6f4d-11e8-9325-f6d518c67da5.PNG)

