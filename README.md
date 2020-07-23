# Area average interpolation: A pixel interpolation method for suppressing loss of image information
Source code for Area average interpolation method.

# Overview
Area average interpolation method is a pixel interpolation method for suppressing loss of image information.

# Features
Conventional pixel interpolation methods (Bi-linear and Bi-cubic interpolation) cause loss of image information bevause some original pixels are not used in the image reduction process.
Since the Area average interpolation method does not cause any loss of image information, it allows for pixel interpolation that is closer to the true value.

# Requirement
c++ version : c++11 and later  
compiler : MSVC2015 and later

# Install
All you have to do is to create your own project and paste the source code.

# Usage
In the 1528th ~ 1534th line of the source code, there are user settings (image path, resolution, rotation angle, etc.), so please edit them if necessary.

# Author
Name:        Isshi Nara

Affiliation: Graduate School of Biomedical Science and Engineering, Hokkaido University.

E-mail:      isshi.nara.01@gmail.com

# License
Copyright (c) 2020 Isshi Nara

All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

