# Area average interpolation: A pixel interpolation method for suppressing loss of image information
Source code for Area average interpolation method.

# Overview
このツールの概要を一言で書く。
Area average補間法は画像情報の損失がない新しい画素補間法である。
Area average interpolation method is a novel pixel interpolation method for suppressing loss of image information.

# Features
AAのセールスポイントや差別化などを説明する。
従来の画素補間法(Bi-linear, Bi-cubic法など)では画像縮小の際に使用されない元画素が存在するため、画像情報の損失が起こる。
Area average補間法では画像情報の損失が起こらないため、より真の値に近い画素補間が可能である。

# Requirement
ツールやライブラリを使うのに依存がある場合。
c++ version : c++11
compiler : MSVC2015

# Install
ソースコードをコピペするだけでよい。

# Usage
使い方。ソースコードの何行目にユーザによる設定値がある。必要があれば編集すること。

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

