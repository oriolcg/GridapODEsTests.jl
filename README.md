<a name="readme-top"></a>
<!--
*** This README file has been created using the [Best-README-Template(https://github.com/othneildrew/Best-README-Template/tree/master)
-->

<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]

# GridapODEsTests.jl

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

`GridapODEsTests.jl` is a [Julia Language](https://julialang.org/) repository that serves as a test bench for time integration schemes implemented in [Gridap.jl](https://github.com/gridap/Gridap.jl). In this repository there are tests to assess the accuracy and performance (time and memory consumption) of different ODE solvers for time-dependent PDEs using a Finite Element framework.

It is authored by [Oriol Colomes](https://oriolcolomes.com).

### Built With

This repository has been created using [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make it a reproducible scientific project

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

To run the tests in this repository you only need to install [Julia](https://julialang.org/downloads/). The tests have been executed with v1.9.2, but earlier versions should also work.

### Installation

To get a local copy up and running follow these simple example steps.

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data and outputs are not included in the git-history. They will be generated when executing the tests.
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.

<!-- USAGE EXAMPLES -->
## Usage

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "GridapODEsTests.jl"
```
which auto-activate the project and enable local path handling from DrWatson.

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/oriolcg/GridapODEsTests.jl.svg?style=for-the-badge
[contributors-url]: https://github.com/oriolcg/GridapODEsTests.jl/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/oriolcg/GridapODEsTests.jl.svg?style=for-the-badge
[forks-url]: https://github.com/oriolcg/GridapODEsTests.jl/network/members
[stars-shield]: https://img.shields.io/github/stars/oriolcg/GridapODEsTests.jl.svg?style=for-the-badge
[stars-url]: https://github.com/oriolcg/GridapODEsTests.jl/stargazers
[issues-shield]: https://img.shields.io/github/issues/oriolcg/GridapODEsTests.jl.svg?style=for-the-badge
[issues-url]: https://github.com/oriolcg/GridapODEsTests.jl/issues
[license-shield]: https://img.shields.io/github/license/oriolcg/GridapODEsTests.jl.svg?style=for-the-badge
[license-url]: https://github.com/oriolcg/GridapODEsTests.jl/blob/master/LICENSE.txt
