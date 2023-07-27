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

It is authored by [Oriol Colomés](https://oriolcolomes.com).

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

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "GridapODEsTests.jl"
```
which auto-activate the project and enable local path handling from DrWatson.

To run a test you simply need to include one of the main files in `scripts/` folder. For example, tu run the tests for manufactured analytical solutions for linear elasticity, execute
```julia
include("mainElasticity.jl")
```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ROADMAP -->
## Roadmap

- [x] Manufactured solutions for linear elasticity
- [ ] Generalize manufactured solution for linear elasticity
- [ ] Manufactured solution fol incompressible Navier-Stokes
    - [ ] Tests with schemes that don't require update
    - [ ] Tests with schemes that require update (requires handling of index 2 DAEs)
- [ ] Beltrami flow test

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this thest bench better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".

Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<!-- CONTACT -->
## Contact

Oriol Colomés - [@OriolCG](https://twitter.com/oriolcg) - j.o.colomesgene@tudelft.nl

<!-- Project Link: [https://github.com/your_username/repo_name](https://github.com/your_username/repo_name) -->


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Best-README-Template](https://github.com/othneildrew/Best-README-Template)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

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
[license-url]: https://github.com/oriolcg/GridapODEsTests.jl/blob/main/LICENSE.txt
