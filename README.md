<!-- Improved compatibility of back to top link: See: https://github.com/dmighty007/SolOrder/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the SolOrder. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
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
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/dmighty007/SolOrder">
    <img src="docs/images/logo22.gif" alt="Logo" width="150" height="150">
  </a>

  <h3 align="center">SolOrder</h3>

  <p align="center">
    A module to calculate various structural order parameter related to water molecules. This project initiated with the aim to 
    distinguish between different polymorphs of ice.
    <br />
    <a href="https://dmighty007.github.io/SolOrder/"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://dmighty007.github.io/SolOrder/guide/first-steps/">Beginner's guide</a>
    ·
    <a href="https://github.com/dmighty007/SolOrder/issues">Report Bug</a>
    ·
    <a href="https://github.com/dmighty007/SolOrder/issues">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

<!--[![Product Name Screen Shot][product-screenshot]](https://example.com)-->

*SolOrder* is purposed to calculate simple structural property based order parameter of water molecules in MD simulation output. This module is validated for gromacs output files.   
SolOrder calculates tetrahedral order, local structure index, F4, entropy
fingerprints, density, coordination numbers, fifth-neighbor distances,
hydrogen-bond counts, Steinhardt bond order, periodic Voronoi measures, angular
distributions, centrosymmetry, and CHILL/CHILL+ classification. Calculations support configurable
Numba/joblib execution. See the [parallel API and definitions](docs/order_parameters.md)
for examples, units, neighbor conventions, and validation details. The
[implementation-status table](docs/order_parameter_status.md) lists completed
parameters and the remaining manuscript scope.
[Extended parameter examples](docs/extended_parameters.md) cover the geometric
methods and optional DScribe/PTM backends. Complex learned models are deferred.
Validation reports cover [GenIce structures](docs/validation/genice.md) and
[366 sampled nucleation configurations](docs/validation/nucleation/README.md).


<p align="right">(<a href="#readme-top">back to top</a>)</p>




<!-- GETTING STARTED -->
## Getting Started


### Prerequisites

SolOrder uses MDAnalysis, NumPy, pandas, Numba, tqdm, and joblib. Runtime
dependencies are declared in `pyproject.toml` and installed automatically by pip.
The development environment uses Python 3.11. NumPy remains below version 2
while the existing numerical implementation is validated for newer versions.

### Installation

Install the published package:

```sh
python -m pip install SolOrder
```

For development, create a conda environment and install the package in editable mode:

```sh
git clone https://github.com/dmighty007/SolOrder.git
cd SolOrder
conda env create --prefix ./.conda-env --file environment.yml
conda activate ./.conda-env
python -m pip install -e ".[dev]"
python -m pytest
```

Library code lives in `src/SolOrder/`. Research notebooks are in
`examples/notebooks/`, sample files in `examples/data/`, and saved results in
`examples/results/`. Notebooks require their original simulation inputs; update
local input paths before running them. Install plotting and notebook tools with
`python -m pip install -e ".[examples]"`.

To build release artifacts:

```sh
python -m build
```


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

```python
import MDAnalysis as mda
from SolOrder import OrderParameter, eFingerprint

u = mda.Universe("liquid.gro")                 # any format MDAnalysis reads
op = OrderParameter(u, cutoff=5.0)              # one frame; oxygens "name OW"

q = op.tetrahedral_order()                      # one value per water molecule
q6 = op.steinhardt(6, averaged=True, cutoff=3.5)
labels = op.chill_plus().names                  # "hexagonal_ice", "cubic_ice", ...
s = eFingerprint(u, cutoff=5.0).localEntropy()  # local entropy fingerprint

op.refresh(frame=10)                            # move to another trajectory frame
```

The [documentation](https://dmighty007.github.io/SolOrder/) has a
**beginner's guide** with runnable examples and small example structures
(`docs/guide/data/`): first steps, a tour of the order parameters, classifying
ice and finding a nucleus, choosing a parameter, and working with trajectories.
It also has the API reference. To read it locally:

```sh
python -m pip install -e ".[docs]"
mkdocs serve        # then open http://127.0.0.1:8000
```

<p align="center">
  <img src="docs/guide/figures/distributions.png" alt="Order-parameter distributions for ice Ih, ice Ic and liquid water" width="600">
</p>

<p align="right">(<a href="#readme-top">back to top</a>)</p>




<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Dibyendu Maity - [@MAITYDIBYENDU99](https://twitter.com/@MAITYDIBYENDU99) - dibyendumaity1999@bose.res.in

Project Link: [https://github.com/dmighty007/SolOrder](https://github.com/dmighty007/SolOrder)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Below are the references, basically from where the ideas and algorithms are taken. We might want to benchmark with their results!!
* [Investigating the quasi-liquid layer on ice surfaces: a comparison of order parameters](https://pubs.rsc.org/en/content/articlelanding/2022/CP/D2CP00752E)
* [Mining of Effective Local Order Parameters to Classify Ice Polymorphs](https://doi.org/10.1021/acs.jpca.1c06685)
* [Characterization of the Local Structure in Liquid Water by Various Order Parameters](https://doi.org/10.1021/acs.jpcb.5b02936)
* [Simulations of the methane hydrate/methane gas interface near hydrate forming conditions conditions](https://doi.org/10.1016/0378-3812(95)02903-6)

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/dmighty007/SolOrder.svg?style=for-the-badge
[contributors-url]: https://github.com/dmighty007/SolOrder/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/dmighty007/SolOrder.svg?style=for-the-badge
[forks-url]: https://github.com/dmighty007/SolOrder/network/members
[stars-shield]: https://img.shields.io/github/stars/dmighty007/SolOrder.svg?style=for-the-badge
[stars-url]: https://github.com/dmighty007/SolOrder/stargazers
[issues-shield]: https://img.shields.io/github/issues/dmighty007/SolOrder.svg?style=for-the-badge
[issues-url]: https://github.com/dmighty007/SolOrder/issues
[license-shield]: https://img.shields.io/github/license/dmighty007/SolOrder.svg?style=for-the-badge
[license-url]: https://github.com/dmighty007/SolOrder/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/dibyendu-maity-a1909b174
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 
