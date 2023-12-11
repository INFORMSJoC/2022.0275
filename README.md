# 2022.0275

[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE.txt).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Regret Minimization and Separation in Multi-Bidder Multi-Item Auctions](https://doi.org/10.1287/ijoc.2022.0275) by Çağıl Koçyiğit, Daniel Kuhn and Napat Rujeerapaiboon.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0275

https://doi.org/10.1287/ijoc.2022.0275.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{kocyigit2023,
  author =        {Çağıl Koçyiğit and Daniel Kuhn and Napat Rujeerapaiboon},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Regret Minimization and Separation in Multi-Bidder Multi-Item Auctions}},
  year =          {2023},
  doi =           {10.1287/ijoc.2022.0275.cd},
  url =           {https://github.com/INFORMSJoC/2022.0275},
}  
```

## Description

This repository aims to showcase the effectiveness of our proposed regret-minimizing mechanism against various benchmark mechanisms.

## Content

The [scripts](scripts) folder contains MATLAB code files for reproducing the numerical results in our paper. [CompRatioComparison.m](scripts/CompRatioComparison.m) contains the code for comparing our mechanism to benchmarks in terms of competitive ratio, while [RegretQuantiles.m](scripts/RegretQuantiles.m) contains the code for assesing regret quantiles, both for 2-bidder case. [HigherDimension.m](scripts/HigherDimension.m) includes the code for comparison with a higher number of bidders. It is separated to avoid the computational challenge of computing the nominal optimal mechanism (one of the benchmarks) in cases with more bidders. 

Note that [CompRatioComparison.m](scripts/CompRatioComparison.m) and [RegretQuantiles.m](scripts/RegretQuantiles.m) both include code to solve optimization problems to determine the nominal optimal benchmark mechanism. To execute the code as provided, it is needed to have the [YALMIP](https://yalmip.github.io/) interface and [GUROBI](https://www.gurobi.com/) solver installed, along with a valid GUROBI license.

The [results](results) folder contains the figures presented in our paper, derived from the conducted numerical experiments.
## Support

For support in using this software, submit an
[issue](https://github.com/INFORMSJoC/2022.0275/issues/new).
