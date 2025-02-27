# Constraint optimisation approaches for designing group-living captive breeding programmes

<div align="center">

<div>
  <a href='https://www.ncl.ac.uk/computing/staff/profile/matthewforshaw.html' target='_blank'>Matthew Forshaw</a><sup>1,2</sup>&emsp;
    Rachel Gray<sup>1</sup>&emsp;
    <a href='https://eeb.princeton.edu/people/bridgett-vonholdt' target='_blank'>Bridgett vonHoldt</a><sup>3</sup>&emsp;
    <a href='https://scholar.google.com/citations?user=qTDTHUwAAAAJ&hl=en' target='_blank'>Alexander Ochoa</a><sup>4</sup>&emsp;
    <a href='https://jmillergenomics.com/' target='_blank'>Joshua M. Miller</a><sup>5</sup>&emsp;
    <a href='https://www.mtu.edu/forest/about/faculty-staff/faculty/brzeski/' target='_blank'>Kristin E. Brzeski</a><sup>6</sup>&emsp;
    <a href='https://caccone.yale.edu/people/director' target='_blank'>Adalgisa Caccone</a><sup>4</sup>&emsp;
    <a href='https://www.ncl.ac.uk/nes/people/profile/evelynjensen.html' target='_blank'>Evelyn L. Jensen</a><sup>1</sup>&emsp;
</div>
<div>

  <sup>1</sup>Newcastle University, Newcastle upon Tyne, NE1 7RU, United Kingdom&emsp;
  <sup>2</sup>The Alan Turing Institute, London, NW1 2DB, United Kingdom&emsp;
  <sup>3</sup>Princeton University, Princeton, NJ 08544, USA&emsp;
  <sup>4</sup>Yale University, New Haven, CT 06520, USA&emsp;
  <sup>5</sup>MacEwan University, Edmonton, Canada&emsp;
  <sup>6</sup>Michigan Technological University, Michigan, USA

</div>
</div>
<p align="center">
<!--   <a href='https://arxiv.org/abs/2412.09441'><img src='https://img.shields.io/badge/Arxiv-2412.09441-b31b1b.svg?logo=arXiv'></a> -->
  <a href=""><img src="https://img.shields.io/github/stars/mattforshaw/AAAI25-CaptiveBreeding?color=4fb5ee"></a>
  <a href=""><img src="https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fmattforshaw%2FAAAI25-CaptiveBreeding&count_bg=%23FFA500&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=visitors&edge_flat=false"></a>
</p>

## Overview
Captive breeding programs play a critical role in combating the ongoing biodiversity crisis by preserving the most endangered species and supporting reintroduction efforts. Maintaining the genetic health of captive populations requires careful management to prevent inbreeding and maximize the effective population size. Decisions about which males and females should be bred together are guided by the principle of minimizing relatedness between pairs. Methods to select breeding pairs are well developed, however, some species' ecology requires them to live in groups, and evaluating optimal groupings of multiple males and females that would be suitable to breed together is a more complex problem. Current computational tools to support the design of group-living captive breeding programs suffer from challenges of scalability and flexibility. In this paper we demonstrate the applicability of constraint programming (CP) approaches to optimize breeding groups to minimize relatedness. We present the example of the Gal\'{a}pagos giant tortoises as the test case used to develop our approach. Exploration of the needs of this captive breeding program has informed the development of our flexible approach to capture the constraints on viable captive breeding program design. Our findings have directly informed the implementation of new group configurations at the captive breeding centre. We further demonstrate that our approach is broadly applicable in other contexts through a second case study, providing multi-objective optimisation of a breeding program of canids. Through these case studies and an ablation study using synthetic datasets, we show that our constraint optimisation approach provides an expressive and generalizable means to support captive breeding program design, including scaling to large captive populations, which are currently intractable using current computational methods.

## Citing this Work
🎉The code repository for "[Constraint optimisation approaches for designing group-living captive breeding programmes](TBC)" (AAAI 2025). If you use any content of this repo for your work, please cite the following bib entry:

```
  @inproceedings{sun2024mos,
    title={{Constraint optimisation approaches for designing group-living captive breeding programmes}},
    author={Forshaw, Matthew and Gray, Rachel and vonHoldt, Bridgett and Ochoa, Alexander and Miller, Joshua and Brzeski, Kristin E. and Caccone, Adalgisa and Jensen, Evelyn L.},
    booktitle={AAAI},
    year={2025}
  }
```
